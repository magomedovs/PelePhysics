#include "SprayCvode.H"
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <cmath>
#include <sstream>
#ifndef AMREX_USE_GPU
#include "SprayRates.H"
#include <nvector/nvector_serial.h>
#endif

namespace pele::physics::spray {
namespace {
struct Options
{
  bool enabled = false;
  amrex::Real rtol = 1.e-6;
  amrex::Real atol_temperature = 1.e-4;
  amrex::Real atol_velocity = 1.e-5 * SprayUnits::len_conv;
  amrex::Real atol_mass = 1.e-18 * SprayUnits::mass_conv;
  int max_order = 2;
  int max_steps = 10000;
} options;
} // namespace

bool
SprayCvode::enabled()
{
  return options.enabled;
}

void
SprayCvode::readParameters()
{
  options = Options{};
  amrex::ParmParse pp("particles");
  std::string method = "explicit";
  pp.query("integrator", method);
  if (method != "explicit" && method != "cvode_bdf") {
    amrex::Abort("particles.integrator must be explicit or cvode_bdf");
  }
  options.enabled = method == "cvode_bdf";
  if (!options.enabled) {
    return;
  }
#ifdef AMREX_USE_GPU
  amrex::Abort("particles.integrator=cvode_bdf currently requires a CPU build");
#else
  amrex::ParmParse cv("particles.cvode");
  cv.query("rtol", options.rtol);
  cv.query("atol_temperature", options.atol_temperature);
  cv.query("atol_velocity", options.atol_velocity);
  cv.query("atol_mass", options.atol_mass);
  cv.query("max_order", options.max_order);
  cv.query("max_steps", options.max_steps);
  if (
    !(std::isfinite(options.rtol) && options.rtol > 0.) ||
    !(std::isfinite(options.atol_temperature) &&
      options.atol_temperature > 0.) ||
    !(std::isfinite(options.atol_velocity) && options.atol_velocity > 0.) ||
    !(std::isfinite(options.atol_mass) && options.atol_mass > 0.) ||
    options.max_order < 1 || options.max_order > 5 || options.max_steps < 1) {
    amrex::Abort(
      "Invalid particles.cvode tolerances, max_order (1..5), or max_steps");
  }
  amrex::Print() << " Spray integration: CVODE BDF, maximum order "
                 << options.max_order << "\n";
#endif
}

#ifndef AMREX_USE_GPU
void
SprayCvode::check(int flag, const char* operation) const
{
  if (flag < 0) {
    amrex::Abort(
      "Spray CVODE particle " + std::to_string(particle_id) + ": " + operation +
      " failed with status " + std::to_string(flag));
  }
}

SprayCvode::SprayCvode()
{
  check(SUNContext_Create(SUN_COMM_NULL, &context), "SUNContext_Create");
  state = N_VNew_Serial(equations, context);
  if (!state) {
    amrex::Abort("Unable to allocate spray CVODE state");
  }
  tolerance = N_VClone(state);
  constraints = N_VClone(state);
  matrix = SUNDenseMatrix(equations, equations, context);
  if (!tolerance || !constraints || !matrix) {
    amrex::Abort("Unable to allocate spray CVODE workspace");
  }
  linear_solver = SUNLinSol_Dense(state, matrix, context);
  memory = CVodeCreate(CV_BDF, context);
  if (
    !state || !tolerance || !constraints || !matrix || !linear_solver ||
    !memory) {
    amrex::Abort("Unable to allocate spray CVODE workspace");
  }
  N_VConst(0., state);
  N_VConst(0., constraints);
  auto* c = N_VGetArrayPointer(constraints);
  c[temperature] = 2.;
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    c[masses + i] = 1.;
  }
  check(CVodeInit(memory, rhs, 0., state), "CVodeInit");
  check(CVodeSetUserData(memory, this), "CVodeSetUserData");
  check(
    CVodeSetLinearSolver(memory, linear_solver, matrix),
    "CVodeSetLinearSolver");
  check(CVodeSetConstraints(memory, constraints), "CVodeSetConstraints");
  check(CVodeSetMaxOrd(memory, options.max_order), "CVodeSetMaxOrd");
}

SprayCvode::~SprayCvode()
{
  CVodeFree(&memory);
  SUNLinSolFree(linear_solver);
  SUNMatDestroy(matrix);
  N_VDestroy(constraints);
  N_VDestroy(tolerance);
  N_VDestroy(state);
  SUNContext_Free(&context);
}

int
SprayCvode::rhs(sunrealtype, N_Vector y, N_Vector dy, void* user)
{
  auto& self = *static_cast<SprayCvode*>(user);
  const auto* s = N_VGetArrayPointer(y);
  auto* d = N_VGetArrayPointer(dy);
  amrex::Real mass = 0.;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> fractions{};
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    if (!std::isfinite(s[masses + i])) {
      return 1;
    }
    // Newton predictors may cross zero before CVODE checks constraints or
    // locates a depletion root. Extend the rates to that trial state using
    // zero liquid mass; never modify y or clip an accepted particle state.
    fractions[i] = amrex::max(s[masses + i], amrex::Real(0.));
    mass += fractions[i];
  }
  if (
    !(mass > 0.) || !(s[temperature] > 0.) || !std::isfinite(s[temperature])) {
    return 1;
  }
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    fractions[i] /= mass;
  }
  const amrex::Real rho = self.data->liqprops.rho_mix(
    fractions, s[temperature], self.boil_temperature);
  if (!(rho > 0.) || !std::isfinite(rho)) {
    return 1;
  }
  const amrex::Real diameter = std::cbrt(6. * mass / (M_PI * rho));
  amrex::RealVect velocity;
  for (int j = 0; j < AMREX_SPACEDIM; ++j) {
    velocity[j] = self.data->fixed_parts ? 0. : s[j];
  }
  SprayRates r;
  evaluateSprayRates(
    *self.gas, *self.data, velocity, s[temperature], diameter, mass, rho,
    fractions, self.distortion, self.boil_temperature, self.transport, r);
  for (int j = 0; j < AMREX_SPACEDIM; ++j) {
    d[j] = self.data->fixed_parts ? 0. : r.part_vel_src[j];
    d[momentum_exchange + j] =
      self.data->mom_trans
        ? r.drag_force * (self.gas->vel_fluid[j] - velocity[j]) +
            velocity[j] * r.m_dot
        : 0.;
  }
  d[temperature] = self.data->mass_trans ? r.part_temp_src : 0.;
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    d[masses + i] = self.data->mass_trans && self.active[i] ? r.mn_dot[i] : 0.;
  }
  d[energy_exchange] = self.data->mass_trans ? r.conv_src : 0.;
#ifndef PELELM_USE_SPRAY
  if (self.data->mom_trans) {
    d[energy_exchange] +=
      r.drag_force * (self.gas->vel_fluid - velocity).dotProduct(velocity);
  }
#endif
#ifndef USE_MANIFOLD_EOS
  if (self.data->mass_trans) {
    amrex::Real enthalpy[NUM_SPECIES];
    auto eos = pele::physics::PhysicsType::eos(self.data->eosparm);
    eos.T2Hi(s[temperature], enthalpy);
    for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
      d[energy_exchange] +=
        d[masses + i] * enthalpy[self.data->dep_indx[i]] * SprayUnits::eng_conv;
    }
#ifndef PELELM_USE_SPRAY
    d[energy_exchange] += 0.5 * velocity.radSquared() * r.m_dot;
#endif
  }
#else
  d[energy_exchange] = 0.;
#endif
  for (int i = 0; i < equations; ++i) {
    if (!std::isfinite(d[i])) {
      return 1;
    }
  }
  return 0;
}

int
SprayCvode::roots(sunrealtype, N_Vector y, sunrealtype* g, void* user)
{
  auto& self = *static_cast<SprayCvode*>(user);
  const auto* s = N_VGetArrayPointer(y);
  amrex::Real mass = 0.;
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    mass += s[masses + i];
    g[i + 1] = self.active[i] ? s[masses + i] - self.depletion_mass[i] : 1.;
  }
  g[0] = mass - SprayUnits::min_mass;
  return 0;
}

void
SprayCvode::transferResidual(amrex::Real* s, int component) const
{
  const amrex::Real residual = s[masses + component];
  if (residual == 0.) {
    return;
  }
  if (data->mom_trans) {
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
      s[momentum_exchange + j] -= residual * s[j];
    }
  }
#ifndef USE_MANIFOLD_EOS
  amrex::Real enthalpy[NUM_SPECIES];
  auto eos = pele::physics::PhysicsType::eos(data->eosparm);
  eos.T2Hi(s[temperature], enthalpy);
  s[energy_exchange] -=
    residual * enthalpy[data->dep_indx[component]] * SprayUnits::eng_conv;
#ifndef PELELM_USE_SPRAY
  amrex::Real ke = 0.;
  for (int j = 0; j < AMREX_SPACEDIM; ++j) {
    ke += 0.5 * s[j] * s[j];
  }
  s[energy_exchange] -= residual * ke;
#endif
#endif
  s[masses + component] = 0.;
}

amrex::Real
SprayCvode::advance(
  amrex::Real flow_dt,
  GasPhaseVals& gas_in,
  const SprayData& data_in,
  amrex::Particle<SprayComps::pstateNum, 0>& p,
  const amrex::Real* boil_in,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* transport_in)
{
  gas = &gas_in;
  data = &data_in;
  boil_temperature = boil_in;
  transport = transport_in;
  particle_id = p.id();
  const amrex::Real interval = data->dtmod * flow_dt;
  if (!(interval > 0.) || !std::isfinite(interval)) {
    amrex::Abort(
      "Spray CVODE requires a finite, positive particle time interval");
  }
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> fractions{}, initial_mass{};
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    fractions[i] = p.rdata(SprayComps::pstateY + i);
    if (!std::isfinite(fractions[i]) || fractions[i] < 0.) {
      amrex::Abort(
        "Spray CVODE particle " + std::to_string(particle_id) +
        ": invalid initial liquid composition");
    }
  }
  const amrex::Real initial_temperature = p.rdata(SprayComps::pstateT);
  const amrex::Real rho =
    data->liqprops.rho_mix(fractions, initial_temperature, boil_in);
  const amrex::Real initial_total =
    M_PI / 6. * rho * std::pow(p.rdata(SprayComps::pstateDia), 3);
  if (
    !(initial_total > 0.) || !std::isfinite(initial_total) ||
    !(initial_temperature > 0.) || !std::isfinite(initial_temperature)) {
    amrex::Abort(
      "Spray CVODE particle " + std::to_string(particle_id) +
      ": invalid initial mass, diameter, density, or temperature");
  }
  N_VConst(0., state);
  auto* s = N_VGetArrayPointer(state);
  auto* atol = N_VGetArrayPointer(tolerance);
  const amrex::Real mass_error =
    SPRAY_FUEL_NUM * options.atol_mass + options.rtol * initial_total;
  for (int j = 0; j < AMREX_SPACEDIM; ++j) {
    s[j] = data->fixed_parts ? 0. : p.rdata(SprayComps::pstateVel + j);
    atol[j] = options.atol_velocity;
    // Exchange starts at zero each stage. Scale its absolute error with the
    // physical state so its relative error test does not over-resolve a small
    // accumulated increment. Propagate both velocity and mass uncertainties.
    const amrex::Real velocity_scale =
      amrex::max(std::abs(s[j]), std::abs(gas->vel_fluid[j]));
    const amrex::Real velocity_error =
      options.atol_velocity + options.rtol * velocity_scale;
    atol[momentum_exchange + j] =
      initial_total * velocity_error + mass_error * velocity_scale;
  }
  s[temperature] = initial_temperature;
  atol[temperature] = options.atol_temperature;
  const amrex::Real cp =
    data->liqprops.cp_mix(fractions, initial_temperature, boil_in);
  if (!(cp > 0.) || !std::isfinite(cp)) {
    amrex::Abort(
      "Spray CVODE requires a finite, positive liquid heat capacity");
  }
  const amrex::Real temperature_scale =
    amrex::max(initial_temperature, gas->T_fluid);
  const amrex::Real temperature_error =
    options.atol_temperature + options.rtol * temperature_scale;
  amrex::Real energy_scale = cp * temperature_scale;
#ifndef USE_MANIFOLD_EOS
  if (data->mass_trans) {
    amrex::Real enthalpy[NUM_SPECIES];
    auto eos = pele::physics::PhysicsType::eos(data->eosparm);
    eos.T2Hi(initial_temperature, enthalpy);
    for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
      energy_scale = amrex::max(
        energy_scale,
        std::abs(enthalpy[data->dep_indx[i]]) * SprayUnits::eng_conv);
    }
  }
#endif
  atol[energy_exchange] =
    initial_total * cp * temperature_error + mass_error * energy_scale;
#ifndef PELELM_USE_SPRAY
  for (int j = 0; j < AMREX_SPACEDIM; ++j) {
    const amrex::Real velocity_scale =
      amrex::max(std::abs(s[j]), std::abs(gas->vel_fluid[j]));
    atol[energy_exchange] +=
      initial_total * velocity_scale *
        (options.atol_velocity + options.rtol * velocity_scale) +
      0.5 * mass_error * velocity_scale * velocity_scale;
  }
#endif
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    initial_mass[i] = fractions[i] * initial_total;
    s[masses + i] = initial_mass[i];
    atol[masses + i] = options.atol_mass;
    depletion_mass[i] = amrex::min(options.atol_mass, initial_mass[i] * 1.e-6);
    active[i] = initial_mass[i] > 0.;
  }
  distortion = data->do_breakup == 1 ? p.rdata(SprayComps::pstateBM1) : 0.;
  check(CVodeReInit(memory, 0., state), "CVodeReInit");
  check(
    CVodeSVtolerances(memory, options.rtol, tolerance), "CVodeSVtolerances");
  amrex::Real time = 0.;
  bool removed = data->mass_trans && initial_total <= SprayUnits::min_mass;
  long int steps = 0;
  while (!removed && time < interval) {
    check(CVodeRootInit(memory, SPRAY_FUEL_NUM + 1, roots), "CVodeRootInit");
    check(CVodeSetStopTime(memory, interval), "CVodeSetStopTime");
    check(
      CVodeSetMaxNumSteps(memory, options.max_steps - steps),
      "CVodeSetMaxNumSteps");
    const int flag = CVode(memory, interval, state, &time, CV_NORMAL);
    if (flag < 0) {
      std::ostringstream detail;
      detail << "advance over h=" << interval
             << ", initial T=" << initial_temperature
             << ", initial m=" << initial_total << ", returned t=" << time;
      for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
        if (s[masses + i] < 0.) {
          detail << ", negative component " << i << "=" << s[masses + i];
        }
      }
      check(flag, detail.str().c_str());
    }
    long int segment_steps = 0;
    check(CVodeGetNumSteps(memory, &segment_steps), "CVodeGetNumSteps");
    steps += segment_steps;
    if (flag != CV_ROOT_RETURN) {
      break;
    }
    int found[SPRAY_FUEL_NUM + 1];
    check(CVodeGetRootInfo(memory, found), "CVodeGetRootInfo");
    removed = found[0] != 0;
    for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
      // Root interpolation can put another nearly depleted component slightly
      // below zero even though accepted steps satisfy CVODE's constraints.
      // Account for simultaneous depletion, including tolerance-sized roundoff,
      // before restarting the solver from the event state.
      if (s[masses + i] < -10. * options.atol_mass) {
        check(CV_CONSTR_FAIL, "negative component mass at depletion root");
      }
      if (removed || found[i + 1] || s[masses + i] <= depletion_mass[i]) {
        transferResidual(s, i);
        active[i] = false;
      }
    }
    if (!removed) {
      if (steps >= options.max_steps) {
        check(CV_TOO_MUCH_WORK, "step budget");
      }
      check(CVodeReInit(memory, time, state), "CVodeReInit after depletion");
    }
  }
  if (removed) {
    for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
      transferResidual(s, i);
    }
  }
  amrex::Real total = 0.;
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    total += s[masses + i];
  }
  if (!removed && !(total > 0.)) {
    check(CV_RHSFUNC_FAIL, "nonpositive final mass");
  }
  const amrex::Real number = p.rdata(SprayComps::pstateNumDens);
  gas_in.fluid_mass_src = number * (total - initial_total) / interval;
  for (int pc = 0; pc < data->N_pc; ++pc) {
    const int spec = data->pc_indx[pc];
    amrex::Real delta = 0.;
    for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
      if (data->dep_indx[i] == spec) {
        delta += s[masses + i] - initial_mass[i];
      }
    }
    gas_in.fluid_Y_dot[spec] = number * delta / interval;
  }
  for (int j = 0; j < AMREX_SPACEDIM; ++j) {
    gas_in.fluid_mom_src[j] = number * s[momentum_exchange + j] / interval;
    p.rdata(SprayComps::pstateVel + j) = s[j];
  }
  gas_in.fluid_eng_src = number * s[energy_exchange] / interval;
  p.rdata(SprayComps::pstateT) = s[temperature];
  if (removed) {
    p.id() = -1;
    return 0.;
  }
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    fractions[i] = s[masses + i] / total;
    p.rdata(SprayComps::pstateY + i) = fractions[i];
  }
  const amrex::Real final_rho =
    data->liqprops.rho_mix(fractions, s[temperature], boil_in);
  const amrex::Real diameter = std::cbrt(6. * total / (M_PI * final_rho));
  p.rdata(SprayComps::pstateDia) = diameter;
  amrex::RealVect velocity;
  for (int j = 0; j < AMREX_SPACEDIM; ++j) {
    velocity[j] = s[j];
  }
  SprayRates rates;
  evaluateSprayRates(
    *gas, *data, velocity, s[temperature], diameter, total, final_rho,
    fractions, distortion, boil_in, transport, rates);
  return rates.Reyn;
}
#endif
} // namespace pele::physics::spray
