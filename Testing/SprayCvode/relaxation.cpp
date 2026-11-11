#include "SprayParticles.H"
#include "SprayCvode.H"
#include "Drag.H"
#include "Transport.H"
#include <AMReX.H>
#include <AMReX_Print.H>
#include <iomanip>

// Exact Stokes relaxation tests exercise the production particle RHS and
// sources.
int
main(int argc, char** argv)
{
  amrex::Initialize(argc, argv);
  {
    int verbose = 0;
    SprayParticleContainer::readSprayParams(verbose);
    amrex::Real force[AMREX_SPACEDIM] = {};
    SprayParticleContainer::spraySetup(force);
    SprayData data = *SprayParticleContainer::getSprayData();
    data.fixed_parts = false;
    data.mass_trans = false;
    data.dtmod = 1.;
    pele::physics::PeleParams<pele::physics::transport::TransParm<
      pele::physics::EosType, pele::physics::TransportType>>
      transport;
    transport.initialize();
    GasPhaseVals ambient;
    ambient.reset();
    ambient.T_fluid = 300.;
    ambient.p_fluid = 1.e5;
    ambient.Y_fluid[N2_ID] = 1.;
    auto eos = pele::physics::PhysicsType::eos();
    eos.molecular_weight(ambient.mw.data());
    eos.Y2WBAR(ambient.Y_fluid.data(), ambient.mw_mix);
    amrex::Real rho;
    eos.PYT2R(
      ambient.p_fluid / SprayUnits::pres_conv, ambient.Y_fluid.data(),
      ambient.T_fluid, rho);
    ambient.rho_fluid = rho * SprayUnits::rho_conv;
    for (auto& w : ambient.mw) {
      w *= SprayUnits::mass_conv;
    }
    ambient.mw_mix *= SprayUnits::mass_conv;
    amrex::Real boil[SPRAY_FUEL_NUM];
    data.calcBoilT(ambient, boil);
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> fractions{};
    fractions[0] = 1.;
    const amrex::Real liquid_rho = data.liqprops.rho_mix(fractions, 300., boil);
    const amrex::Real diameter = 1.e-4 * SprayUnits::len_conv; // one micrometre
    const amrex::Real mass = M_PI / 6. * liquid_rho * std::pow(diameter, 3);
    const amrex::Real initial_velocity = 0.1 * SprayUnits::len_conv;
    amrex::RealVect velocity = amrex::RealVect::TheZeroVector();
    velocity[0] = initial_velocity;
    SprayRates rates;
    evaluateSprayRates(
      ambient, data, velocity, 300., diameter, mass, liquid_rho, fractions, 0.,
      boil, transport.device_parm(), rates);
    const amrex::Real inverse_tau = rates.drag_force_p;
    AMREX_ALWAYS_ASSERT(inverse_tau > 0. && rates.Reyn < 1.);
    pele::physics::spray::SprayCvode solver;
    for (amrex::Real stiffness : {3., 1000.}) {
      SprayParticleContainer::ParticleType p{};
      p.id() = 1;
      p.cpu() = 0;
      p.rdata(SprayComps::pstateVel) = initial_velocity;
      p.rdata(SprayComps::pstateT) = 300.;
      p.rdata(SprayComps::pstateDia) = diameter;
      p.rdata(SprayComps::pstateY) = 1.;
      p.rdata(SprayComps::pstateNumDens) = 1.;
      auto explicit_particle = p;
      GasPhaseVals gas = ambient;
      const amrex::Real dt = stiffness / inverse_tau;
      solver.advance(dt, gas, data, p, boil, transport.device_parm());
      AMREX_ALWAYS_ASSERT(p.id() > 0);
      const amrex::Real exact = initial_velocity * std::exp(-stiffness);
      const amrex::Real error = std::abs(p.rdata(0) - exact);
      const amrex::Real balance =
        std::abs(
          mass * (p.rdata(0) - initial_velocity) - dt * gas.fluid_mom_src[0]) /
        (mass * initial_velocity);
      GasPhaseVals explicit_gas = ambient;
      calculateSpraySource(
        dt, explicit_gas, data, explicit_particle, boil,
        transport.device_parm());
      const amrex::Real explicit_error =
        std::abs(explicit_particle.rdata(0) - exact);
      amrex::Print() << std::setprecision(17) << "Stokes h/tau=" << stiffness
                     << " CVODE_error=" << error
                     << " explicit_error=" << explicit_error
                     << " momentum_balance=" << balance << "\n";
      AMREX_ALWAYS_ASSERT(error < 1.e-7 * SprayUnits::len_conv);
      AMREX_ALWAYS_ASSERT(balance < 1.e-9);
      if (stiffness > 100.) {
        AMREX_ALWAYS_ASSERT(
          !std::isfinite(explicit_error) || explicit_error > 100. * error);
      }
    }
    // The moving SprayTest droplets are 1 mm at 15 m/s. In the intermediate
    // Reynolds regime, dv/dt = -a v (1 + b v^(2/3)) has a closed-form solution.
    // This checks the revised exchange scales against a nonlinear drag problem,
    // independently of a second numerical integrator.
    {
      const amrex::Real large_diameter = 0.1 * SprayUnits::len_conv;
      const amrex::Real speed = 1500. * SprayUnits::len_conv;
      const amrex::Real large_mass =
        M_PI / 6. * liquid_rho * std::pow(large_diameter, 3);
      auto gas = ambient;
      gas.T_fluid = 500.;
      eos.PYT2R(
        gas.p_fluid / SprayUnits::pres_conv, gas.Y_fluid.data(), gas.T_fluid,
        rho);
      gas.rho_fluid = rho * SprayUnits::rho_conv;
      velocity[0] = speed;
      evaluateSprayRates(
        gas, data, velocity, 300., large_diameter, large_mass, liquid_rho,
        fractions, 0., boil, transport.device_parm(), rates);
      AMREX_ALWAYS_ASSERT(rates.Reyn > 1. && rates.Reyn < 1000.);
      const amrex::Real correction = std::cbrt(rates.Reyn * rates.Reyn) / 6.;
      const amrex::Real a = rates.drag_force_p / (1. + correction);
      const amrex::Real dt = 0.005208333333;
      const amrex::Real exact =
        speed *
        std::pow(1. + (1. + correction) * std::expm1(2. * a * dt / 3.), -1.5);
      AMREX_ALWAYS_ASSERT(rates.Reyn * exact / speed > 1.);
      SprayParticleContainer::ParticleType p{};
      p.id() = 4;
      p.rdata(SprayComps::pstateVel) = speed;
      p.rdata(SprayComps::pstateT) = 300.;
      p.rdata(SprayComps::pstateDia) = large_diameter;
      p.rdata(SprayComps::pstateY) = 1.;
      p.rdata(SprayComps::pstateNumDens) = 1.;
      solver.advance(dt, gas, data, p, boil, transport.device_parm());
      const amrex::Real error =
        std::abs(p.rdata(SprayComps::pstateVel) - exact);
      const amrex::Real balance =
        std::abs(
          large_mass * (p.rdata(SprayComps::pstateVel) - speed) -
          dt * gas.fluid_mom_src[0]) /
        (large_mass * speed);
      amrex::Print() << "Nonlinear drag CVODE_error=" << error
                     << " momentum_balance=" << balance << "\n";
      AMREX_ALWAYS_ASSERT(error < 1.e-4 * SprayUnits::len_conv);
      AMREX_ALWAYS_ASSERT(balance < 1.e-9);
    }
    // Saturated gas suppresses evaporation. With the constant MP heat capacity,
    // the integrated heat exchange must equal m cp (T_new - T_old).
    data.mass_trans = true;
    auto heating = ambient;
    heating.T_fluid = 350.;
    heating.Y_fluid.fill(0.);
    heating.Y_fluid[data.dep_indx[0]] = 1.;
    eos.Y2WBAR(heating.Y_fluid.data(), heating.mw_mix);
    heating.mw_mix *= SprayUnits::mass_conv;
    SprayParticleContainer::ParticleType heated{};
    heated.id() = 3;
    heated.cpu() = 0;
    heated.rdata(SprayComps::pstateT) = 300.;
    heated.rdata(SprayComps::pstateDia) = diameter;
    heated.rdata(SprayComps::pstateY) = 1.;
    heated.rdata(SprayComps::pstateNumDens) = 1.;
    solver.advance(0.01, heating, data, heated, boil, transport.device_parm());
    const amrex::Real sensible_heat =
      mass * data.liqprops.cp_mix(fractions, 300., boil) *
      (heated.rdata(SprayComps::pstateT) - 300.);
    const amrex::Real energy_balance =
      std::abs(sensible_heat - 0.01 * heating.fluid_eng_src) / sensible_heat;
    amrex::Print() << "Heating energy_balance=" << energy_balance
                   << " final_T=" << heated.rdata(SprayComps::pstateT) << "\n";
    AMREX_ALWAYS_ASSERT(heating.fluid_mass_src == 0.);
    AMREX_ALWAYS_ASSERT(energy_balance < 1.e-9);
    AMREX_ALWAYS_ASSERT(
      std::abs(heated.rdata(SprayComps::pstateT) - 350.) < 1.e-3);

    // Complete evaporation inside a longer requested stage must transfer the
    // full initial mass exactly once, normalized by that requested interval.
    ambient.T_fluid = 1000.;
    SprayParticleContainer::ParticleType evaporating{};
    evaporating.id() = 2;
    evaporating.cpu() = 0;
    evaporating.rdata(SprayComps::pstateT) = 300.;
    evaporating.rdata(SprayComps::pstateDia) = 2. * diameter;
    evaporating.rdata(SprayComps::pstateY) = 1.;
    evaporating.rdata(SprayComps::pstateNumDens) = 1.;
    const amrex::Real evaporation_dt = 0.01;
    solver.advance(
      evaporation_dt, ambient, data, evaporating, boil,
      transport.device_parm());
    const amrex::Real mass_balance =
      std::abs(evaporation_dt * ambient.fluid_mass_src + 8. * mass) /
      (8. * mass);
    amrex::Print() << "Evaporation mass_balance=" << mass_balance
                   << " final_T=" << evaporating.rdata(SprayComps::pstateT)
                   << "\n";
    AMREX_ALWAYS_ASSERT(evaporating.id() < 0);
    AMREX_ALWAYS_ASSERT(mass_balance < 1.e-12);
    AMREX_ALWAYS_ASSERT(std::isfinite(ambient.fluid_eng_src));
    transport.deallocate();
    SprayParticleContainer::SprayCleanUp();
  }
  amrex::Finalize();
}
