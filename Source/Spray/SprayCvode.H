#ifndef PELE_SPRAY_CVODE_H
#define PELE_SPRAY_CVODE_H

#include "SprayFuelData.H"
#include <AMReX_Particle.H>
#ifndef AMREX_USE_GPU
#include <cvode/cvode.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#endif

namespace pele::physics::spray {

// One host workspace per worker, reused across independent particle advances.
class SprayCvode
{
public:
  static void readParameters();
  static bool enabled();
#ifndef AMREX_USE_GPU
  SprayCvode();
  ~SprayCvode();
  SprayCvode(const SprayCvode&) = delete;
  SprayCvode& operator=(const SprayCvode&) = delete;

  amrex::Real advance(
    amrex::Real flow_dt,
    GasPhaseVals& gas,
    const SprayData& data,
    amrex::Particle<SprayComps::pstateNum, 0>& particle,
    const amrex::Real* boil_temperature,
    pele::physics::transport::TransParm<
      pele::physics::EosType,
      pele::physics::TransportType> const* transport);

private:
  static constexpr int temperature = AMREX_SPACEDIM;
  static constexpr int masses = temperature + 1;
  static constexpr int momentum_exchange = masses + SPRAY_FUEL_NUM;
  static constexpr int energy_exchange = momentum_exchange + AMREX_SPACEDIM;
  static constexpr int equations = energy_exchange + 1;
  static int rhs(sunrealtype, N_Vector, N_Vector, void*);
  static int roots(sunrealtype, N_Vector, sunrealtype*, void*);
  void check(int flag, const char* operation) const;
  void transferResidual(amrex::Real* state, int component) const;

  SUNContext context = nullptr;
  N_Vector state = nullptr, tolerance = nullptr, constraints = nullptr;
  SUNMatrix matrix = nullptr;
  SUNLinearSolver linear_solver = nullptr;
  void* memory = nullptr;
  const GasPhaseVals* gas = nullptr;
  const SprayData* data = nullptr;
  const amrex::Real* boil_temperature = nullptr;
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* transport = nullptr;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> depletion_mass{};
  amrex::GpuArray<bool, SPRAY_FUEL_NUM> active{};
  amrex::Real distortion = 0.;
  amrex::Long particle_id = 0;
#endif
};

} // namespace pele::physics::spray
#endif
