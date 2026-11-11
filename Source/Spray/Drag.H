#ifndef DRAG_H
#define DRAG_H

#include "SprayRates.H"

// Compute source terms and update particles
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calculateSpraySource(
  const amrex::Real flow_dt,
  GasPhaseVals& gpv,
  const SprayData& fdat,
  SprayParticleContainer::ParticleType& p,
  amrex::Real* cBoilT,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* trans_parm)
{
  auto eos = pele::physics::PhysicsType::eos(fdat.eosparm);
  const amrex::Real C_eps = 1.E-15;
  const amrex::Real min_mass = SprayUnits::min_mass;
  const int nSubMax = 100;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_part{}, oldY_part{};
  SprayRates rates;
  const amrex::Real num_ppp = p.rdata(SprayComps::pstateNumDens);
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SprayComps::pstateVel), p.rdata(SprayComps::pstateVel + 1),
    p.rdata(SprayComps::pstateVel + 2)));
  // TAB model distortion
  amrex::Real y_tab = 0.;
  if (fdat.do_breakup == 1) {
    y_tab = p.rdata(SprayComps::pstateBM1);
  }
  // If particle is fixed in place, make velocity zero
  if (fdat.fixed_parts) {
    vel_part = amrex::RealVect::TheZeroVector();
  }
  amrex::Real T_part = p.rdata(SprayComps::pstateT);
  amrex::Real dia_part = p.rdata(SprayComps::pstateDia);

  // Storing Y_{d,n}
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_part[spf] = p.rdata(SprayComps::pstateY + spf);
  }
  // Calculating \rho_{d}
  amrex::Real rho_part = fdat.liqprops.rho_mix(Y_part, T_part, cBoilT);

  amrex::Real dt = flow_dt;
  int isub = 1;
  int nsub = 1;
  amrex::Real pmass = M_PI / 6. * rho_part * std::pow(dia_part, 3);
  amrex::Real startmass = pmass;
  amrex::Real Reyn;
  amrex::RealVect part_vel_src;
  while (isub <= nsub) {
    evaluateSprayRates(
      gpv, fdat, vel_part, T_part, dia_part, pmass, rho_part, Y_part, y_tab,
      cBoilT, trans_parm, rates);
    Reyn = rates.Reyn;
    part_vel_src = rates.part_vel_src;
    const auto& mn_dot = rates.mn_dot;
    const amrex::Real part_temp_src = rates.part_temp_src;
    const amrex::Real m_dot = rates.m_dot;
    if (fdat.mom_trans) {
      gpv.fluid_mom_src +=
        num_ppp * rates.drag_force * (gpv.vel_fluid - vel_part);
#ifndef PELELM_USE_SPRAY
      amrex::Real S_dmu_dot_u = (gpv.vel_fluid - vel_part).dotProduct(vel_part);
      gpv.fluid_eng_src += num_ppp * rates.drag_force * S_dmu_dot_u;
#endif
      if (isub == 1) {
        amrex::Real inv_tau_var = rates.drag_force_p;
        nsub = amrex::min<int>(
          amrex::max<int>(nsub, static_cast<int>(flow_dt * inv_tau_var) + 1),
          nSubMax);
      }
    }
    if (fdat.mass_trans) {
      gpv.fluid_eng_src += num_ppp * rates.conv_src;
      if (isub == 1 && rates.delT > C_eps) {
        amrex::Real inv_tau_d = -m_dot / (0.2 * pmass);
        amrex::Real inv_tau_T =
          rates.conv_src * (1. / pmass) / (rates.cp_part * rates.delT);
        nsub = amrex::min<int>(
          amrex::max<int>(
            nsub, amrex::max<int>(
                    static_cast<int>(flow_dt * inv_tau_d),
                    static_cast<int>(flow_dt * inv_tau_T)) +
                    1),
          nSubMax);
      }
    }

    if (isub == 1) {
      dt = flow_dt / static_cast<amrex::Real>(nsub);
    }
    const amrex::Real part_dt = fdat.dtmod * dt;
    if (!fdat.fixed_parts) {
      // Update particle velocity
      AMREX_D_TERM(
        vel_part[0] += part_dt * part_vel_src[0];
        , vel_part[1] += part_dt * part_vel_src[1];
        , vel_part[2] += part_dt * part_vel_src[2];);
    }
    if (fdat.mass_trans) {
      amrex::Real newT = T_part + part_dt * part_temp_src;
      amrex::Real new_mass = pmass + m_dot * part_dt;
      // If droplet is still reasonable size
      if (new_mass > min_mass) {
        T_part = newT;
#if SPRAY_FUEL_NUM > 1
        amrex::Real sumY = 0.;
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          Y_part[spf] = amrex::min<amrex::Real>(
            1.,
            amrex::max<amrex::Real>(
              0., (Y_part[spf] * pmass + mn_dot[spf] * part_dt) / new_mass));
          oldY_part[spf] = p.rdata(SprayComps::pstateY + spf);
          if (Y_part[spf] <= 1.E-12) {
            Y_part[spf] = 0.;
          }
          sumY += Y_part[spf];
        }
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          Y_part[spf] /= sumY;
        }
        rho_part = fdat.liqprops.rho_mix(Y_part, T_part, cBoilT);
#else
        oldY_part[0] = p.rdata(SprayComps::pstateY);
        rho_part =
          fdat.liqprops.rho_i(amrex::min<amrex::Real>(T_part, cBoilT[0]), 0);
#endif
        pmass = new_mass;
        dia_part = std::cbrt(6. * pmass / (M_PI * rho_part));
      } else {
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          oldY_part[spf] = p.rdata(SprayComps::pstateY + spf);
        }
        pmass = 0.;
        p.id() = -1;
        nsub = isub;
        isub = nsub + 1;
      }
    }
    ++isub;
  }
  if (nsub > 1) {
    gpv.fluid_eng_src /= static_cast<amrex::Real>(nsub);
    gpv.fluid_mom_src /= static_cast<amrex::Real>(nsub);
  }
  // Must add any mass related sources at the end in case
  // some species disappear completely
  amrex::Real mdot_total = (pmass - startmass) / (fdat.dtmod * flow_dt);
  gpv.fluid_mass_src = num_ppp * mdot_total;
#ifndef PELELM_USE_SPRAY
  amrex::Real part_ke = 0.5 * vel_part.radSquared();
  gpv.fluid_eng_src += num_ppp * part_ke * mdot_total;
#endif
  if (fdat.mom_trans) {
    gpv.fluid_mom_src += num_ppp * vel_part * mdot_total;
  }

#ifndef USE_MANIFOLD_EOS
  amrex::Real h_gas[NUM_SPECIES];
  eos.T2Hi(T_part, h_gas);
#endif
  for (int pc = 0; pc < fdat.N_pc; ++pc) {
    const int pcspec = fdat.pc_indx[pc];
    amrex::Real oldY_part_pc =
      fdat.binary_csr_rowdotprod(pcspec, fdat.L_row, fdat.L_col, oldY_part);
    amrex::Real newY =
      fdat.binary_csr_rowdotprod(pcspec, fdat.L_row, fdat.L_col, Y_part);
    amrex::Real midot =
      (newY * pmass - oldY_part_pc * startmass) / (fdat.dtmod * flow_dt);
    gpv.fluid_Y_dot[pcspec] = num_ppp * midot;
#ifdef USE_MANIFOLD_EOS
    gpv.fluid_eng_src = 0.0;
#else
    gpv.fluid_eng_src += num_ppp * midot * h_gas[pcspec] * SprayUnits::eng_conv;
#endif
  }
  AMREX_D_TERM(p.rdata(SprayComps::pstateVel) = vel_part[0];
               , p.rdata(SprayComps::pstateVel + 1) = vel_part[1];
               , p.rdata(SprayComps::pstateVel + 2) = vel_part[2];);
  p.rdata(SprayComps::pstateT) = T_part;
  p.rdata(SprayComps::pstateDia) = dia_part;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    p.rdata(SprayComps::pstateY + spf) = Y_part[spf];
  }
  return Reyn;
}

#endif
