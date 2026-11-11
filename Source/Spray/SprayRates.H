#ifndef PELE_SPRAY_RATES_H
#define PELE_SPRAY_RATES_H

#include "SprayFuelData.H"
#include "Transport.H"

// Compute the heat transfer coefficient using the
// corrected Nusselt number and B_T value
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calcHeatCoeff(
  const amrex::Real& ratio,
  const amrex::Real& B_M,
  const amrex::Real& B_eps,
  const amrex::Real& C_eps,
  const amrex::Real& Nu_0)
{
  if (B_M <= C_eps) {
    return 0.;
  }
  const int maxIter = 100;
  const amrex::Real NU2 = Nu_0 - 2.;
  const amrex::Real BM1 = 1. + B_M;
  amrex::Real phi = ratio / Nu_0;
  amrex::Real B_T_old = std::pow(BM1, phi) - 1.;
  amrex::Real logB = std::log1p(B_T_old);
  amrex::Real invFT = B_T_old / (logB * std::pow(1. + B_T_old, 0.7));
  amrex::Real Nu_num = 2. + NU2 * invFT;
  phi = ratio / Nu_num;
  amrex::Real B_T = std::pow(BM1, phi) - 1.;
  amrex::Real error = std::abs(B_T - B_T_old);
  int k = 0;
  while (k < maxIter && error > B_eps) {
    B_T_old = B_T;
    logB = std::log1p(B_T);
    invFT = B_T / (logB * std::pow(1. + B_T, 0.7));
    Nu_num = 2. + NU2 * invFT;
    phi = ratio / Nu_num;
    B_T = std::pow(BM1, phi) - 1.;
    error = std::abs(B_T - B_T_old);
    k++;
  }
  logB = std::log1p(B_T);
  invFT = B_T / (logB * std::pow(1. + B_T, 0.7));
  Nu_num = 2. + NU2 * invFT;
  return Nu_num * logB / B_T;
}

// Compute the state in the vapor and skin phase
//
template <typename EOSType>
struct calcVaporState
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    const SprayData& fdat,
    const GasPhaseVals& gpv,
    const amrex::Real& rule,
    const amrex::Real& T_in,
    const amrex::Real& C_eps,
    const amrex::Real& mw_part,
    const amrex::Real* Y_dn,
    const amrex::Real T_skin,
    const amrex::Real* cBoilT,
    amrex::Real* Y_skin,
    amrex::Real* X_vn,
    amrex::Real* L_n,
    amrex::Real& B_M,
    amrex::Real& sumXVap,
    amrex::Real& cp_skin,
    amrex::Real& mw_skin)
  {
    auto eos = pele::physics::PhysicsType::eos(fdat.eosparm);
    amrex::Real h_gas[NUM_SPECIES];
    eos.T2Hi(T_in, h_gas);

    // Mole fraction of the droplet, need vector for many-to-one mapping later
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> X_dn{};
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const amrex::Real mw_ln = fdat.liqprops.mw_i(spf);
      X_dn[spf] = Y_dn[spf] * mw_part / mw_ln;
    }

    // Estimate X_{v,n} from the liquid state
    amrex::Real mbar_vap = 0.; // Average molar mass of evaporated fuel in vapor
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const int fdspec = fdat.dep_indx[spf];

      // Compute latent heat
      const amrex::Real T_part =
        amrex::min<amrex::Real>(T_in, 0.99 * cBoilT[spf]);
      const amrex::Real part_latent = fdat.liqprops.latent_i(
        T_part, h_gas[fdspec] * SprayUnits::eng_conv, spf);
      L_n[spf] = part_latent;
      const amrex::Real pres_sat =
        fdat.liqprops.psat_i(T_part, spf, part_latent);

      // Mole fraction in gas phase
      const amrex::Real X_gi =
        gpv.Y_fluid[fdspec] * gpv.mw_mix / gpv.mw[fdspec];

      // Estimate X_gn by distributing X_gi across all the liquid species
      const amrex::Real sumXdn_SLi =
        fdat.binary_csr_rowdotprod(fdspec, fdat.L_row, fdat.L_col, X_dn);
      amrex::Real X_gn = X_dn[spf] * X_gi;
      if (sumXdn_SLi > X_dn[spf]) {
        // Multiple species deposit to fdspec, scale accordingly
        X_gn /= sumXdn_SLi;
      }

      if (pres_sat > X_gn * gpv.p_fluid) {
        const amrex::Real mw_ln = fdat.liqprops.mw_i(spf);
        X_vn[spf] = X_dn[spf] * pres_sat / gpv.p_fluid;
        mbar_vap += X_vn[spf] * mw_ln;
        sumXVap += X_vn[spf];
      } else {
        X_vn[spf] = 0.;
      }
    }

    // Collapse vapor state X_{v,n} on to dependent gas species X_{v,i}
    amrex::GpuArray<amrex::Real, NUM_SPECIES> X_vi{};
    for (int pc = 0; pc < fdat.N_pc; ++pc) {
      const int pcspec = fdat.pc_indx[pc];
      X_vi[pcspec] =
        fdat.binary_csr_rowdotprod(pcspec, fdat.L_row, fdat.L_col, X_vn);
    }

    // Mass fractions in vapor state
    const amrex::Real totalmwx = mbar_vap + (1. - sumXVap) * gpv.mw_mix;
    amrex::Real sumYSkin = 0.;   // Sum of mass frac of fuel in the skin phase
    amrex::Real sumYfFluid = 0.; // Sum of mass frac of fuel in the gas phase
    amrex::Real sumYVap = 0.;    // Sum of mass frac of fuel in the vapor phase
    amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_vi{};
    for (int pc = 0; pc < fdat.N_pc; ++pc) {
      const int pcspec = fdat.pc_indx[pc];
      if (X_vi[pcspec] > 0.) {
        Y_vi[pcspec] = X_vi[pcspec] * gpv.mw[pcspec] / totalmwx;
        sumYVap += Y_vi[pcspec];
      }
    }

    for (int pc = 0; pc < fdat.N_pc; ++pc) {
      const int pcspec = fdat.pc_indx[pc];
      if (X_vi[pcspec] > 0.) {
        const amrex::Real Y_gi = gpv.Y_fluid[pcspec];
        sumYfFluid += Y_gi;
        const amrex::Real Y_ri = Y_vi[pcspec] + rule * (Y_gi - Y_vi[pcspec]);
        Y_skin[pcspec] = Y_ri;
        sumYSkin += Y_ri;
      }
    }

    // Normalize skin mass fractions to ensure they sum to 1
    amrex::Real renorm = (1. - sumYSkin) / (1. - sumYfFluid);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      if (Y_skin[ns] == 0.) {
        Y_skin[ns] = gpv.Y_fluid[ns] * renorm;
      }
    }
    // eos and unit conversions
    eos.Y2WBAR(Y_skin, mw_skin);
    eos.TY2Cp(T_skin, Y_skin, cp_skin);
    mw_skin *= SprayUnits::mass_conv;
    cp_skin *= SprayUnits::eng_conv;
    B_M =
      (sumYVap - sumYfFluid) / amrex::max<amrex::Real>(C_eps, (1. - sumYVap));
    B_M = amrex::min<amrex::Real>(20., B_M);
  }
};

template <>
struct calcVaporState<pele::physics::eos::Manifold>
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    const SprayData& fdat,
    const GasPhaseVals& gpv,
    const amrex::Real& rule,
    const amrex::Real& T_in,
    const amrex::Real& C_eps,
    const amrex::Real& mw_part,
    const amrex::Real* Y_dn,
    const amrex::Real T_skin,
    const amrex::Real* cBoilT,
    amrex::Real* Y_skin,
    amrex::Real* X_vn,
    amrex::Real* L_n,
    amrex::Real& B_M,
    amrex::Real& sumXVap,
    amrex::Real& cp_skin,
    amrex::Real& mw_skin)
  {
    amrex::Real mbar_vap = 0.; // mean molar mass of evaporated fuel in vapor
    auto eos = pele::physics::PhysicsType::eos(fdat.eosparm);

    // Mole fraction of the droplet, need vector for many-to-one mapping later
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> X_dn{};
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const amrex::Real mw_ln = fdat.liqprops.mw_i(spf);
      X_dn[spf] = Y_dn[spf] * mw_part / mw_ln;
    }

    // Estimate X_{v,n} from the liquid state
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const int fdspec = fdat.dep_indx[spf];

      // Compute latent heat using correlation
      amrex::Real T_part = amrex::min<amrex::Real>(T_in, 0.99 * cBoilT[spf]);
      amrex::Real part_latent =
        fdat.liqprops.latent_correlation_i(fdat.liqprops.boilT_i(spf), spf);
      L_n[spf] = part_latent;

      // psat - any approach should work because we populated part_latent
      amrex::Real pres_sat = fdat.liqprops.psat_i(T_part, spf, part_latent);

      // Mole fraction in gas phase
      const amrex::Real X_gi =
        gpv.Yspec_fluid[fdspec] * gpv.mw_mix / gpv.mw[fdspec];

      // Estimate X_gn by distributing X_gi across all the liquid species
      const amrex::Real sumXdn_SLi =
        fdat.binary_csr_rowdotprod(fdspec, fdat.L_row, fdat.L_col, X_dn);
      amrex::Real X_gn = X_dn[spf] * X_gi;
      if (sumXdn_SLi > X_dn[spf]) {
        // Multiple species deposit to fdspec, scale accordingly
        X_gn /= sumXdn_SLi;
      }

      if (pres_sat > X_gn * gpv.p_fluid) {
        const amrex::Real mw_ln = fdat.liqprops.mw_i(spf);
        X_vn[spf] = X_dn[spf] * pres_sat / gpv.p_fluid;
        mbar_vap += X_vn[spf] * mw_ln;
        sumXVap += X_vn[spf];
      } else {
        X_vn[spf] = 0.;
      }
    }

    // Collapse vapor state X_{v,n} on to dependent gas species X_{v,i}
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> X_vi{};
    for (int pc = 0; pc < fdat.N_pc; ++pc) {
      const int pcspec = fdat.pc_indx[pc];
      X_vi[pcspec] =
        fdat.binary_csr_rowdotprod(pcspec, fdat.L_row, fdat.L_col, X_vn);
    }

    // Mass fractions in vapor state
    // Y_{g,i}^{pc} is already stored in gpv.Yspec_fluid
    // Calculate Y_{v,i}, Y_{r,i}^{pc}, and sum Y_{r,i}^{pc}
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_vi{0.0};
    amrex::Real sumYVap = 0.0;
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_ri_pc{0.0};
    amrex::Real sumYripc = 0.0;
    amrex::Real sumYgipc = 0.0;
    amrex::Real totalmwx = mbar_vap + (1. - sumXVap) * gpv.mw_mix;
    for (int pc = 0; pc < fdat.N_pc; ++pc) {
      const int pcspec = fdat.pc_indx[pc];
      if (X_vi[pcspec] > 0.) {
        Y_vi[pcspec] = X_vi[pcspec] * gpv.mw[pcspec] / totalmwx;
        sumYVap += Y_vi[pcspec];
        Y_ri_pc[pcspec] =
          Y_vi[pcspec] + rule * (gpv.Yspec_fluid[pcspec] - Y_vi[pcspec]);
        sumYripc += Y_ri_pc[pcspec];         //\sigma Y_{r,k}^{pc} in the docs
        sumYgipc += gpv.Yspec_fluid[pcspec]; //\sigma Y_{g,k}^{pc} in the docs
      }
    }

    // Now calculate theta
    amrex::Real theta = (1.0 - sumYripc) / (1.0 - sumYgipc);

    // Now calculate \zi_{r,j}
    amrex::GpuArray<amrex::Real, NUM_SPECIES - 1> WYgi{
      0.0}; // W_ij Y_{g,i}^{pc}
    amrex::GpuArray<amrex::Real, NUM_SPECIES - 1> WYri{
      0.0}; // W_ij Y_{r,i}^{pc}

    // Alternative way of doing above without using W_ij assuming W_ij is some
    // row exchanged form of the identity matrix
    for (int pc = 0; pc < fdat.N_pc; ++pc) {
      const int pcspec = fdat.pc_indx[pc];
      if (X_vi[pcspec] > 0.) {
        const int pcmanspec = fdat.pc_manifold_indx[pc];
        WYgi[pcmanspec] = gpv.Yspec_fluid[pcspec];
        WYri[pcmanspec] = Y_ri_pc[pcspec];
      }
    }

    for (int n = 0; n < NUM_SPECIES - 1; n++) {
      const amrex::Real zi_rj = WYri[n] + theta * (gpv.Y_fluid[n] - WYgi[n]);
      Y_skin[n] = zi_rj;
    }

    eos.Y2WBAR(Y_skin, mw_skin);
    mw_skin *= SprayUnits::mass_conv;

    // T_skin is just a dummy here
    eos.TY2Cp(T_skin, Y_skin, cp_skin);
    cp_skin *= SprayUnits::eng_conv;

    B_M = (sumYVap - sumYgipc) / amrex::max<amrex::Real>(C_eps, (1. - sumYVap));
    B_M = amrex::min<amrex::Real>(20., B_M);
  }
};

struct SprayRates
{
  amrex::RealVect part_vel_src;
  amrex::Real part_temp_src = 0.;
  amrex::Real m_dot = 0.;
  amrex::Real drag_force = 0.;
  amrex::Real drag_force_p = 0.;
  amrex::Real conv_src = 0.;
  amrex::Real cp_part = 0.;
  amrex::Real delT = 0.;
  amrex::Real Reyn = 0.;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> mn_dot{};
};

// Evaluate instantaneous rates only. The caller owns time integration and
// sources. mn_dot is caller-owned to preserve the legacy explicit substep
// behavior.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
evaluateSprayRates(
  const GasPhaseVals& gpv,
  const SprayData& fdat,
  const amrex::RealVect& vel_part,
  const amrex::Real T_part,
  const amrex::Real dia_part,
  const amrex::Real pmass,
  const amrex::Real rho_part,
  const amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM>& Y_part,
  const amrex::Real y_tab,
  const amrex::Real* cBoilT,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* trans_parm,
  SprayRates& rates)
{
  auto eos = pele::physics::PhysicsType::eos(fdat.eosparm);
  const amrex::Real rule = 1. / 3.;
  const amrex::Real C_eps = 1.E-15;
  const amrex::Real B_eps = 1.E-7;
  bool get_xi = false, get_Ddiag = fdat.mass_trans;
  bool get_lambda = fdat.mass_trans, get_mu = true, get_chi = false;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_skin{};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> rhoD_ri{};
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> L_n{}, X_vn{};
  auto& mn_dot = rates.mn_dot;
  const amrex::Real cp_part = fdat.liqprops.cp_mix(Y_part, T_part, cBoilT);
  // Average molar mass of liquid droplet
  const amrex::Real mw_part = fdat.liqprops.mw_mix(Y_part);
  // Model the fuel vapor using the one-third rule
  amrex::Real delT = amrex::max<amrex::Real>(gpv.T_fluid - T_part, 0.);
  amrex::Real T_skin = gpv.T_fluid;
  if (fdat.mass_trans) {
    T_skin = T_part + rule * delT;
  }
  Y_skin.fill(0.0);
  amrex::Real cp_skin = 0.; // Average C_p in modeled skin phase
  amrex::Real mw_skin = 0.; // Average molar mass of skin phase
  amrex::Real B_M = 0.;     // Mass Spalding number
  amrex::Real sumXVap = 0.; // Sum of X_vn

  if (fdat.mass_trans) {
    calcVaporState<pele::physics::PhysicsType::eos_type>()(
      fdat, gpv, rule, T_part, C_eps, mw_part, Y_part.data(), T_skin, cBoilT,
      Y_skin.data(), X_vn.data(), L_n.data(), B_M, sumXVap, cp_skin, mw_skin);
  } else {
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_skin[n] = gpv.Y_fluid[n];
    }
    eos.TY2Cp(T_skin, Y_skin.data(), cp_skin);
    cp_skin *= SprayUnits::eng_conv;
    mw_skin = gpv.mw_mix;
  }
  amrex::Real lambda_skin = 0.;
  amrex::Real mu_skin = 0.;
  amrex::Real xi_skin = 0.;
  amrex::Real rho_skin, rho_cgs;
  if (fdat.mass_trans) {
    eos.PYT2R(
      gpv.p_fluid / SprayUnits::pres_conv, Y_skin.data(), T_skin, rho_cgs);
    rho_skin = rho_cgs * SprayUnits::rho_conv;
  } else {
    rho_skin = gpv.rho_fluid;
    rho_cgs = rho_skin / SprayUnits::rho_conv;
  }
  auto trans = pele::physics::PhysicsType::transport();

  trans.transport(
    get_xi, get_mu, get_lambda, get_Ddiag, get_chi, T_skin, rho_cgs,
    Y_skin.data(), rhoD_ri.data(), nullptr, mu_skin, xi_skin, lambda_skin,
    trans_parm);

  mu_skin *= SprayUnits::mu_conv;
  lambda_skin *= SprayUnits::lambda_conv;
#ifdef USE_MANIFOLD_EOS
  // Hack: manifold model may inaccurately predict temperature in skin layer
  // to improve predictions, adjust computed transport coefficients to the
  // skin layer temperature using a simple temperature scaling based on
  // Sutherland's law. Rho also needs to be modified (uses ideal gas law)
  amrex::Real T_skin_predicted = 300.0;
  eos.HY2T(0.0, Y_skin.data(), T_skin_predicted);
  const amrex::Real temp_ratio = T_skin / T_skin_predicted;
  const amrex::Real sutherland_temp = 110.4;
  const amrex::Real sutherland_ratio = std::sqrt(temp_ratio) * temp_ratio *
                                       (T_skin_predicted + sutherland_temp) /
                                       (T_skin + sutherland_temp);
  mu_skin *= sutherland_ratio;
  lambda_skin *= sutherland_ratio;
  for (int pc = 0; pc < fdat.N_pc; ++pc) {
    const int pcspec = fdat.pc_indx[pc];
    rhoD_ri[pcspec] *= sutherland_ratio;
  }
  rho_skin /= temp_ratio;
#endif
  amrex::RealVect diff_vel = gpv.vel_fluid - vel_part;
  amrex::Real diff_vel_mag = diff_vel.vectorLength();
  // Local Reynolds number
  amrex::Real Reyn = rho_skin * diff_vel_mag * dia_part / mu_skin;

  // Solve momentum source terms
  const amrex::Real inv_pmass = 1. / pmass;

  amrex::Real drag_coef = 0.;
  if (Reyn > 1000.) {
    drag_coef = 0.424;
  } else if (Reyn > 1.) {
    drag_coef = 24. / Reyn * (1. + std::cbrt(Reyn * Reyn) / 6.);
  } else if (Reyn > 0.) {
    drag_coef = 24. / Reyn;
  }
  drag_coef *= (1. + 2.6432 * y_tab);
  amrex::Real drag_force =
    0.125 * rho_skin * drag_coef * M_PI * dia_part * dia_part * diff_vel_mag;
  amrex::Real drag_force_p =
    0.75 * rho_skin * drag_coef * diff_vel_mag / (dia_part * rho_part);
  amrex::RealVect part_vel_src = drag_force_p * diff_vel + fdat.body_force;

  // Solve mass and energy transfer source terms
  amrex::Real part_temp_src = 0.;
  amrex::Real conv_src = 0.;
  amrex::Real m_dot = 0.;
  if (fdat.mass_trans) {
    amrex::Real Pr_skin = mu_skin * cp_skin / lambda_skin;

    amrex::Real powR = amrex::max<amrex::Real>(
      std::pow(amrex::min<amrex::Real>(400., Reyn), 0.077), 1.);
    amrex::Real Nu_num = 0.;
    amrex::Real sumL = 0.;
    amrex::Real Nu_0 = 1. + powR * std::cbrt(1. + Reyn * Pr_skin);
    //  If gas phase is not saturated
    if (sumXVap > 0.) {
      amrex::Real rhoDtotal = 0.;
      amrex::Real cp_fuel = 0.;
      amrex::Real sumYSkin = 0.;
      for (int pc = 0; pc < fdat.N_pc; ++pc) {
        const int pcspec = fdat.pc_indx[pc];
        // Convert mass diffusion coefficient from mixture average
        // to binary for fuel only, not concerned with other species
        rhoD_ri[pcspec] *= mw_skin / gpv.mw[pcspec] * SprayUnits::rhod_conv;
        // Get normalized vapor molar fraction
        amrex::Real X_vi = fdat.binary_csr_rowdotprod(
          pcspec, fdat.L_row, fdat.L_col, X_vn.data());
        amrex::Real X_norm = X_vi / sumXVap;
        rhoD_ri[pcspec] *= X_norm;
        rhoDtotal += rhoD_ri[pcspec];
        sumYSkin += Y_skin[pcspec];
      }

#ifndef USE_MANIFOLD_EOS
      amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_fuel_skin;
      Y_fuel_skin.fill(0.0);
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        Y_fuel_skin[fdat.dep_indx[spf]] = Y_skin[fdat.dep_indx[spf]] / sumYSkin;
      }
      eos.TY2Cp(T_skin, Y_fuel_skin.data(), cp_fuel);
#else
      eos.Y2GenericManifoldData(
        Y_skin.data(), &(fdat.eosparm->idx_Cp_fuel), &cp_fuel, 1);
      // TODO: account for variation between table temperature and skin
      // temperature
      //  cp_fuel = 1e7 + (cp_fuel - 1e7) * T_skin / T_skin_predicted;
#endif
      cp_fuel *= SprayUnits::eng_conv;

      amrex::Real Sc_skin = mu_skin / rhoDtotal;
      amrex::Real logB = std::log1p(B_M);
      amrex::Real invFM = B_M / (logB * std::pow(1. + B_M, 0.7));
      amrex::Real Sh_0 = 1. + powR * std::cbrt(1. + Reyn * Sc_skin);
      amrex::Real Sh_num = 2. + (Sh_0 - 2.) * invFM;
      amrex::Real mdotcoeff = M_PI * dia_part * Sh_num * logB;
      amrex::Real ratio = cp_fuel * Sh_num * rhoDtotal / lambda_skin;

      Nu_num = calcHeatCoeff(ratio, B_M, B_eps, C_eps, Nu_0);

      amrex::GpuArray<amrex::Real, NUM_SPECIES> L_n_pc{};
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        // Species index
        const int fdspec = fdat.dep_indx[spf];
        if (X_vn[spf] > 0.) {
          amrex::Real X_vi =
            fdat.binary_csr_rowdotprod(fdspec, fdat.L_row, fdat.L_col, X_vn);
          amrex::Real Xnorm = 1.;
          if (X_vi > X_vn[spf]) {
            // Multiple species deposit to fdspec, scale rhoD_{r,n}
            // accordingly
            Xnorm = X_vn[spf] / X_vi;
          }
          mn_dot[spf] =
            -amrex::max<amrex::Real>(Xnorm * rhoD_ri[fdspec] * mdotcoeff, 0.);
          L_n_pc[fdspec] = L_n[spf];
          m_dot += mn_dot[spf];
          sumL += mn_dot[spf] * L_n[spf];
        }
      }
    } else {
      // Still do convective heat transfer for saturated gas
      Nu_num = Nu_0;
    }

    conv_src = M_PI * lambda_skin * dia_part * delT * Nu_num;
    part_temp_src = (sumL + conv_src) * inv_pmass / cp_part;
  }
  rates.part_vel_src = part_vel_src;
  rates.part_temp_src = part_temp_src;
  rates.m_dot = m_dot;
  rates.drag_force = drag_force;
  rates.drag_force_p = drag_force_p;
  rates.conv_src = conv_src;
  rates.cp_part = cp_part;
  rates.delT = delT;
  rates.Reyn = Reyn;
}

#endif
