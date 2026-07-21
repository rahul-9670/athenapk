//========================================================================================
// AthenaPK - WS-2 gas thermochemistry (heating/cooling) for the M1-RT collapse envelope.
//
// Provides Gamma - Lambda [erg cm^-3 s^-1] as a device/host function of (n_H, T, abundances,
// T_dust), following the network-struct pattern of network_gow17_reduced.hpp. Rates are the
// standard molecular-cloud set (Gong, Ostriker & Wolfire 2017 and refs therein):
//   Heating:  CR ionization, H2 grain-formation, photoelectric (A_V-shielded).
//   Cooling:  [C II] 158um, [O I] 63um, CO J=1-0 (simplified), gas-grain (dust) exchange.
// The gas-grain term Lambda_gd ~ n_H^2 (T - T_dust) DOMINATES for n >~ 10^4.5 and pins the
// gas temperature to T_dust in the core; C+/O/CO line cooling set the ~10-15 K envelope.
//
// UNITS: this file works in cgs (n_H in cm^-3, T in K). The caller converts rho_code->n_H
// (network nH()) and e_code<->cgs. Design constraint: WS-2 is only active when RT is on
// (RT owns e_th; the barotropic overwrite would erase any Gamma-Lambda work).
//========================================================================================
#ifndef CHEMISTRY_THERMO_HPP_
#define CHEMISTRY_THERMO_HPP_

#ifdef KOKKOS_INLINE_FUNCTION
#define CHEMT_FN KOKKOS_INLINE_FUNCTION
#else
#define CHEMT_FN inline
#endif

#include <cmath>

namespace Chemistry {

namespace thermo_cgs {
constexpr double eV = 1.602176634e-12; // erg
constexpr double kB = 1.380649e-16;    // erg/K
} // namespace thermo_cgs

//----------------------------------------------------------------------------------------
//! Two-level-atom volumetric cooling per unit collider density, times n(species). Returns
//! Lambda [erg cm^-3 s^-1]. n_coll is the collider density (~ n_H for neutral gas).
//!   Lambda = n(sp) * dE * A_ul * (gu/gl)e^{-T0/T} / (1 + (gu/gl)e^{-T0/T} + n_cr/n_coll)
//! with n_cr = A_ul/gamma_ul. Captures Boltzmann suppression (cold gas) and the LTE
//! saturation above the critical density.
CHEMT_FN double TwoLevelCool(double n_sp, double n_coll, double T, double T0, double A_ul,
                             double dE, double gu_gl, double gamma_ul) {
  if (n_sp <= 0.0 || T <= 0.0) return 0.0;
  const double exc = gu_gl * std::exp(-T0 / T);
  const double n_cr = A_ul / gamma_ul;
  return n_sp * dE * A_ul * exc / (1.0 + exc + n_cr / n_coll);
}

//----------------------------------------------------------------------------------------
//! Thermochemistry parameters (cgs). Defaults are the FHC molecular-cloud set.
struct ThermoParams {
  // --- Heating ---
  double q_cr_erg = 20.0 * thermo_cgs::eV; // heat deposited per CR ionization [erg]
  double zeta = 1.0e-17;                   // primary CR ionization rate [s^-1] (= network)
  double h2_heat_erg = 0.2 * thermo_cgs::eV; // grain-borne H2 formation heat [erg]
  double kgr = 3.0e-17;                      // H2 grain formation [cm^3 s^-1] (= network)
  double pe_rate0 = 1.3e-24;   // photoelectric heating coeff [erg s^-1] (Bakes-Tielens)
  double pe_eps = 0.05;        // PE efficiency (envelope value; ~const for FHC)
  double G0 = 1.7;             // FUV field in Habing units (Draine ~1.7 G_Habing)
  // A_V(n_H) shielding proxy for the standalone curve: A_V = av0*(n_H/av_n0)^av_exp.
  // G_eff = G0 exp(-2.5 A_V) -> PE heating shuts off in the shielded (dense) gas.
  double av0 = 0.5, av_n0 = 1.0e3, av_exp = 2.0 / 3.0;
  // --- Cooling: [C II] 158um ---
  double cII_T0 = 91.2, cII_A = 2.4e-6, cII_dE = 1.259e-14, cII_gugl = 2.0, cII_gamma = 8.0e-10;
  // --- Cooling: [O I] 63um (fixed abundance) ---
  double x_O = 3.2e-4;
  double oI_T0 = 227.7, oI_A = 8.9e-5, oI_dE = 3.144e-14, oI_gugl = 0.6, oI_gamma = 5.0e-10;
  // --- Cooling: CO J=1-0 (simplified 2-level; higher-J ladder omitted) ---
  double co_T0 = 5.53, co_A = 7.2e-8, co_dE = 7.63e-16, co_gugl = 3.0, co_gamma = 3.3e-11;
  // CO line-trapping (escape probability): CO becomes optically THICK in dense gas, which
  // strongly suppresses its cooling (the optically-thin 2-level above would over-cool and
  // drive T below T_dust at n~10^4). tau_CO = co_tau_coeff * x_CO * n_H / sqrt(T/10);
  // beta = (1-e^{-tau})/tau. The coefficient encodes the (0-D unmodeled) CO column and is
  // calibrated so the FHC envelope sits at ~10-15 K at n=10^4 (Goldsmith 2001-style).
  double co_tau_coeff = 5.5;
  // --- Gas-grain (dust) energy exchange ---
  double alpha_gd = 3.2e-34; // [erg cm^3 K^-3/2]  Lambda_gd = alpha_gd n_H^2 sqrt(T)(T-Td)
};

//----------------------------------------------------------------------------------------
//! Net volumetric heating rate Gamma - Lambda [erg cm^-3 s^-1] at (n_H, T, abundances,
//! T_dust). Abundances x_* = n_*/n_H. x_H = 1 - 2 x_H2 - x_Hp (neutral H nuclei fraction).
CHEMT_FN double NetHeatCool(const ThermoParams &p, double n_H, double T, double x_H,
                            double x_Cp, double x_CO, double T_dust) {
  // --- Heating ---
  const double G_cr = p.q_cr_erg * p.zeta * n_H;
  const double G_h2 = p.h2_heat_erg * p.kgr * n_H * n_H * (x_H > 0.0 ? x_H : 0.0);
  const double A_V = p.av0 * std::pow(n_H / p.av_n0, p.av_exp);
  const double G_eff = p.G0 * std::exp(-2.5 * A_V);
  const double G_pe = p.pe_rate0 * n_H * p.pe_eps * G_eff;
  const double Gamma = G_cr + G_h2 + G_pe;

  // --- Line cooling (2-level atoms/molecule) ---
  const double L_cII =
      TwoLevelCool(x_Cp * n_H, n_H, T, p.cII_T0, p.cII_A, p.cII_dE, p.cII_gugl, p.cII_gamma);
  const double L_oI =
      TwoLevelCool(p.x_O * n_H, n_H, T, p.oI_T0, p.oI_A, p.oI_dE, p.oI_gugl, p.oI_gamma);
  double L_CO =
      TwoLevelCool(x_CO * n_H, n_H, T, p.co_T0, p.co_A, p.co_dE, p.co_gugl, p.co_gamma);
  // CO escape probability (line trapping) -> optically-thick suppression at high column.
  const double tau_CO = p.co_tau_coeff * x_CO * n_H / std::sqrt(T / 10.0);
  const double beta_CO = (tau_CO > 1.0e-4) ? (1.0 - std::exp(-tau_CO)) / tau_CO : 1.0;
  L_CO *= beta_CO;
  // --- Gas-grain (drives T -> T_dust at high n) ---
  const double L_gd = p.alpha_gd * n_H * n_H * std::sqrt(T) * (T - T_dust);
  const double Lambda = L_cII + L_oI + L_CO + L_gd;

  return Gamma - Lambda;
}

//----------------------------------------------------------------------------------------
//! Advance the gas internal-energy density e_code [code units] over dt_code under the
//! gas thermochemistry Gamma-Lambda, at fixed rho and abundances. SEMI-IMPLICIT sub-cycler
//! (numerical-Jacobian backward Euler): e_new = e + dt*R/(1 - dt*dR/de) with dR/de < 0 (the
//! cooling self-stabilizes), so it is unconditionally stable and relaxes toward T_eq without
//! overshoot. Sub-step is accuracy-limited (|de| <= cfl*e) and floored at dt/nsub_max so the
//! full dt is always covered in <= nsub_max steps. Conversions: T_K = gm1*e/rho * T_unit;
//! rate_to_code converts Gamma-Lambda [erg cm^-3 s^-1] to code energy-density rate
//! (= t_unit / (rho_unit*v_unit^2)). Returns new e_code; *ntrunc += 1 if nsub_max exhausted.
CHEMT_FN double AdvanceThermoEnergy(const ThermoParams &p, double rho_code, double e_code,
                                    double n_H, double x_H, double x_Cp, double x_CO,
                                    double T_dust, double gm1, double T_unit,
                                    double rate_to_code, double dt_code, int nsub_max,
                                    double cfl, double e_floor, int *ntrunc) {
  double e = (e_code > e_floor) ? e_code : e_floor;
  const double dt_floor = dt_code / static_cast<double>(nsub_max);
  double t = 0.0;
  int nsub = 0;
  while (t < dt_code && nsub < nsub_max) {
    double T_K = gm1 * e / rho_code * T_unit;
    if (T_K < 1.0e-4) T_K = 1.0e-4;
    const double R = NetHeatCool(p, n_H, T_K, x_H, x_Cp, x_CO, T_dust) * rate_to_code;
    // numerical Jacobian dR/de (one extra eval); dR/de < 0 near/away from equilibrium.
    const double de = 1.0e-4 * e + 1.0e-30;
    double T2 = gm1 * (e + de) / rho_code * T_unit;
    if (T2 < 1.0e-4) T2 = 1.0e-4;
    const double R2 = NetHeatCool(p, n_H, T2, x_H, x_Cp, x_CO, T_dust) * rate_to_code;
    double dRde = (R2 - R) / de;
    if (dRde > 0.0) dRde = 0.0; // guard: keep the implicit step contractive
    // accuracy-limited sub-step (|Delta e| <= cfl*e); near equilibrium R->0 => take the
    // full remaining dt. Floored so stiffness cannot stall progress.
    double dt_sub = dt_code - t;
    const double absR = (R < 0.0 ? -R : R);
    if (absR > 0.0) {
      const double dt_acc = cfl * e / absR;
      if (dt_acc < dt_sub) dt_sub = dt_acc;
    }
    if (dt_sub < dt_floor) dt_sub = dt_floor;
    if (t + dt_sub > dt_code) dt_sub = dt_code - t;
    e += dt_sub * R / (1.0 - dt_sub * dRde);
    if (e < e_floor) e = e_floor;
    t += dt_sub;
    ++nsub;
  }
  if (t < dt_code && ntrunc != nullptr) *ntrunc += 1;
  return e;
}

//----------------------------------------------------------------------------------------
//! Equilibrium temperature T_eq(n_H): bisection root of NetHeatCool = 0 in [Tlo, Thi] K.
//! Net is monotonically DECREASING in T over the relevant range (cooling rises with T),
//! so a sign change brackets a unique root. Returns Tlo/Thi if no bracket (pinned).
CHEMT_FN double SolveTeq(const ThermoParams &p, double n_H, double x_H, double x_Cp,
                         double x_CO, double T_dust, double Tlo = 3.0, double Thi = 1.0e4) {
  double flo = NetHeatCool(p, n_H, Tlo, x_H, x_Cp, x_CO, T_dust);
  double fhi = NetHeatCool(p, n_H, Thi, x_H, x_Cp, x_CO, T_dust);
  if (flo <= 0.0) return Tlo; // cooling wins even at Tlo -> pinned near/below T_dust
  if (fhi >= 0.0) return Thi; // heating wins even at Thi (should not happen for FHC)
  for (int it = 0; it < 100; ++it) {
    const double Tm = 0.5 * (Tlo + Thi);
    const double fm = NetHeatCool(p, n_H, Tm, x_H, x_Cp, x_CO, T_dust);
    if (fm > 0.0) {
      Tlo = Tm;
      flo = fm;
    } else {
      Thi = Tm;
      fhi = fm;
    }
    if (Thi - Tlo < 1.0e-4 * Tm) break;
  }
  return 0.5 * (Tlo + Thi);
}

} // namespace Chemistry

#endif // CHEMISTRY_THERMO_HPP_
