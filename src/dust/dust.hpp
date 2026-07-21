//========================================================================================
// AthenaPK - WS-4 single-moment dust evolution (device/host).
//
// Evolves two passive scalars: the dust-to-gas mass ratio f_dg and the mass-weighted
// characteristic grain radius a_c [cm]. Grains with a_c < um have Stokes number << 1 at
// core densities => they co-move with the gas (no drift); what matters is (i) GROWTH
// (raises a_c -> lowers total cross-section -> lowers kappa, raises the non-ideal etas) and
// (ii) SUBLIMATION above T_subl. This is the monodisperse sweep-up model (NOT the Artemis
// two-fluid dust + coagulation, which is the drift-limited fallback).
//
//   da_c/dt = (f_dg * rho_gas / rho_grain) * v_rel / 4          (monodisperse sweep-up)
//   v_rel   = sqrt(v_Brownian^2 + v_turb^2)
//     v_Brownian = sqrt(8 k_B T / (pi m_red)),  m_red = m_grain/2, m_grain=(4/3)pi a^3 rho_gr
//     v_turb     = sqrt(alpha_turb) * c_s * St,  St = (rho_grain a_c)/(rho_gas v_th) * Omega_ref
//   f_dg -> 0 smoothly above T_subl (sublimation), restored below (no re-formation lag, v1).
//
// UNITS: this header works in cgs (rho_gas g/cm^3, T K, a_c cm, c_s cm/s). The caller
// converts code<->cgs. All growth is capped at <=10%/step (the a_c->ionization->dt coupling
// risk) inside integrate_cell.
//========================================================================================
#ifndef DUST_DUST_HPP_
#define DUST_DUST_HPP_

#ifdef KOKKOS_INLINE_FUNCTION
#define DUST_FN KOKKOS_INLINE_FUNCTION
#else
#define DUST_FN inline
#endif

#include <cmath>

namespace Dust {

namespace dcgs {
constexpr double kB = 1.380649e-16;    // erg/K
constexpr double pi = 3.14159265358979323846;
} // namespace dcgs

//----------------------------------------------------------------------------------------
//! Single-moment dust parameters (cgs unless noted). Defaults are FHC-representative.
struct DustModel {
  double rho_grain = 3.0;    // internal grain density [g/cm^3] (silicate ~3)
  double a_ref = 1.0e-5;     // reference grain radius [cm] (0.1 um, MRN a_max-ish)
  double f_dg_ref = 0.01;    // canonical dust-to-gas ratio (ISM)
  double alpha_turb = 1.0e-2; // turbulent alpha for v_turb
  double omega_ref = 1.0e-13; // reference frequency [s^-1] setting St (~ 1/t_ff at n~1e6)
  double mu_gas = 2.33;      // mean molecular weight (for v_th, c_s)
  double m_H = 1.6726219e-24;
  double T_subl = 1.5e3;     // sublimation temperature [K]
  double subl_width = 100.0; // K, smoothness of the sublimation switch
  double a_floor = 1.0e-7;   // minimum grain radius [cm] (1 nm)
  double growth_cap = 0.10;  // max fractional a_c change per sub-step (dt-thrash guard)
  int    nsub_max = 200;
  // code-unit conversions (for the in-code path; the standalone test works in cgs directly)
  double rho_unit = 5.467e-19;
  double t_unit = 1.476998822551e12;

  //! Relative grain-grain velocity [cm/s] at (rho_gas, T, a_c).
  DUST_FN double v_rel(double rho_gas, double T, double a_c) const {
    const double cs = std::sqrt(dcgs::kB * T / (mu_gas * m_H)); // isothermal sound speed
    const double v_th = std::sqrt(8.0 / dcgs::pi) * cs;         // mean thermal speed
    const double m_grain = (4.0 / 3.0) * dcgs::pi * a_c * a_c * a_c * rho_grain;
    const double m_red = 0.5 * m_grain;
    const double v_B = std::sqrt(8.0 * dcgs::kB * T / (dcgs::pi * m_red));
    // Epstein Stokes number: St = t_stop * Omega_ref, t_stop = rho_grain a_c/(rho_gas v_th).
    const double St = (rho_grain * a_c / (rho_gas * v_th)) * omega_ref;
    const double v_t = std::sqrt(alpha_turb) * cs * St;
    return std::sqrt(v_B * v_B + v_t * v_t);
  }

  //! da_c/dt [cm/s] (monodisperse sweep-up) at fixed gas state and f_dg.
  DUST_FN double dadt(double rho_gas, double T, double a_c, double f_dg) const {
    return (f_dg * rho_gas / rho_grain) * v_rel(rho_gas, T, a_c) / 4.0;
  }

  //! Sublimation factor in [0,1]: 1 below T_subl, ->0 above (smooth tanh switch).
  DUST_FN double subl_factor(double T) const {
    return 0.5 * (1.0 - std::tanh((T - T_subl) / subl_width));
  }

  //! Advance (a_c, f_dg) over dt_s [s] at fixed gas state. a_c grows by monodisperse sweep-up
  //! (RK2 midpoint; sub-step capped at growth_cap/step to bound the a_c->ionization->dt
  //! coupling). f_dg is slaved to f_dg_ref*subl_factor(T) (sublimation/restore, no re-formation
  //! lag = v1). Writes a_c, f_dg in place; returns nsub used.
  DUST_FN int integrate_cell(double &a_c, double &f_dg, double rho_gas, double T,
                             double dt_s) const {
    // f_dg sublimation/restore (algebraic, no lag): 0.01 below T_subl, ->0 above.
    f_dg = f_dg_ref * subl_factor(T);
    if (a_c < a_floor) a_c = a_floor; // guard v_Brownian ~ 1/sqrt(a^3) against a_c -> 0
    double t = 0.0;
    int nsub = 0;
    const double dt_floor = dt_s / static_cast<double>(nsub_max);
    while (t < dt_s && nsub < nsub_max) {
      const double rate = dadt(rho_gas, T, a_c, f_dg);
      double dt_sub = dt_s - t;
      if (rate > 0.0) {
        const double dt_cap = growth_cap * a_c / rate; // |da|/a <= growth_cap
        if (dt_cap < dt_sub) dt_sub = dt_cap;
      }
      if (dt_sub < dt_floor) dt_sub = dt_floor;
      if (t + dt_sub > dt_s) dt_sub = dt_s - t;
      // RK2 midpoint (2nd order -> tracks exponential/power-law growth accurately).
      const double a_mid = a_c + 0.5 * dt_sub * rate;
      const double rate_mid = dadt(rho_gas, T, a_mid, f_dg);
      a_c += dt_sub * rate_mid;
      if (a_c < a_floor) a_c = a_floor;
      t += dt_sub;
      ++nsub;
    }
    return nsub;
  }
};

//----------------------------------------------------------------------------------------
//! Dust consumer multiplier for opacity / grain cross-section (WS-4 inc2). In the geometric
//! limit kappa ∝ total grain area ∝ (dust mass) / a  => factor = (f_dg/f_ref)*(a_ref/a_c).
//! At the reference state (f_dg=f_ref, a_c=a_ref) this is EXACTLY 1.0 => wiring it into
//! AbsorptionOpacity/RosselandOpacity and the ionization grain bins is bit-identical to the
//! static-grain code when grains are frozen at reference. Guarded against a_c -> 0.
DUST_FN double DustFactor(double f_dg, double a_c, double f_ref, double a_ref) {
  const double a = (a_c > 1.0e-30) ? a_c : 1.0e-30;
  return (f_dg / f_ref) * (a_ref / a);
}

} // namespace Dust

#endif // DUST_DUST_HPP_
