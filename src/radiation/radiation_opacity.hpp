//========================================================================================
// AthenaPK - M1 moment radiation transport package.
// Gray opacity model + code-unit calibration for the FHC collapse problem.
//
// Two responsibilities:
//   (1) Convert the physical radiation constants / opacities into the run's SHARED
//       code units (the same rho0/v0/l0 normalization the two codes use), so that
//       arad and the absorption coefficient are dimensionally consistent with the
//       gas energy density (= pressure unit rho0*v0^2) and the M1 transport.
//   (2) Provide a device-callable gray opacity kappa(rho,T): either a constant
//       placeholder or a Bell & Lin (1994) low-temperature dust law, which is what
//       produces the radiative trapping that sets the first-hydrostatic-core entropy.
//========================================================================================
#ifndef RADIATION_RADIATION_OPACITY_HPP_
#define RADIATION_RADIATION_OPACITY_HPP_

#include <cmath>

#include <Kokkos_Core.hpp>

#include <basic_types.hpp> // parthenon::Real

#include "radiation_closure.hpp" // parthenon::Real alias, RadFuzz

namespace Radiation {

// CGS physical constants used only for the one-time code-unit calibration (host side).
namespace cgs {
constexpr Real a_rad = 7.5657e-15;     // radiation constant a_R [erg cm^-3 K^-4]
constexpr Real k_B = 1.380649e-16;     // Boltzmann [erg/K]
constexpr Real m_H = 1.6726219e-24;    // hydrogen mass [g]
constexpr Real c_light = 2.99792458e10; // speed of light [cm/s]
} // namespace cgs

enum class OpacityModel { constant, dust, belllin };

//----------------------------------------------------------------------------------------
//! Device-side opacity parameters (all in CODE units / conversion factors). Captured by
//! value into the matter-coupling kernel; no host pointers, safe on GPU.
struct OpacityParams {
  OpacityModel model = OpacityModel::constant;
  // constant model: code-unit opacities (per unit mass), coefficient = rho_code*kappa.
  Real kappa_a0 = 0.0;
  Real kappa_s0 = 0.0;
  // dust model (Bell & Lin 1994 ice-grain regime kappa = k0 * T_phys^2 [cm^2/g]):
  Real dust_k0_cgs = 2.0e-4; // [cm^2 g^-1 K^-2]
  Real dust_Tmax = 1.5e3;    // sublimation cap on T_phys [K] (kappa frozen above)
  Real kappa_s_dust = 0.0;   // optional gray scattering opacity (code units)
  // unit conversions:
  Real T_unit = 1.0;     // K per unit code temperature (T_phys = T_unit * T_code)
  Real rho_unit = 1.0;   // g/cm^3 per unit code density (rho_phys = rho_unit * rho_code)
  Real kappa_conv = 1.0; // kappa_code = kappa_phys[cm^2/g] * kappa_conv (= rho0*l0)
  // Planck/Rosseland split: Bell&Lin/dust laws are ROSSELAND fits; the emission (matter-
  // coupling) term should use the PLANCK mean, which for dust exceeds the Rosseland mean
  // (kappa_P/kappa_R ~ 2-5). Modeled as kappa_P = planck_ross_ratio * kappa_R. Default 1.0
  // => bit-identical to the pre-split single-opacity behavior. A full (log rho, log T) ->
  // (kappa_P, kappa_R) table (Semenov et al. 2003) is the follow-on increment.
  Real planck_ross_ratio = 1.0;
  // Regime-skip-robust Bell&Lin walk (fixes the low-rho/high-T kappa discontinuity). Default
  // false -> the plain walk (bit-identical to the shipped behavior AND to the fixed walk on the
  // collapse track; the two differ only in the off-track low-rho/high-T corner). Opt in via
  // <radiation> bell_lin_fix_regime_skip = true.
  bool bell_lin_fix_regime_skip = false;
};

//----------------------------------------------------------------------------------------
//! Bell & Lin (1994) Rosseland-mean gray opacity [cm^2/g] (rho,T in cgs). Eight piecewise
//! power-law regimes kappa = k0 * rho^a * T^b: ice grains, ice evaporation, metal/dust
//! grains, DUST SUBLIMATION (the opacity gap ~1500-2000 K where kappa plunges T^-24),
//! molecular, H- (negative hydrogen ion), bound-free+free-free (Kramers), and electron
//! scattering. This is the physical opacity through the SECOND-core regime (T>1500 K),
//! replacing the frozen-dust cap. Regime is chosen by walking up the density-dependent
//! transition temperatures (where adjacent regimes' opacities cross).
KOKKOS_INLINE_FUNCTION
Real BellLinKappa(const Real rho, const Real T) {
  constexpr Real k0[8] = {2.0e-4, 2.0e16, 0.1, 2.0e81, 1.0e-8, 1.0e-36, 1.5e20, 0.348};
  constexpr Real aa[8] = {0.0, 0.0, 0.0, 1.0, 2.0 / 3.0, 1.0 / 3.0, 1.0, 0.0};
  constexpr Real bb[8] = {2.0, -7.0, 0.5, -24.0, 3.0, 10.0, -2.5, 0.0};
  int i = 0;
  for (; i < 7; ++i) {
    // transition T where regime i meets i+1: k0[i] rho^aa[i] T^bb[i] = k0[i+1] rho^.. T^..
    const Real Tt = std::pow((k0[i] / k0[i + 1]) * std::pow(rho, aa[i] - aa[i + 1]),
                             1.0 / (bb[i + 1] - bb[i]));
    if (T < Tt) break;
  }
  return k0[i] * std::pow(rho, aa[i]) * std::pow(T, bb[i]);
}

//----------------------------------------------------------------------------------------
//! Regime-skip-robust Bell & Lin walk. The plain BellLinKappa walk assumes the adjacent-regime
//! transition temperatures are encountered in order; at low rho (<~1e-11) they go OUT OF ORDER
//! (the Kramers regime's window closes) and the walk jumps across a skipped regime, giving a
//! kappa DISCONTINUITY of up to ~4 decades at T~4-9 kK. This version skips any regime whose
//! forward window has closed [cross(j,j+1) <= cross(i,j)] and bridges regime i directly to the
//! next ACTIVE regime, so kappa is CONTINUOUS everywhere. VERIFIED: identical to BellLinKappa
//! on the collapse track (rho>=1e-10, rel diff 0) and gives the correct canonical values
//! (0.02 @ 10 K, 2.0 @ 100 K, 0.348 e-scatter); only the off-track corner changes.
KOKKOS_INLINE_FUNCTION
Real BellLinKappaFixed(const Real rho, const Real T) {
  constexpr Real k0[8] = {2.0e-4, 2.0e16, 0.1, 2.0e81, 1.0e-8, 1.0e-36, 1.5e20, 0.348};
  constexpr Real aa[8] = {0.0, 0.0, 0.0, 1.0, 2.0 / 3.0, 1.0 / 3.0, 1.0, 0.0};
  constexpr Real bb[8] = {2.0, -7.0, 0.5, -24.0, 3.0, 10.0, -2.5, 0.0};
  auto Tcross = [&](int p, int q) { // T where regime p == regime q at this rho
    return std::pow((k0[p] / k0[q]) * std::pow(rho, aa[p] - aa[q]), 1.0 / (bb[q] - bb[p]));
  };
  int i = 0;
  while (i < 7) {
    int j = i + 1;
    while (j < 7 && Tcross(j, j + 1) <= Tcross(i, j)) ++j; // skip closed-window regimes
    if (T < Tcross(i, j)) return k0[i] * std::pow(rho, aa[i]) * std::pow(T, bb[i]);
    i = j;
  }
  return k0[7] * std::pow(rho, aa[7]) * std::pow(T, bb[7]);
}

//----------------------------------------------------------------------------------------
//! ROSSELAND-mean absorption opacity in CODE units (per unit mass), kappa_R(rho,T). Used
//! for flux attenuation / the diffusion limit. Bell & Lin is itself a Rosseland fit, so
//! this is the base opacity; the coefficient in the flux term is rho_code * (this + scat).
KOKKOS_INLINE_FUNCTION
Real RosselandOpacity(const OpacityParams &op, const Real rho_code, const Real T_code) {
  if (op.model == OpacityModel::constant) return op.kappa_a0;
  const Real T_phys = op.T_unit * T_code;
  if (op.model == OpacityModel::belllin) {
    const Real rho_phys = op.rho_unit * rho_code;
    const Real kap = op.bell_lin_fix_regime_skip ? BellLinKappaFixed(rho_phys, T_phys)
                                                 : BellLinKappa(rho_phys, T_phys);
    return kap * op.kappa_conv;
  }
  // dust: Bell & Lin low-T ice-grain law only, frozen above sublimation.
  const Real Td = std::min(T_phys, op.dust_Tmax);
  const Real kappa_phys = op.dust_k0_cgs * Td * Td; // [cm^2/g]
  return kappa_phys * op.kappa_conv;
}

//----------------------------------------------------------------------------------------
//! PLANCK-mean absorption opacity in CODE units (per unit mass), kappa_P(rho,T). Used for
//! the emission/absorption matter coupling S_E = chat*rho*kappa_P*(E-B). Modeled as
//! kappa_P = planck_ross_ratio * kappa_R (ratio=1 -> identical to kappa_R). See
//! OpacityParams::planck_ross_ratio.
KOKKOS_INLINE_FUNCTION
Real PlanckOpacity(const OpacityParams &op, const Real rho_code, const Real T_code) {
  return op.planck_ross_ratio * RosselandOpacity(op, rho_code, T_code);
}

//----------------------------------------------------------------------------------------
//! Gray SCATTERING opacity in CODE units (per unit mass). Constant for both models.
KOKKOS_INLINE_FUNCTION
Real ScatteringOpacity(const OpacityParams &op, const Real rho_code, const Real T_code) {
  return (op.model == OpacityModel::constant) ? op.kappa_s0 : op.kappa_s_dust;
}

} // namespace Radiation

#endif // RADIATION_RADIATION_OPACITY_HPP_
