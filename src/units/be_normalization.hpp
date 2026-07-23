//========================================================================================
// AthenaPK - Bonnor-Ebert collapse normalization (flagship Phase 1, Option A).
//
// THE single derivation of the FHC code->cgs base scales from the physical IC
// (total mass [Msun], temperature [K], central-density enhancement f). It is shared by:
//   - the collapse_be problem generator (which sets the initial conditions), and
//   - PhysicalUnits::BuildPhysicalUnits (which the microphysics packages consume),
// so the code->cgs mapping used by the dynamics and by chemistry/radiation/non-ideal-MHD
// cannot drift apart. Historically the pgen used this exact derivation while PhysicalUnits
// fell back to rounded <units> defaults (rho 5.467e-19 vs 5.4668e-19, l 2.81e16 vs 2.8063e16),
// a ~0.13% mismatch in diffusivity_unit/opacity_unit. Routing both through this function
// drives that drift to zero.
//
// Dependency-free (only <cmath>) so it can be unit-tested host-only. The constants and the
// exact formula sequence match the historical collapse_be.cpp derivation bit-for-bit, so the
// pgen's initial conditions and sidecar are unchanged by adopting it.
//========================================================================================
#ifndef UNITS_BE_NORMALIZATION_HPP_
#define UNITS_BE_NORMALIZATION_HPP_

#include <cmath>

namespace PhysUnits {

namespace be_cgs {
constexpr double bemass_code = 197.561; // total mass of the critical BE sphere [code units]
constexpr double cs10 = 1.9e4;          // isothermal sound speed at 10 K [cm/s]
constexpr double msun = 1.9891e33;      // solar mass [g]
constexpr double G = 6.67259e-8;        // gravitational constant [cgs]
} // namespace be_cgs

//! Base cgs scales per code unit for the BE collapse IC. rho0*l0^3 == m0 and l0/v0 == t0 to
//! round-off (the normalization is self-consistent with four_pi_G = 1).
struct BENormalization {
  double rho0; // g/cm^3 per code density
  double v0;   // cm/s   per code velocity (== isothermal c_s)
  double l0;   // cm     per code length
  double t0;   // s      per code time
  double m0;   // g      per code mass
};

//! Derive the exact BE normalization from the physical IC. Same operation sequence as the
//! historical collapse_be.cpp block (mass = bemass*f in code units, c_s = 1 at temperature).
inline BENormalization DeriveBENormalization(double mass_msun, double temp_K, double f) {
  using namespace be_cgs;
  BENormalization b;
  b.m0 = mass_msun * msun / (bemass_code * f);
  b.v0 = cs10 * std::sqrt(temp_K / 10.0);
  b.rho0 = std::pow(b.v0, 6) / (b.m0 * b.m0) / (64.0 * M_PI * M_PI * M_PI * G * G * G);
  b.t0 = 1.0 / std::sqrt(4.0 * M_PI * G * b.rho0);
  b.l0 = b.v0 * b.t0;
  return b;
}

} // namespace PhysUnits

#endif // UNITS_BE_NORMALIZATION_HPP_
