//========================================================================================
// AthenaPK - M1 moment radiation transport package.
// M1 closure relations ported faithfully from Artemis (LANL, BSD)
// src/radiation/moments/moments.hpp. Licensed BSD 3-Clause.
//
// Pure device-callable math (no AthenaPK/Artemis data structures): the Eddington
// factor/tensor, M1 signal speeds, and flux normalization that the transport
// (increment 2b) and matter coupling (increment 3) build on.
//========================================================================================
#ifndef RADIATION_RADIATION_CLOSURE_HPP_
#define RADIATION_RADIATION_CLOSURE_HPP_

#include <array>
#include <cmath>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include <basic_types.hpp> // parthenon::Real

namespace Radiation {

using parthenon::Real;

enum class Closure { P1, M1 };

// Square + a tiny floor to guard 0/0 in flux-direction normalization (Artemis "Fuzz").
KOKKOS_INLINE_FUNCTION Real RadSqr(const Real x) { return x * x; }
KOKKOS_FORCEINLINE_FUNCTION Real RadFuzz() { return 1.0e-100; }

//----------------------------------------------------------------------------------------
//! 3x the Eddington factor chi(f) for the given closure (f = reduced flux |F|/(c E)).
template <Closure CTYP>
KOKKOS_INLINE_FUNCTION Real ThriceEddingtonFactor(const Real f) {
  if constexpr (CTYP == Closure::P1) {
    return 1.0;
  } else { // M1
    const Real f2 = f * f;
    return 3.0 * (3.0 + 4.0 * f2) / (5.0 + 2.0 * std::sqrt(4.0 - 3.0 * f2));
  }
}

//----------------------------------------------------------------------------------------
//! Symmetric Eddington tensor P/E given the reduced-flux vector fred[3].
//! Returns {xx, yy, zz, yz, xz, xy}.
template <Closure CTYP>
KOKKOS_INLINE_FUNCTION std::array<Real, 6>
EddingtonTensor(const std::array<Real, 3> &fred) {
  if constexpr (CTYP == Closure::P1) {
    return {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0};
  } else { // M1
    Real fmag = std::sqrt(RadSqr(fred[0]) + RadSqr(fred[1]) + RadSqr(fred[2]));
    const std::array<Real, 3> n{fred[0] / (fmag + RadFuzz()),
                                fred[1] / (fmag + RadFuzz()),
                                fred[2] / (fmag + RadFuzz())};
    fmag = std::min(1.0, fmag);
    const Real chi = ThriceEddingtonFactor<CTYP>(fmag);
    const Real ca = (3.0 - chi) / 6.0;
    const Real cb = 0.5 * (chi - 1.0);
    return {ca + cb * n[0] * n[0], ca + cb * n[1] * n[1], ca + cb * n[2] * n[2],
            cb * n[1] * n[2],      cb * n[0] * n[2],      cb * n[0] * n[1]};
  }
}

//----------------------------------------------------------------------------------------
//! M1 signal speeds (min,max) along a face, given mu = n.nhat and f = |reduced flux|.
//! Returned in units of chat (multiply by chat downstream). Kokkos::pair is device-safe.
template <Closure CTYP>
KOKKOS_INLINE_FUNCTION Kokkos::pair<Real, Real> WaveSpeed(const Real mu, const Real f) {
  if constexpr (CTYP == Closure::P1) {
    const Real val = std::sqrt(1.0 / 3.0);
    return {-val, val};
  } else { // M1
    const Real f2 = f * f;
    const Real det = 4.0 - 3.0 * f2;
    const Real sdet = std::sqrt(det);
    // The radical is analytically >= 0 for f in [0,1] but can dip slightly negative
    // by roundoff near f=1 (beam limit); clamp to avoid NaN signal speeds.
    const Real arg = 2.0 / 3.0 * (det - sdet) + 2.0 * mu * mu * (2.0 - f2 - sdet);
    const Real fac = std::sqrt(std::max(0.0, arg));
    const Real norm = 1.0 / (sdet + RadFuzz());
    return {norm * (mu * f - fac), norm * (mu * f + fac)};
  }
}

//----------------------------------------------------------------------------------------
//! Clamp the reduced-flux vector to the causal limit |f| <= 1, preserving direction.
KOKKOS_INLINE_FUNCTION
std::array<Real, 3> NormalizeFlux(const Real fx1, const Real fx2, const Real fx3) {
  Real f = std::sqrt(RadSqr(fx1) + RadSqr(fx2) + RadSqr(fx3));
  const Real nx1 = fx1 / (f + RadFuzz());
  const Real nx2 = fx2 / (f + RadFuzz());
  const Real nx3 = fx3 / (f + RadFuzz());
  f = std::min(1.0, f);
  return {nx1 * f, nx2 * f, nx3 * f};
}

//----------------------------------------------------------------------------------------
//! Fleck factor dB/dE = 4 a T^3 / cv (used by the implicit matter coupling, incr. 3).
KOKKOS_INLINE_FUNCTION
Real FleckFactor(const Real ar, const Real T, const Real cv) {
  return 4.0 * ar * T * T * T / cv;
}

} // namespace Radiation

#endif // RADIATION_RADIATION_CLOSURE_HPP_
