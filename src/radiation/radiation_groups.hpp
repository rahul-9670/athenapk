//========================================================================================
// AthenaPK - M1 radiation: multigroup frequency structure (Phase 4 scaffold, increment 1).
//
// Foundation for multigroup RHD. This increment adds ONLY the frequency-group structure and
// the group-integrated Planck function (no transport changes yet). It is designed so that
// N_GROUP = 1 reduces EXACTLY to the current gray behavior (a single group spanning [0, inf)
// has Planck fraction 1, so B_group = a_R T^4), keeping production bit-identical until the
// multigroup path is wired into the moments update + opacity in later increments.
//
// Design (later increments): per-group moments (Er_g, Fr_g), per-group M1 closure (reuses
// radiation_closure.hpp unchanged -- the closure is frequency-independent), per-group Planck/
// Rosseland opacity means from a monochromatic dataset, and IMEX coupling summing the group
// emission/absorption. The M1 closure and signal speeds are per-group; only the source term
// couples groups (through the matter temperature).
//========================================================================================
#ifndef RADIATION_RADIATION_GROUPS_HPP_
#define RADIATION_RADIATION_GROUPS_HPP_

#include <cmath>

#include <Kokkos_Core.hpp>

#include <basic_types.hpp> // parthenon::Real

namespace Radiation {

using parthenon::Real;

// Max compile-time group count (POD-friendly; N_GROUP=1 default keeps the gray path).
constexpr int MAX_GROUP = 8;

//----------------------------------------------------------------------------------------
//! Cumulative Planck energy fraction F(x) = (15/pi^4) * Int_0^x u^3/(e^u - 1) du, x = h*nu/kT.
//! F(0)=0, F(inf)=1; the fraction in a group [x1,x2] is F(x2)-F(x1). The tail integral
//! Int_x^inf u^3/(e^u-1) du = sum_{n>=1} e^{-n x}( x^3/n + 3 x^2/n^2 + 6 x/n^3 + 6/n^4 )
//! (Clark 1965), and Int_0^inf = pi^4/15, so F(x) = 1 - (15/pi^4)*tail. The series converges
//! fast; at x=0 the tail sums to pi^4/15 (F=0), at large x it -> 0 (F=1).
KOKKOS_INLINE_FUNCTION Real PlanckCumFraction(const Real x) {
  if (x <= 0.0) return 0.0;
  constexpr Real k15_pi4 = 15.0 / (M_PI * M_PI * M_PI * M_PI);
  const Real x2 = x * x, x3 = x2 * x;
  Real tail = 0.0;
  for (int n = 1; n <= 30; ++n) {
    const Real nn = static_cast<Real>(n);
    tail += std::exp(-nn * x) * (x3 / nn + 3.0 * x2 / (nn * nn) + 6.0 * x / (nn * nn * nn) +
                                 6.0 / (nn * nn * nn * nn));
    if (nn * x > 60.0) break; // remaining terms negligible
  }
  Real f = 1.0 - k15_pi4 * tail;
  return (f < 0.0) ? 0.0 : ((f > 1.0) ? 1.0 : f);
}

//----------------------------------------------------------------------------------------
//! Frequency-group structure. POD captured by value into device kernels. Group boundaries are
//! stored as x = h*nu/(k*T_ref) reference values are NOT used -- instead the group boundaries
//! are physical frequencies [Hz]; the Planck fraction below uses the LOCAL matter T. Edges are
//! nu[0]=0 .. nu[n_group]=inf(large). N_GROUP=1 with nu=[0,inf) => gray.
struct RadGroups {
  int n_group = 1;
  // group frequency edges [Hz], size n_group+1; nu_edge[0]=0, nu_edge[n_group]=+inf (HUGE).
  Real nu_edge[MAX_GROUP + 1] = {0.0};
  Real h_over_k = 4.799243e-11; // h/k_B [s*K]  (Planck const / Boltzmann)

  //! Equilibrium Planck energy fraction in group g at matter temperature T [K].
  //! sum_g PlanckFraction(g,T) == 1 exactly (edges span [0,inf)).
  KOKKOS_INLINE_FUNCTION Real PlanckFraction(const int g, const Real T) const {
    if (n_group == 1) return 1.0; // gray: the single group holds all the energy
    const Real xlo = h_over_k * nu_edge[g] / T;
    const Real xhi = (g == n_group - 1) ? 1.0e30 : h_over_k * nu_edge[g + 1] / T;
    return PlanckCumFraction(xhi) - PlanckCumFraction(xlo);
  }
};

//----------------------------------------------------------------------------------------
//! Build a default (gray) or log-spaced multigroup structure. n_group=1 -> gray (bit-identical).
inline RadGroups BuildRadGroups(int n_group, Real nu_min_hz, Real nu_max_hz) {
  RadGroups g;
  if (n_group < 1) n_group = 1;
  if (n_group > MAX_GROUP) n_group = MAX_GROUP;
  g.n_group = n_group;
  g.nu_edge[0] = 0.0;
  if (n_group == 1) {
    g.nu_edge[1] = 1.0e30; // [0, inf) -> gray
    return g;
  }
  // interior edges log-spaced in [nu_min, nu_max]; first edge 0, last +inf so nothing is lost.
  const Real ll = std::log10(nu_min_hz), lh = std::log10(nu_max_hz);
  for (int i = 1; i < n_group; ++i)
    g.nu_edge[i] = std::pow(10.0, ll + (lh - ll) * (i - 1) / (n_group - 2 + 1e-30));
  g.nu_edge[n_group] = 1.0e30;
  return g;
}

} // namespace Radiation

#endif // RADIATION_RADIATION_GROUPS_HPP_
