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
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

#include <basic_types.hpp> // parthenon::Real

namespace Radiation {

using parthenon::Real;

// Max compile-time group count (POD-friendly; N_GROUP=1 default keeps the gray path).
constexpr int MAX_GROUP = 8;

//----------------------------------------------------------------------------------------
//! Band-limited moment of the Planck weight:  Int_a^b x^s / (e^x - 1) dx  (s = 3 gives the
//! Planck energy weight; s = 3+beta weights a kappa_nu ~ nu^beta opacity for the Planck mean).
//! Composite Simpson in x; the integrand ~ x^{s-1} as x->0 (finite for s>=1) and decays like
//! x^s e^-x for large x. Device-callable.
//! IMPORTANT: the moments are used ONLY inside same-band ratios (num/den with the same lower
//! edge `a`), so a common factor e^{-a} is deliberately FACTORED OUT: the integrand carries
//! e^{a-x} (<=1 on [a,b]) instead of e^{-x}. This makes the ratio robust even in the deep Wien
//! tail (a >> 1) where e^{-x} underflows to 0 and would otherwise give 0/0 -> a spurious 0
//! opacity for the hardest group. The returned value equals e^{a} * (true moment); every use
//! divides two of these with the same a, so the e^{a} cancels exactly.
KOKKOS_INLINE_FUNCTION Real PlanckMoment(const Real s, Real a, Real b, const int npts) {
  if (a < 1.0e-8) a = 1.0e-8; // avoid the removable 1/x singularity of the weight at x=0
  if (b <= a) return 0.0;
  const int n = (npts % 2 == 0) ? npts : npts + 1; // Simpson needs an even number of panels
  const Real h = (b - a) / n;
  Real sum = 0.0;
  for (int k = 0; k <= n; ++k) {
    const Real x = a + k * h;
    const Real w = (k == 0 || k == n) ? 1.0 : ((k % 2) ? 4.0 : 2.0);
    // x^s/(e^x-1), with e^{-a} factored out: x^s e^{a-x}/(1-e^{-x}).
    sum += w * std::pow(x, s) * std::exp(a - x) / (1.0 - std::exp(-x));
  }
  return sum * h / 3.0;
}

//! Band-limited moment of the Rosseland weight:  Int_a^b x^s e^x/(e^x-1)^2 dx  (s = 4 gives the
//! Rosseland weight dB/dT; s = 4-beta forms the harmonic (Rosseland) mean of kappa_nu ~ nu^beta).
//! Same e^{-a}-factored-out ratio convention as PlanckMoment (device-callable).
KOKKOS_INLINE_FUNCTION Real RossMoment(const Real s, Real a, Real b, const int npts) {
  if (a < 1.0e-8) a = 1.0e-8;
  if (b <= a) return 0.0;
  const int n = (npts % 2 == 0) ? npts : npts + 1;
  const Real h = (b - a) / n;
  Real sum = 0.0;
  for (int k = 0; k <= n; ++k) {
    const Real x = a + k * h;
    const Real w = (k == 0 || k == n) ? 1.0 : ((k % 2) ? 4.0 : 2.0);
    const Real d = 1.0 - std::exp(-x);
    // x^s e^x/(e^x-1)^2, with e^{-a} factored out: x^s e^{a-x}/(1-e^{-x})^2.
    sum += w * std::pow(x, s) * std::exp(a - x) / (d * d);
  }
  return sum * h / 3.0;
}

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
  // Monochromatic dust opacity spectral index: kappa_nu(rho,T) = kappa_gray(rho,T) * phi(nu),
  // phi(nu) = (nu/nu_ref)^beta, normalized so the FULL-SPECTRUM Planck (resp. Rosseland) mean
  // of kappa_nu equals kappa_gray. Then the per-group Planck/Rosseland means are the band-limited
  // averages below. beta=0 => phi=1 => every group mean == kappa_gray (gray-per-group; the
  // equivalence gate holds exactly). Filled by FillGroupOpacity (host, at package init).
  Real beta = 0.0;
  Real cP_full = 1.0; // full-band Planck normalization  I_{3+beta}(0,inf)/I_3(0,inf)
  Real cR_full = 1.0; // full-band Rosseland normalization J_4(0,inf)/J_{4-beta}(0,inf)
  int qpts = 64;      // Simpson panels per band moment
  Real nquad_hi = 60.0; // integrate at most x_lo+60 in each band (weight ~ e^-x beyond => lost part negligible)

  //! Equilibrium Planck energy fraction in group g at matter temperature T [K].
  //! sum_g PlanckFraction(g,T) == 1 exactly (edges span [0,inf)).
  KOKKOS_INLINE_FUNCTION Real PlanckFraction(const int g, const Real T) const {
    if (n_group == 1) return 1.0; // gray: the single group holds all the energy
    const Real xlo = h_over_k * nu_edge[g] / T;
    const Real xhi = (g == n_group - 1) ? 1.0e30 : h_over_k * nu_edge[g + 1] / T;
    return PlanckCumFraction(xhi) - PlanckCumFraction(xlo);
  }

  //! group g's band edges in x = h*nu/(k*T) at matter temperature Tk [K]. The upper edge is
  //! capped at x_lo + nquad_hi so the fixed-panel Simpson keeps its resolution where the weight
  //! lives (the integrand decays like e^-x, so a very wide band [x_lo, x_hi>>x_lo] would waste
  //! all panels on the negligible tail and misestimate the ratio). The lost tail is ~e^-nquad_hi.
  KOKKOS_INLINE_FUNCTION void BandX(const int g, const Real Tk, Real &a, Real &b) const {
    a = h_over_k * nu_edge[g] / Tk;
    const Real bphys = (g == n_group - 1) ? (a + nquad_hi) : h_over_k * nu_edge[g + 1] / Tk;
    b = (bphys < a + nquad_hi) ? bphys : (a + nquad_hi);
  }

  //! PLANCK-mean opacity multiplier for group g at matter temperature Tk [K]:
  //!   kappa_P,g / kappa_gray = [ Int_g x^{3+beta}/(e^x-1) dx / Int_g x^3/(e^x-1) dx ] / cP_full.
  //! Sum-weighted by the group Planck fractions this reproduces the gray Planck emission exactly.
  KOKKOS_INLINE_FUNCTION Real PlanckBandMult(const int g, const Real Tk) const {
    if (n_group == 1 || beta == 0.0) return 1.0;
    Real a, b;
    BandX(g, Tk, a, b);
    const Real num = PlanckMoment(3.0 + beta, a, b, qpts);
    const Real den = PlanckMoment(3.0, a, b, qpts);
    return (num / (den + 1.0e-300)) / cP_full;
  }

  //! ROSSELAND-mean opacity multiplier for group g at matter temperature Tk [K]:
  //!   kappa_R,g / kappa_gray = [ Int_g x^4 w dx / Int_g x^{4-beta} w dx ] / cR_full,
  //! w = e^x/(e^x-1)^2 (the harmonic/Rosseland average of kappa_nu ~ nu^beta over the band).
  KOKKOS_INLINE_FUNCTION Real RossBandMult(const int g, const Real Tk) const {
    if (n_group == 1 || beta == 0.0) return 1.0;
    Real a, b;
    BandX(g, Tk, a, b);
    const Real num = RossMoment(4.0, a, b, qpts);
    const Real den = RossMoment(4.0 - beta, a, b, qpts);
    return (num / (den + 1.0e-300)) / cR_full;
  }
};

//----------------------------------------------------------------------------------------
//! Configure the monochromatic-dust opacity model kappa_nu ~ nu^beta on the group structure:
//! store beta and the full-spectrum normalizations so the band means reduce to kappa_gray when
//! summed over all frequencies (=> beta=0 is exactly gray-per-group). Host-side, once at init.
inline void FillGroupOpacity(RadGroups &g, const Real beta) {
  g.beta = beta;
  if (beta == 0.0) {
    g.cP_full = 1.0;
    g.cR_full = 1.0;
    return;
  }
  // full-band [~0, X] moments (X large enough that the e^-x tail is negligible).
  const Real X = 80.0;
  g.cP_full = PlanckMoment(3.0 + beta, 1.0e-8, X, 4000) /
              PlanckMoment(3.0, 1.0e-8, X, 4000);
  g.cR_full = RossMoment(4.0, 1.0e-8, X, 4000) /
              RossMoment(4.0 - beta, 1.0e-8, X, 4000);
}

//----------------------------------------------------------------------------------------
//! HDF5/pack field names for group g's four M1 moments (Er, Fr1, Fr2, Fr3). Group 0 keeps
//! the ORIGINAL gray names ("rad.Er", ...) so a gray (n_group=1) run's output/restart is
//! byte-for-byte unchanged; groups g>0 get a "_gG" suffix. Order is always {Er,Fr1,Fr2,Fr3}.
inline std::vector<std::string> GroupFieldNames(const int g) {
  if (g == 0) return {"rad.Er", "rad.Fr1", "rad.Fr2", "rad.Fr3"};
  const std::string s = "_g" + std::to_string(g);
  return {"rad.Er" + s, "rad.Fr1" + s, "rad.Fr2" + s, "rad.Fr3" + s};
}

//! All radiation moment field names across n_group groups (group 0 first). For n_group=1
//! this is exactly the four gray names, in the original order.
inline std::vector<std::string> AllRadFieldNames(const int n_group) {
  std::vector<std::string> v;
  v.reserve(4 * n_group);
  for (int g = 0; g < n_group; ++g)
    for (const auto &nm : GroupFieldNames(g)) v.push_back(nm);
  return v;
}

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
