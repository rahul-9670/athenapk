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
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

#include <basic_types.hpp>      // parthenon::Real
#include <parthenon_arrays.hpp> // parthenon::ParArray2D (device tables)

#include "opacity_table_format.hpp" // shared table header parser (audit N13)

namespace Radiation {

using parthenon::ParArray2D;
using parthenon::Real;

// Max compile-time group count (POD-friendly; N_GROUP=1 default keeps the gray path).
constexpr int MAX_GROUP = 8;

//----------------------------------------------------------------------------------------
// NOTE (audit N9, 2026-08-05): two device-callable band-moment helpers, PlanckMoment() and
// RossMoment(), used to live here. They were DEAD -- nothing in src/ or the tests ever called
// them. They were superseded by the `phimean` quadrature inside GroupMultsAtT() below, which
// computes the same band means but in the form actually needed (a ratio of a phi-weighted to
// an unweighted moment) and on the host, where the table is built. Removed rather than kept
// "for later": a device-callable helper with a subtle e^{-a}-factored-out return convention
// is exactly the kind of thing a future caller uses without reading the convention. The
// convention itself survives, documented at its one real use site in `phimean`.
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
//! Semenov-class monochromatic dust+gas opacity SHAPE phi(nu,T): the frequency dependence of
//! kappa_nu(rho,T) = kappa_gray(rho,T) * phi(nu,T), where kappa_gray is the existing gray
//! Bell&Lin/dust law (which already carries the rho,T magnitude AND the sublimation drop). phi
//! only redistributes that magnitude across frequency, normalized (by the ratio anchoring in the
//! table build) so the full-spectrum Planck/Rosseland means of kappa_nu equal kappa_gray exactly.
//!
//! Physics captured (protostellar dust, Semenov et al. 2003 / Draine grain emissivity):
//!  - Dust emissivity rises as nu^beta in the far-IR (Q_abs ~ nu^beta, beta~1.5-2 for
//!    silicate/carbon+ice aggregates) and FLATTENS to the geometric (gray) limit above a
//!    grain-size break nu_break (lambda ~ 2*pi*a ~ few micron => nu_break ~ 1e14 Hz).
//!  - Above dust sublimation (~1500 K) grains are destroyed: the opacity becomes gas-dominated
//!    (molecular / H- / Kramers continuum), taken here as ~gray (phi -> 1). survival(T) ramps
//!    the dust frequency structure off across [T_sub_lo, T_sub_hi].
//! beta=0 => phi==1 everywhere => gray-per-group (equivalence gate exact).
struct DustOpacityModel {
  Real beta = 0.0;             // far-IR dust emissivity index (0 => gray)
  Real nu_break_hz = 1.0e14;   // grain-size turnover: Q_abs flattens to gray above this
  Real T_sub_lo = 1400.0;      // dust sublimation ramp start [K]
  Real T_sub_hi = 1600.0;      // dust sublimation ramp end   [K] (phi->1, gray gas)
  Real h_over_k = 4.799243e-11;

  //! Dust survival fraction (1 cold, 0 above sublimation).
  Real survival(const Real Tk) const {
    if (Tk <= T_sub_lo) return 1.0;
    if (Tk >= T_sub_hi) return 0.0;
    return (T_sub_hi - Tk) / (T_sub_hi - T_sub_lo);
  }
  //! Opacity frequency shape at x = h*nu/(k*T), matter temperature Tk. Host-side (table build).
  Real phi(const Real x, const Real Tk) const {
    if (beta == 0.0) return 1.0;
    const Real S = survival(Tk);
    const Real xb = h_over_k * nu_break_hz / Tk;         // break in x-units at this T
    Real psi = std::pow(x / (xb + 1.0e-300), beta);      // nu^beta far-IR rise
    if (psi > 1.0) psi = 1.0;                            // geometric (gray) limit above break
    return S * psi + (1.0 - S);                          // dust shape blended toward gray gas
  }
};

//----------------------------------------------------------------------------------------
//! Ratio-anchored per-group Planck & Rosseland opacity multipliers at matter temperature Tk:
//!   mult_P,g = <phi>_{Planck,band g} / <phi>_{Planck,full}       (=> sum_g mult_P,g*frac_g = 1)
//!   mult_R,g = <1/phi>_{Ross,full}  / <1/phi>_{Ross,band g}      (harmonic; mult_R,full = 1)
//! <.>_Planck weights by B_nu (x^3/(e^x-1)); <.>_Ross by dB/dT (x^4 e^x/(e^x-1)^2). The e^{-a}
//! factor is pulled out of each band ratio (integrand carries e^{a-x}) so deep-Wien bands stay
//! finite. Host-side, called once per T-grid point in the table build.
inline void GroupMultsAtT(const DustOpacityModel &m, const RadGroups &g, const Real Tk,
                          Real *mP, Real *mR) {
  const int NQ = 240;
  const Real X = 80.0;
  // <phi^{inv?}>_weight over [a,b]; weight: 0=Planck (x^3), 1=Rosseland (x^4 * e^x/(e^x-1))
  auto phimean = [&](Real a, Real b, bool ross, bool inv) -> Real {
    if (a < 1.0e-8) a = 1.0e-8;
    if (b <= a) return 1.0;
    const int n = NQ;
    const Real h = (b - a) / n;
    Real num = 0.0, den = 0.0;
    for (int k = 0; k <= n; ++k) {
      const Real x = a + k * h;
      const Real w = (k == 0 || k == n) ? 1.0 : ((k % 2) ? 4.0 : 2.0);
      const Real ex = std::exp(-x);
      Real base = std::exp(a - x); // e^{a-x}, common factor pulled out (cancels in the ratio)
      base *= ross ? (x * x * x * x / ((1.0 - ex) * (1.0 - ex))) : (x * x * x / (1.0 - ex));
      Real ph = m.phi(x, Tk);
      if (inv) ph = 1.0 / ph;
      num += w * ph * base;
      den += w * base;
    }
    return num / (den + 1.0e-300);
  };
  const Real phiP_full = phimean(1.0e-8, X, false, false);
  const Real invphiR_full = phimean(1.0e-8, X, true, true);
  for (int gi = 0; gi < g.n_group; ++gi) {
    Real a = g.h_over_k * g.nu_edge[gi] / Tk;
    Real b = (gi == g.n_group - 1) ? (a + 60.0) : g.h_over_k * g.nu_edge[gi + 1] / Tk;
    if (b > a + 60.0) b = a + 60.0;
    mP[gi] = phimean(a, b, false, false) / (phiP_full + 1.0e-300);
    mR[gi] = invphiR_full / (phimean(a, b, true, true) + 1.0e-300);
  }
}

//----------------------------------------------------------------------------------------
//! Tabulated per-group opacity multipliers vs matter temperature. The band means are a
//! per-T Simpson quadrature of the Semenov-class shape phi -- precomputed here on a log10(T[K])
//! grid once at init and O(1) linear-interpolated on device (the by-value-Views-into-kernel
//! pattern of EosTable). active_=false (beta=0 or gray) => mult=1, exact gray-per-group.
struct GroupOpacityTable {
  ParArray2D<Real> wP_, wR_; // [n_group][nT] Planck- and Rosseland-mean multipliers
  Real lT0_ = 0.0, dlT_ = 1.0; // log10(T[K]) grid origin + spacing
  int ng_ = 1, nT_ = 0;
  bool active_ = false;

  KOKKOS_INLINE_FUNCTION
  Real Lookup(const ParArray2D<Real> &W, const int g, const Real Tk) const {
    if (!active_) return 1.0;
    Real f = (std::log10(Tk) - lT0_) / dlT_;
    int i = static_cast<int>(f);
    // AUDIT 2026-08-05 (N10): clamp HIGH first, then LOW. The old order (low then high) left
    // i = -1 whenever nT_ < 2, because the second clamp `i > nT_-2` then pushed i BELOW the
    // zero it had just been raised to -- an out-of-bounds device read of W(g,-1). nT_ < 2 is
    // rejected at construction now (see BuildGroupOpacityTable*), so this is belt-and-braces,
    // but the clamp order was wrong on its own terms and cost nothing to fix.
    if (i > nT_ - 2) i = nT_ - 2;
    if (i < 0) i = 0;
    Real t = f - i;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    return (1.0 - t) * W(g, i) + t * W(g, i + 1);
  }
  KOKKOS_INLINE_FUNCTION Real PlanckMult(const int g, const Real Tk) const {
    return Lookup(wP_, g, Tk);
  }
  KOKKOS_INLINE_FUNCTION Real RossMult(const int g, const Real Tk) const {
    return Lookup(wR_, g, Tk);
  }
};

//! Build the tabulated multipliers from the Semenov-class band means over [Tmin,Tmax] K,
//! log-spaced with nT points. beta=0 / gray => inactive (mult==1, exact gray equivalence).
inline GroupOpacityTable BuildGroupOpacityTable(const DustOpacityModel &m, const RadGroups &g,
                                                const Real Tmin_K, const Real Tmax_K,
                                                const int nT) {
  GroupOpacityTable tab;
  tab.ng_ = g.n_group;
  tab.nT_ = nT;
  if (g.n_group == 1 || m.beta == 0.0) return tab; // inactive => Lookup returns 1
  // AUDIT 2026-08-05 (N10): validate the T grid before it can reach a device kernel.
  // radiation/opacity_table_nT, _Tmin_K, _Tmax_K are free deck knobs with no checks: nT < 2
  // made dlT_ divide by zero AND drove Lookup's index negative; Tmin <= 0 makes lT0_ = -inf
  // (every lookup then NaN); Tmax <= Tmin makes dlT_ <= 0 so the interpolation runs backwards.
  // All three are silent on the host and only show up as garbage opacity multipliers.
  if (nT < 2)
    throw std::runtime_error("GroupOpacityTable: radiation/opacity_table_nT must be >= 2 (got " +
                             std::to_string(nT) + ")");
  if (!(Tmin_K > 0.0) || !(Tmax_K > Tmin_K))
    throw std::runtime_error(
        "GroupOpacityTable: need 0 < radiation/opacity_table_Tmin_K < opacity_table_Tmax_K (got " +
        std::to_string(Tmin_K) + ", " + std::to_string(Tmax_K) + ")");
  tab.active_ = true;
  tab.lT0_ = std::log10(Tmin_K);
  tab.dlT_ = (std::log10(Tmax_K) - tab.lT0_) / (nT - 1);
  tab.wP_ = ParArray2D<Real>("rad_wP", g.n_group, nT);
  tab.wR_ = ParArray2D<Real>("rad_wR", g.n_group, nT);
  auto hP = Kokkos::create_mirror_view(tab.wP_);
  auto hR = Kokkos::create_mirror_view(tab.wR_);
  Real mP[MAX_GROUP], mR[MAX_GROUP];
  for (int it = 0; it < nT; ++it) {
    const Real Tk = std::pow(10.0, tab.lT0_ + it * tab.dlT_);
    GroupMultsAtT(m, g, Tk, mP, mR);
    for (int ig = 0; ig < g.n_group; ++ig) {
      hP(ig, it) = mP[ig];
      hR(ig, it) = mR[ig];
    }
  }
  Kokkos::deep_copy(tab.wP_, hP);
  Kokkos::deep_copy(tab.wR_, hR);
  return tab;
}

//! Build the per-group multiplier table by reading the TABULATED opacity file (produced by
//! gen_opacity_table.py). Binary layout: int64 [ng,nr,nT]; double [lr0,dlr,lT0,dlT]; float64
//! gray kP[nr,nT], kR, ks (skipped here); then multipliers mP[ng,nT], mR[ng,nT] (T-only,
//! rho-independent). This REPLACES the analytic DustOpacityModel band means with the tabulated
//! frequency-resolved values. n_group must match the file's ng.
//!
//! The multipliers are RATIOS of opacities, so they are dimensionless and carry no unit
//! error in either format version -- this reader only needs the header to know where the
//! payload starts. It parses it through the SAME helper as OpacityTable::Load (audit N13);
//! the two used to carry independent copies of the layout.
inline GroupOpacityTable BuildGroupOpacityTableFromFile(const std::string &path,
                                                        const int n_group) {
  GroupOpacityTable tab;
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("GroupOpacityTable: cannot open " + path);
  const auto hdr = ReadOpacityHeader(f, path);
  const int ng = hdr.ng;
  const int nr = hdr.nr;
  const int nT = hdr.nT;
  if (ng != n_group)
    throw std::runtime_error("GroupOpacityTable: file n_group mismatch (" +
                             std::to_string(ng) + " != " + std::to_string(n_group) + ")");
  tab.ng_ = ng;
  tab.nT_ = nT;
  tab.lT0_ = hdr.lT0;
  tab.dlT_ = hdr.dlT;
  if (ng == 1) return tab; // gray => inactive (Lookup returns 1)
  // audit N10: the T grid comes straight out of a binary header here, so it is even less
  // checked than the analytic path. Same three conditions, same reasons.
  if (nT < 2 || !(tab.dlT_ > 0.0))
    throw std::runtime_error("GroupOpacityTable: bad T grid in " + path + " (nT=" +
                             std::to_string(nT) + ", dlogT=" + std::to_string(tab.dlT_) +
                             "); need nT >= 2 and dlogT > 0");
  // skip the 3 gray arrays kP,kR,ks (each nr*nT doubles)
  f.seekg(static_cast<std::streamoff>(3) * nr * nT * sizeof(double), std::ios::cur);
  tab.active_ = true;
  tab.wP_ = ParArray2D<Real>("rad_wP", ng, nT);
  tab.wR_ = ParArray2D<Real>("rad_wR", ng, nT);
  auto hP = Kokkos::create_mirror_view(tab.wP_);
  auto hR = Kokkos::create_mirror_view(tab.wR_);
  std::vector<double> buf(static_cast<size_t>(ng) * nT);
  auto read_mult = [&](decltype(hP) &h) {
    f.read(reinterpret_cast<char *>(buf.data()),
           static_cast<std::streamsize>(buf.size() * sizeof(double)));
    for (int ig = 0; ig < ng; ++ig)
      for (int it = 0; it < nT; ++it) h(ig, it) = static_cast<Real>(buf[ig * nT + it]);
  };
  read_mult(hP); // mP
  read_mult(hR); // mR
  if (!f) throw std::runtime_error("GroupOpacityTable: truncated file " + path);
  Kokkos::deep_copy(tab.wP_, hP);
  Kokkos::deep_copy(tab.wR_, hR);
  return tab;
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
  // There are n_group-1 interior edges, so the log-spacing denominator is n_group-2.
  // AUDIT 2026-08-05 (N11): n_group == 2 has exactly ONE interior edge, so it cannot span
  // [nu_min, nu_max] -- the loop places that edge at nu_min and nu_max is silently DROPPED.
  // (The `+1e-30` below existed only to keep that degenerate case from dividing 0 by 0; it
  // gave the right answer for the wrong-looking reason.) Say so out loud instead: a deck that
  // sets nu_max and gets two groups deserves to know its upper edge was ignored.
  const Real ll = std::log10(nu_min_hz), lh = std::log10(nu_max_hz);
  if (n_group == 2) {
    g.nu_edge[1] = nu_min_hz; // the single split point
    g.nu_edge[2] = 1.0e30;
    return g;
  }
  for (int i = 1; i < n_group; ++i)
    g.nu_edge[i] = std::pow(10.0, ll + (lh - ll) * (i - 1) / (n_group - 2));
  g.nu_edge[n_group] = 1.0e30;
  return g;
}

} // namespace Radiation

#endif // RADIATION_RADIATION_GROUPS_HPP_
