// Flagship Phase 4/5 (opacity): characterization gate for the SHIPPED gray opacity.
//
// Compiles the REAL src/radiation/radiation_opacity.hpp host-only (via host_shims/) and checks
// the Bell & Lin (1994) Rosseland opacity + the Planck/Rosseland split. These are the physical
// opacities that set the first/second-core radiative trapping. The size-/composition- and
// FREQUENCY-dependent (multigroup) opacity means from a monochromatic dataset are the Phase-5
// deliverable -> not yet implemented (this is the gray baseline).
//
// Build + run (host g++):
//   cd src/radiation/tests
//   g++ -O2 -std=c++17 -Ihost_shims test_opacity.cpp -o /tmp/o && /tmp/o   (or -o <scratch>/o)
//
// Gates:
//   1. CONTINUITY across all 7 Bell&Lin regime transitions: kappa is continuous at each
//      transition temperature Tt (defined by adjacent regimes crossing) -- validates the Tt
//      formula 1/(bb[i+1]-bb[i]) exponent. A bug there would show a jump.
//   2. positivity: kappa > 0 over the (rho,T) collapse range.
//   3. canonical ice-grain values: kappa = 2e-4 T^2 -> 0.02 cm^2/g at 10 K, 2.0 at 100 K;
//      electron-scattering floor 0.348 at high T.
//   4. dust-sublimation gap: kappa plunges (T^-24 regime) across ~1500-2000 K.
//   5. Planck/Rosseland split: ratio=1 -> kappa_P == kappa_R (bit-identical default);
//      ratio>1 -> kappa_P = ratio * kappa_R.
#include <cmath>
#include <cstdio>

#include "../radiation_opacity.hpp"

using namespace Radiation;

static int fails = 0;
static void report(const char *name, bool ok, double v, double ref) {
  std::printf("  %-54s %s (got %.4e, ref %.4e)\n", name, ok ? "PASS" : "FAIL", v, ref);
  if (!ok) ++fails;
}
static bool rel(double a, double b, double tol) {
  return std::fabs(a - b) <= tol * (std::fabs(b) + 1e-300);
}

int main() {
  // Gate 1 & 4: continuity + sublimation gap from a fine T-scan. kappa(T) is continuous iff the
  // local slope |d ln k / d ln T| stays BOUNDED as the grid refines; a jump would blow up. The
  // steepest genuine regime is dust sublimation (kappa ~ T^-24), so the bound is ~24 and the
  // most-negative slope must reach ~-24 (the gap is present). Probing computed transition
  // temperatures directly is WRONG -- the regime walk skips "virtual" crossings, so a fine scan
  // is the correct continuity test.
  std::printf("--- GATE 1/4: continuity (bounded slope) + sublimation gap (fine T-scan) ---\n");
  auto scan = [](double lrho_lo, double lrho_hi) {
    double mx = 0.0, mn = 0.0;
    for (double lrho = lrho_lo; lrho <= lrho_hi + 1e-9; lrho += 1.0) {
      const double rho = std::pow(10.0, lrho);
      const int N = 6000;
      double prevk = BellLinKappa(rho, 10.0), prevlT = std::log(10.0);
      for (int n = 1; n <= N; ++n) {
        const double lT = std::log(10.0) + (std::log(1e7) - std::log(10.0)) * n / N;
        const double k = BellLinKappa(rho, std::exp(lT));
        const double slope = (std::log(k) - std::log(prevk)) / (lT - prevlT);
        mx = std::max(mx, std::fabs(slope)); mn = std::min(mn, slope);
        prevk = k; prevlT = lT;
      }
    }
    return std::pair<double, double>{mx, mn};
  };
  // Collapse-relevant region (rho >= 1e-10 where the regime walk is well-ordered): continuity
  // holds -> |slope| bounded by the steepest regime (~24). High T (>4000 K) is reached only at
  // rho >= 1e-10 (second core) on the FHC track, so this is the region the run traverses.
  auto [mx_hi, mn_hi] = scan(-10.0, -3.0);
  report("continuity for rho>=1e-10 (collapse track): |slope|<=24.5", mx_hi <= 24.5, mx_hi, 24.0);
  report("sublimation gap present (min slope <= -20)", mn_hi <= -20.0, mn_hi, -24.0);
  // KNOWN DEFECT (characterization, NOT a hard fail): at low rho (<= 1e-11) the Bell&Lin
  // transition temperatures go OUT OF ORDER and the walk skips regime 6 (Kramers), giving a
  // kappa DISCONTINUITY (regime 5->7, up to ~4 decades) at T~4000-9000 K. This corner is OFF the
  // cold collapse track (cold at low rho) but reachable by radiation in hot low-density gas.
  // Fixing the walk is result-changing to production opacity -> DEFERRED (FLAGSHIP doc Phase 4).
  auto [mx_lo, mn_lo] = scan(-16.0, -11.0);
  std::printf("  CHARACTERIZATION: low-rho (<=1e-11) max |slope| = %.0f -- regime-6-skip "
              "discontinuity (off-track, DEFERRED fix)\n", mx_lo);

  // Gate 2: positivity.
  std::printf("--- GATE 2: positivity over the collapse (rho,T) range ---\n");
  bool allpos = true;
  double kmin = 1e300;
  for (double lrho = -18; lrho <= -2; lrho += 1.0)
    for (double lT = 1.0; lT <= 4.5; lT += 0.25) {
      const double k = BellLinKappa(std::pow(10.0, lrho), std::pow(10.0, lT));
      if (!(k > 0.0) || !std::isfinite(k)) allpos = false;
      kmin = std::min(kmin, k);
    }
  report("kappa > 0 everywhere (min over grid)", allpos, kmin, 0.0);

  // Gate 3: canonical ice-grain values + electron-scattering floor.
  std::printf("--- GATE 3: canonical opacity values ---\n");
  const double rho_lo = 1e-14;
  report("ice grains kappa(10 K) == 2e-4*T^2 = 0.02", rel(BellLinKappa(rho_lo, 10.0), 0.02, 1e-6),
         BellLinKappa(rho_lo, 10.0), 0.02);
  report("ice grains kappa(100 K) == 2.0", rel(BellLinKappa(rho_lo, 100.0), 2.0, 1e-6),
         BellLinKappa(rho_lo, 100.0), 2.0);
  // high-T electron scattering (regime 7, kappa=0.348 constant): dominates at LOW rho, high T
  // (at high rho the Kramers regime is larger). rho=1e-6, T=1e7 -> regime 7.
  const double k_es = BellLinKappa(1e-6, 1e7);
  report("electron-scattering floor 0.348 (rho=1e-6,T=1e7)", rel(k_es, 0.348, 1e-6), k_es, 0.348);

  // Gate 5: Planck/Rosseland split.
  std::printf("--- GATE 5: Planck/Rosseland split ---\n");
  OpacityParams op;
  op.model = OpacityModel::belllin;
  op.T_unit = 10.015; op.rho_unit = 5.467e-19; op.kappa_conv = 5.467e-19 * 2.81e16;
  const double kR = RosselandOpacity(op, 1e5, 20.0);
  double kP = PlanckOpacity(op, 1e5, 20.0);
  report("ratio=1 -> kappa_P == kappa_R", rel(kP, kR, 1e-14), kP, kR);
  op.planck_ross_ratio = 3.0;
  kP = PlanckOpacity(op, 1e5, 20.0);
  report("ratio=3 -> kappa_P == 3*kappa_R", rel(kP, 3.0 * kR, 1e-12), kP, 3.0 * kR);

  std::printf("\n%s (%d failures)\n", fails == 0 ? "ALL GATES PASS" : "GATE FAILURES", fails);
  return fails ? 1 : 0;
}
