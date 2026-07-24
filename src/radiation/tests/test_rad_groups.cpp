// Flagship Phase 4 (multigroup RHD scaffold): gate for the frequency-group structure.
//
// Compiles the REAL src/radiation/radiation_groups.hpp host-only (via host_shims/) and checks
// the group-integrated Planck machinery that later multigroup increments build on.
//
// Build + run:  cd src/radiation/tests
//   g++ -O2 -std=c++17 -Ihost_shims test_rad_groups.cpp -o <scratch>/t && <scratch>/t
//
// Gates:
//   1. Planck cumulative fraction series F(x) matches an INDEPENDENT numerical integral of
//      (15/pi^4) x^3/(e^x-1); F(0)=0, F(inf)=1, monotonic.
//   2. group-fraction conservation: sum_g PlanckFraction(g,T) == 1 for any multigroup structure
//      spanning [0,inf), at every T (radiation energy is partitioned, not lost).
//   3. gray reduction: n_group=1 -> PlanckFraction(0,T) == 1 exactly (bit-identical to gray).
#include <cmath>
#include <cstdio>
#include <initializer_list>

#include "../radiation_groups.hpp"

using namespace Radiation;

static int fails = 0;
static void rep(const char *n, bool ok, double v, double r) {
  std::printf("  %-50s %s (got %.6e, ref %.6e)\n", n, ok ? "PASS" : "FAIL", v, r);
  if (!ok) ++fails;
}

// Independent numerical F(x) = (15/pi^4) Int_0^x u^3/(e^u-1) du by fine trapezoid.
static double F_numeric(double x) {
  if (x <= 0) return 0.0;
  const int N = 200000;
  const double h = x / N;
  double s = 0.0;
  for (int i = 0; i <= N; ++i) {
    const double u = i * h;
    double f;
    if (u < 1e-8) f = u * u; // u^3/(e^u-1) -> u^2 as u->0
    else f = u * u * u / (std::exp(u) - 1.0);
    s += (i == 0 || i == N) ? 0.5 * f : f;
  }
  return 15.0 / (M_PI * M_PI * M_PI * M_PI) * s * h;
}

int main() {
  std::printf("--- GATE 1: Planck cumulative-fraction series vs numerical integral ---\n");
  double worst = 0.0;
  for (double x : {0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0}) {
    const double se = PlanckCumFraction(x), nu = F_numeric(x);
    worst = std::max(worst, std::fabs(se - nu));
  }
  rep("max |series - numeric| over x in [0.5,20]", worst < 1e-6, worst, 0.0);
  rep("F(0)=0", PlanckCumFraction(0.0) == 0.0, PlanckCumFraction(0.0), 0.0);
  rep("F(60)=1 (all energy)", std::fabs(PlanckCumFraction(60.0) - 1.0) < 1e-9,
      PlanckCumFraction(60.0), 1.0);
  bool mono = true;
  double prev = 0.0;
  for (double x = 0.1; x <= 40; x += 0.1) { double f = PlanckCumFraction(x); if (f < prev - 1e-12) mono = false; prev = f; }
  rep("F(x) monotonic non-decreasing", mono, 1.0, 1.0);

  std::printf("--- GATE 2: sum_g PlanckFraction(g,T) == 1 (energy partition) ---\n");
  // A 5-group structure over IR..UV; edges [0, 1e12, 3e13, 1e14, 1e15, inf) Hz.
  RadGroups g5 = BuildRadGroups(5, 1e12, 1e15);
  double worst_sum = 0.0;
  for (double T : {10.0, 30.0, 100.0, 300.0, 1000.0, 5000.0}) {
    double sum = 0.0;
    for (int gi = 0; gi < g5.n_group; ++gi) sum += g5.PlanckFraction(gi, T);
    worst_sum = std::max(worst_sum, std::fabs(sum - 1.0));
  }
  rep("max |sum_g frac - 1| over T", worst_sum < 1e-6, worst_sum, 0.0);

  std::printf("--- GATE 3: gray reduction (n_group=1 -> fraction 1) ---\n");
  RadGroups g1 = BuildRadGroups(1, 0, 0);
  bool gray_ok = true;
  for (double T : {10.0, 100.0, 1000.0}) gray_ok = gray_ok && (g1.PlanckFraction(0, T) == 1.0);
  rep("n_group=1: PlanckFraction(0,T) == 1 exactly", gray_ok, g1.PlanckFraction(0, 100.0), 1.0);

  std::printf("\n%s (%d failures)\n", fails == 0 ? "ALL GATES PASS" : "GATE FAILURES", fails);
  return fails ? 1 : 0;
}
