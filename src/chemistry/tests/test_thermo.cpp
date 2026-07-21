// WS-2 increment 1: standalone gas-thermochemistry equilibrium curve T_eq(n_H).
// For n_H = 1e2..1e8 cm^-3, evolve the reduced gow17 network to abundance equilibrium
// (self-consistently at T_eq) and solve Gamma - Lambda = 0. Gate: the canonical
// molecular-cloud curve -- ~10-15 K at n=1e4, converging to T_dust for n >~ 1e5.
#include "../network_gow17_reduced.hpp"
#include "../thermo.hpp"
#include <cstdio>
#include <cmath>

using namespace Chemistry;

int main() {
  Gow17ReducedNetwork net;
  ThermoParams tp;
  tp.zeta = net.zeta; // keep CR rate consistent between chemistry and heating
  tp.kgr = net.kgr;
  tp.x_O = 3.2e-4;
  const double T_dust = 10.0; // representative FHC envelope dust temperature [K]

  const double nH_per_code = net.rho_unit / (net.mu_n * net.m_H); // n_H = rho_code * this

  // gate reference points (approx gow17 / canonical cloud): n_H -> T_eq [K]
  struct Ref { double nH, Tref, tol; };
  const Ref refs[] = {{1.0e4, 13.0, 0.30}, {1.0e6, T_dust, 0.30}};

  printf("  n_H[cm^-3]   T_eq[K]   x_H       x_C+       x_CO      A_V     G_eff\n");
  double Teq_at[16];
  double nH_list[16];
  int nN = 0;
  for (double lognH = 2.0; lognH <= 8.001; lognH += 0.5) {
    const double n_H = std::pow(10.0, lognH);
    const double rho_code = n_H / nH_per_code;

    // Self-consistent T_eq + abundance equilibrium (abundances depend weakly on T via
    // recombination): iterate integrate-to-equilibrium <-> SolveTeq a few times.
    double y[NSPEC_GOW] = {1.0e-4, net.xe_floor, 1.0e-6, net.x_floor, net.xe_floor};
    double T = 15.0;
    double x_H = 1.0, x_Cp = 0.0, x_CO = 0.0;
    for (int outer = 0; outer < 4; ++outer) {
      // relax abundances to equilibrium at the current T
      for (int it = 0; it < 40000; ++it) net.integrate_cell(y, rho_code, T, 0.05);
      x_H = 1.0 - 2.0 * y[gH2] - y[gHp];
      if (x_H < 0.0) x_H = 0.0;
      x_Cp = y[gCp];
      x_CO = y[gCO];
      T = SolveTeq(tp, n_H, x_H, x_Cp, x_CO, T_dust);
    }
    const double A_V = tp.av0 * std::pow(n_H / tp.av_n0, tp.av_exp);
    const double G_eff = tp.G0 * std::exp(-2.5 * A_V);
    printf("  %.3e   %7.3f   %.3e  %.3e  %.3e  %5.2f   %.2e\n", n_H, T, x_H, x_Cp, x_CO, A_V,
           G_eff);
    nH_list[nN] = n_H;
    Teq_at[nN] = T;
    ++nN;
  }

  // --- Gate checks ---
  printf("\n--- GATE (T_eq within +-30%% of canonical values) ---\n");
  int fail = 0;
  for (const auto &r : refs) {
    // nearest tabulated n_H
    double best = 1e99, Tbest = 0;
    for (int i = 0; i < nN; ++i)
      if (std::fabs(std::log10(nH_list[i] / r.nH)) < best) {
        best = std::fabs(std::log10(nH_list[i] / r.nH));
        Tbest = Teq_at[i];
      }
    const double rel = std::fabs(Tbest - r.Tref) / r.Tref;
    const bool ok = rel <= r.tol;
    printf("  n_H=%.1e: T_eq=%.2f K  (ref %.1f K, rel %.1f%%)  %s\n", r.nH, Tbest, r.Tref,
           100 * rel, ok ? "PASS" : "FAIL");
    if (!ok) ++fail;
  }
  // monotonic convergence to T_dust at high n
  const double Thi = Teq_at[nN - 1];
  const bool conv = std::fabs(Thi - T_dust) / T_dust < 0.15;
  printf("  n_H=%.1e: T_eq=%.2f K -> T_dust=%.1f K (converged: %s)\n", nH_list[nN - 1], Thi,
         T_dust, conv ? "PASS" : "FAIL");
  if (!conv) ++fail;

  printf("\n%s (%d failures)\n", fail == 0 ? "ALL GATES PASS" : "GATE FAILURES", fail);
  return fail == 0 ? 0 : 1;
}
