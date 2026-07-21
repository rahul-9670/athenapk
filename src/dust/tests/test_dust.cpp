// WS-4 increment 1: single-cell dust-growth test vs the analytic monodisperse solution.
//   (a) Brownian-dominated (alpha_turb=0): v_B ~ a^{-3/2} => a^{5/2} = a0^{5/2} + (5/2)K t.
//   (b) Turbulent-dominated (St ~ a): da/dt = C a => a(t) = a0 exp(C t) [3 decades].
//   (c) Sublimation switch round-trips f_dg (hot -> ~0, cool -> f_dg_ref) without mass error.
// Gate: model integrate_cell matches the analytic solution to < 5%.
#include "../dust.hpp"
#include <cstdio>
#include <cmath>
using namespace Dust;

int main() {
  int fail = 0;
  const double rho_gas = 1.0e-13; // g/cm^3 (~ n_H 1e11? just a fixed test state)
  const double T = 20.0;          // K
  const double f_dg = 0.01;

  // ---- (a) Brownian regime: alpha_turb = 0 -> pure a^{-3/2} growth ----
  {
    DustModel d;
    d.alpha_turb = 0.0; // kill turbulent term
    d.growth_cap = 0.05;
    d.nsub_max = 100000;
    // analytic K: da/dt = (f_dg rho/(4 rho_gr)) * V0 * a^{-3/2}, V0 = sqrt(12 kB T/(pi^2 rho_gr))
    const double V0 = std::sqrt(12.0 * dcgs::kB * T / (dcgs::pi * dcgs::pi * d.rho_grain));
    const double K = (f_dg * rho_gas / (4.0 * d.rho_grain)) * V0;
    double a = 1.0e-7, fdg = f_dg;
    const double a0 = a;
    const double tmax = 3.0e12; // s
    d.integrate_cell(a, fdg, rho_gas, T, tmax);
    const double a_ana = std::pow(std::pow(a0, 2.5) + 2.5 * K * tmax, 0.4);
    const double rel = std::fabs(a - a_ana) / a_ana;
    printf("(a) Brownian: a_num=%.4e a_ana=%.4e cm  rel=%.2f%%  %s\n", a, a_ana, 100 * rel,
           rel < 0.05 ? "PASS" : "FAIL");
    if (rel >= 0.05) ++fail;
  }

  // ---- (b) Turbulent regime: large omega_ref + large a so St-driven exp growth dominates ----
  {
    DustModel d;
    d.alpha_turb = 0.1;
    d.omega_ref = 1.0e-8; // boost v_turb so St term dominates over Brownian across the range
    d.growth_cap = 0.02;
    d.nsub_max = 100000;
    const double cs = std::sqrt(dcgs::kB * T / (d.mu_gas * d.m_H));
    const double v_th = std::sqrt(8.0 / dcgs::pi) * cs;
    const double C = (f_dg / 4.0) * std::sqrt(d.alpha_turb) * (cs / v_th) * d.omega_ref;
    double a = 1.0e-4, fdg = f_dg;
    const double a0 = a;
    const double tmax = std::log(1000.0) / C; // grow ~3 decades
    d.integrate_cell(a, fdg, rho_gas, T, tmax);
    const double a_ana = a0 * std::exp(C * tmax);
    const double rel = std::fabs(a - a_ana) / a_ana;
    // check v_turb truly dominates at the start (else the exp form is contaminated)
    const double m_grain = (4.0 / 3.0) * dcgs::pi * a0 * a0 * a0 * d.rho_grain;
    const double vB = std::sqrt(8.0 * dcgs::kB * T / (dcgs::pi * 0.5 * m_grain));
    const double St = (d.rho_grain * a0 / (rho_gas * v_th)) * d.omega_ref;
    const double vt = std::sqrt(d.alpha_turb) * cs * St;
    printf("(b) Turbulent: a_num=%.4e a_ana=%.4e cm (%.1f decades)  rel=%.2f%%  vt/vB@start=%.1f  %s\n",
           a, a_ana, std::log10(a / a0), 100 * rel, vt / vB, rel < 0.05 ? "PASS" : "FAIL");
    if (rel >= 0.05) ++fail;
  }

  // ---- (c) Sublimation round-trip: hot destroys f_dg, cool restores; no spurious mass ----
  {
    DustModel d;
    double a = 1.0e-5, fdg = 0.01;
    d.integrate_cell(a, fdg, rho_gas, 2.0e3, 1.0e11); // T > T_subl -> f_dg -> ~0
    const double fdg_hot = fdg;
    d.integrate_cell(a, fdg, rho_gas, 10.0, 1.0e11); // T < T_subl -> restored
    const double fdg_cool = fdg;
    const bool ok = (fdg_hot < 1.0e-3) && (std::fabs(fdg_cool - 0.01) / 0.01 < 0.02);
    printf("(c) Sublimation: f_dg hot=%.3e (->0) cool=%.4f (->0.01)  %s\n", fdg_hot, fdg_cool,
           ok ? "PASS" : "FAIL");
    if (!ok) ++fail;
  }

  // ---- (d) Consumer factor: bit-identical at reference, correct geometric scaling ----
  {
    const double f_ref = 0.01, a_ref = 1.0e-5;
    const double k1 = DustFactor(f_ref, a_ref, f_ref, a_ref);      // reference -> 1
    const double k2 = DustFactor(f_ref, 2.0e-5, f_ref, a_ref);     // 2x grain -> 0.5
    const double k3 = DustFactor(0.02, a_ref, f_ref, a_ref);       // 2x dust -> 2
    const bool ok = std::fabs(k1 - 1.0) < 1e-12 && std::fabs(k2 - 0.5) < 1e-12 &&
                    std::fabs(k3 - 2.0) < 1e-12;
    printf("(d) DustFactor: ref=%.6f (=1) bigger-grain=%.3f (=0.5) more-dust=%.3f (=2)  %s\n",
           k1, k2, k3, ok ? "PASS" : "FAIL");
    if (!ok) ++fail;
  }

  printf("\n%s (%d failures)\n", fail == 0 ? "ALL GATES PASS" : "GATE FAILURES", fail);
  return fail == 0 ? 0 : 1;
}
