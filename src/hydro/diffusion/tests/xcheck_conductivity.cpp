// Flagship Phase 3: numerical C++-vs-Python cross-check of the Wardle conductivity tensor.
//
// Compiles the REAL src/hydro/diffusion/ionization.hpp host-only (its compute path is POD +
// stack arrays, only KOKKOS_INLINE_FUNCTION + <cmath>) via the tiny build shims in
// host_shims/ (which stub <Kokkos_Core.hpp> and <basic_types.hpp>), and dumps the charge
// state (n_e, n_i, per-bin grain charge Z_k) AND eta_{O,H,A} in CGS at sample (rho,T,B).
//
// Two modes, selected by argv[1] = dust-to-gas ratio f_dg:
//   f_dg = 0    -> gas-phase only; matches the original independent Python reference
//                  (eta_O,eta_H to ~1e-7; eta_A to ~2e-3, cancellation-limited).
//   f_dg > 0    -> GRAIN-INCLUSIVE (default 0.01), the Phase-3 remaining item: exercises
//                  SolveCharges' MRN grain-charge fixed point + the charged-grain terms in
//                  the tensor sum. test_conductivity.py reimplements that model independently
//                  (different algorithm: 1-D bisection in r = n_e/n_i, vs the C++ relaxed
//                  fixed point with inner Newton) and compares against this dump.
//
// Build + run (host g++, no Kokkos/MPI):
//   cd src/hydro/diffusion/tests
//   g++ -O2 -std=c++17 -Ihost_shims xcheck_conductivity.cpp -o <out>/xcheck
//   <out>/xcheck 0.01 > <out>/xcheck_grains.csv
#include <cstdio>
#include <cstdlib>

#include "../ionization.hpp"

using namespace Ionization;

int main(int argc, char **argv) {
  const double f_dg = (argc > 1) ? std::atof(argv[1]) : 0.01;

  IonizationModel m;
  m.rho_unit = 1.0; m.T_unit = 1.0; m.B_unit = 1.0; m.eta_unit = 1.0; // evaluate in CGS
  m.f_dg = f_dg;
  m.zeta = 1.0e-16;                      // production CR rate (matches the Python ref)
  m.ad_closure = ADClosure::tensor;      // compare the TENSOR eta_A (default is single_fluid)
  SetupGrainBins(m, 5.0e-7, 2.5e-5, 3.5); // MRN bins over 0.005-0.25 micron

  // Sample grid. Includes T above T_subl=1500 K so the grain-sublimation branch (ng->0,
  // thermal Saha takes over) is exercised on both sides.
  const double rhos[] = {1e-16, 1e-14, 1e-12, 1e-10, 1e-8};
  const double Ts[] = {15.0, 50.0, 300.0, 2000.0};
  const double Bs[] = {1e-5, 1e-3};

  // header: the a_k/mw_k bin setup, so Python can use the IDENTICAL MRN discretization
  // (the bin quadrature is a modelling choice, not the physics under test).
  std::printf("# N_BIN=%d\n", N_BIN);
  for (int k = 0; k < N_BIN; ++k)
    std::printf("# bin %d a_k=%.12e mw_k=%.12e\n", k, m.a_k[k], m.mw_k[k]);
  std::printf("# f_dg=%.6e zeta=%.6e T_subl=%.6e\n", m.f_dg, m.zeta, m.T_subl);

  std::printf("rho,T,B,n_e,n_i");
  for (int k = 0; k < N_BIN; ++k) std::printf(",Z%d", k);
  for (int k = 0; k < N_BIN; ++k) std::printf(",ng%d", k);
  std::printf(",eta_O,eta_H,eta_A\n");

  for (double rho : rhos)
    for (double T : Ts)
      for (double B : Bs) {
        // charge state (same call Diffusivities makes internally)
        const double n_n = rho / (m.mu_n * cgs::m_H);
        double n_e, n_i, Zk[N_BIN], ng[N_BIN];
        SolveCharges(m, n_n, T, n_e, n_i, Zk, ng);
        double eO, eH, eA;
        Diffusivities(m, rho, T, B, eO, eH, eA);
        std::printf("%.6e,%.6e,%.6e,%.10e,%.10e", rho, T, B, n_e, n_i);
        for (int k = 0; k < N_BIN; ++k) std::printf(",%.10e", Zk[k]);
        for (int k = 0; k < N_BIN; ++k) std::printf(",%.10e", ng[k]);
        std::printf(",%.10e,%.10e,%.10e\n", eO, eH, eA);
      }
  return 0;
}
