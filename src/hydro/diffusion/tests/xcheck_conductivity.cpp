// Flagship Phase 3: numerical C++-vs-Python cross-check of the Wardle conductivity tensor.
//
// Compiles the REAL src/hydro/diffusion/ionization.hpp host-only (its compute path is POD +
// stack arrays, only KOKKOS_INLINE_FUNCTION + <cmath>) via the tiny build shims in
// host_shims/ (which stub <Kokkos_Core.hpp> and <basic_types.hpp>), and dumps eta_{O,H,A}
// in CGS at sample (rho,T,B). test_conductivity.py's companion comparison confirms the shipped
// device tensor matches the INDEPENDENT Python reimplementation to ~1e-7 (eta_O,eta_H); eta_A
// agrees to ~2e-3, limited by the catastrophic cancellation in its tensor closure
// (eta_A = c^2/4pi sigma_P/sigma_perp^2 - eta_O, with eta_A << eta_O where the worst point sits)
// -- still < the Phase-3 <1% criterion. Grains are OFF here (f_dg=0) so this matches the
// gas-phase Python reference; the grain-inclusive cross-check needs SolveCharges' grain-charge
// model reimplemented in Python (remaining item).
//
// Build + run (host g++, no Kokkos/MPI):
//   cd src/hydro/diffusion/tests
//   g++ -O2 -std=c++17 -Ihost_shims xcheck_conductivity.cpp -o /tmp/xcheck && /tmp/xcheck
#include <cstdio>

#include "../ionization.hpp"

using namespace Ionization;

int main() {
  IonizationModel m;
  m.rho_unit = 1.0; m.T_unit = 1.0; m.B_unit = 1.0; m.eta_unit = 1.0; // evaluate in CGS
  m.f_dg = 0.0;                          // grains OFF -> gas-phase (matches the Python ref)
  m.zeta = 1.0e-16;                      // production CR rate (matches the Python ref)
  m.ad_closure = ADClosure::tensor;      // compare the TENSOR eta_A (default is single_fluid)
  SetupGrainBins(m, 5.0e-7, 2.5e-5, 3.5); // fill a_k so m_g != 0 (grains still off via f_dg=0)

  const double rhos[] = {1e-16, 1e-14, 1e-12, 1e-10, 1e-8};
  const double Ts[] = {15.0, 50.0};
  const double Bs[] = {1e-5, 1e-3};
  std::printf("rho,T,B,eta_O,eta_H,eta_A\n");
  for (double rho : rhos)
    for (double T : Ts)
      for (double B : Bs) {
        double eO, eH, eA;
        Diffusivities(m, rho, T, B, eO, eH, eA);
        std::printf("%.6e,%.6e,%.6e,%.10e,%.10e,%.10e\n", rho, T, B, eO, eH, eA);
      }
  return 0;
}
