// Helper for the WS-2 increment-2 gate: given (n_H, x_H, x_Cp, x_CO, T_dust) on argv,
// print T_eq from the SAME thermo.hpp the simulation uses. Lets the single-cell run be
// compared to its own self-consistent equilibrium without re-deriving the physics.
#include "../thermo.hpp"
#include <cstdio>
#include <cstdlib>
using namespace Chemistry;
int main(int argc, char **argv) {
  if (argc < 6) {
    fprintf(stderr, "usage: teq_checker n_H x_H x_Cp x_CO T_dust\n");
    return 2;
  }
  ThermoParams tp; // defaults must match the <chemistry> input block
  tp.zeta = 1.0e-16; // single_cell.in sets <chemistry> zeta_cr_cgs = 1e-16 (10x the struct default)
  const double n_H = atof(argv[1]);
  const double x_H = atof(argv[2]);
  const double x_Cp = atof(argv[3]);
  const double x_CO = atof(argv[4]);
  const double T_dust = atof(argv[5]);
  printf("%.6f\n", SolveTeq(tp, n_H, x_H, x_Cp, x_CO, T_dust));
  return 0;
}
