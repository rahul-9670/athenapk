#include "network_h2.hpp"
#include <cstdio>
using namespace Chemistry;
int main(){
  H2Network net;
  // test at several code densities (rho_code; rho_cgs = rho_code*rho_unit)
  double rhos[] = {1.0, 1e2, 1e4, 1e6};
  for(double rho: rhos){
    // start fully atomic
    double y[2] = {1.0, net.x_floor};
    double nH = net.nH(rho);
    // integrate many code-time steps to reach equilibrium (dt ~ a fraction of t_unit)
    double dt = 0.05;
    for(int it=0; it<20000; ++it) net.integrate_cell(y, rho, dt);
    double cons = y[iH] + 2.0*y[iH2];      // should be ~1
    double xeq = net.x_H2_eq(rho);
    double err = (y[iH2]>0)? (y[iH2]-xeq)/xeq : 1.0;
    printf("rho_code=%-8g  nH=%.3e cm^-3  x_H2=%.6e  x_H2_eq=%.6e  rel.err=%+.2e  (x_H+2x_H2=%.6f)\n",
           rho, nH, y[iH2], xeq, err, cons);
  }
  return 0;
}
