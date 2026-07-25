// Flagship Phase 3: production impact of the 2026-07-25 charge-solver fix.
//
// Prints eta_{O,H,A} along the fiducial collapse track (barotropic T(rho), flux-frozen B)
// for BOTH the legacy and the fixed charge-state solvers, so the result-changing region is
// explicit and reproducible. See test_conductivity.py gates 6-8 for the correctness proof
// (the legacy solvers violate charge neutrality by up to ~1e7 and inject spurious thermal
// electrons in cold gas; the fixed pair agrees with an independent Python solve to ~5e-11).
//
// Result (2026-07-25): identical (0.0%) for rho <= 1e-13 -- the whole pre-first-core
// collapse -- and identical again for rho >= 1e-9 where grains have sublimated. The fix
// matters only in the grain-dominated first-core band rho ~ 1e-12..1e-10, where it changes
// eta_O by up to ~375% and FLIPS THE SIGN of eta_H (the Hall sign sets the direction of
// field transport, so this is a physics-level correction, not a tolerance tweak).
//
// Build + run (host g++, no Kokkos/MPI):
//   cd src/hydro/diffusion/tests
//   g++ -O2 -std=c++17 -Ihost_shims impact_charge_solver_fix.cpp -o <out>/impact && <out>/impact
#include <cstdio>
#include <cmath>
#include "../ionization.hpp"
using namespace Ionization;
int main(){
  // Collapse track: barotropic T(rho) -- 10 K isothermal until the first core becomes
  // optically thick (~1e-13 g/cm^3), then adiabatic gamma=5/3.
  const double rho_crit=1e-13, T0=10.0;
  const double B0=4.98e-5*0.15, rho0=5.467e-19;   // flux-frozen B = B0 (rho/rho0)^(1/2)
  std::printf("%10s %8s %11s | %-34s | %-34s\n","rho","T","B",
              "LEGACY  eta_O / eta_H / eta_A","FIXED   eta_O / eta_H / eta_A");
  for (int i=0;i<=10;++i){
    double rho=std::pow(10.0,-18.0+i);
    double T=(rho<rho_crit)?T0:T0*std::pow(rho/rho_crit,2.0/3.0);
    double B=B0*std::sqrt(rho/rho0);
    double eO[2],eH[2],eA[2];
    for(int leg=0;leg<2;++leg){
      IonizationModel m;
      m.rho_unit=1;m.T_unit=1;m.B_unit=1;m.eta_unit=1;
      m.f_dg=0.01; m.zeta=1e-16; m.legacy_charge_solver=(leg==1);
      SetupGrainBins(m,5e-7,2.5e-5,3.5);
      Diffusivities(m,rho,T,B,eO[leg],eH[leg],eA[leg]);
    }
    auto rel=[](double a,double b){return (b!=0.0)?std::fabs(a-b)/std::fabs(b):0.0;};
    std::printf("%10.1e %8.1f %11.3e | %10.3e %10.3e %10.3e | %10.3e %10.3e %10.3e  "
                "dO=%7.1f%% dH=%7.1f%% dA=%7.1f%%\n",
      rho,T,B,eO[1],eH[1],eA[1],eO[0],eH[0],eA[0],
      100*rel(eO[0],eO[1]),100*rel(eH[0],eH[1]),100*rel(eA[0],eA[1]));
  }
  return 0;
}
