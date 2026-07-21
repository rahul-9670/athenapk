#ifndef CHEMISTRY_NETWORK_H2_HPP_
#define CHEMISTRY_NETWORK_H2_HPP_
//========================================================================================
// AthenaPK GPU chemistry -- minimal H2 formation/destruction network (greenfield port of
// athena++/src/chemistry/network/H2.cpp). Device-callable; the species are carried as
// AthenaPK passive scalars (abundances x_i = n_i / n_H), evolved operator-split per cell.
//
// Species (NSPEC = 2):  0 = H,  1 = H2   (H-nucleus conservation: x_H + 2 x_H2 = 1)
// Reactions:
//   H2 formation on grains:   rate_gr = kgr * n_H * x_H      [s^-1]
//   H2 destruction by CRs:    rate_cr = kcr * x_H2           [s^-1]
//   d(x_H2)/dt =  rate_gr - rate_cr
//   d(x_H )/dt = -2 rate_gr + 2 rate_cr
// Equilibrium:  kgr n_H x_H = kcr x_H2,  with x_H + 2 x_H2 = 1.
//
// All abundances are dimensionless; rates are evaluated in cgs and converted to the code
// time unit. n_H [cm^-3] = rho_code * rho_unit / (mu_n m_H).
//========================================================================================

// Becomes KOKKOS_INLINE_FUNCTION when compiled inside AthenaPK; plain inline for the
// standalone unit test.
#ifdef KOKKOS_INLINE_FUNCTION
#define CHEM_FN KOKKOS_INLINE_FUNCTION
#else
#define CHEM_FN inline
#endif

#include <cmath>

namespace Chemistry {

constexpr int NSPEC_H2 = 2;
constexpr int iH  = 0;
constexpr int iH2 = 1;

struct H2Network {
  // cgs rate coefficients
  double kgr = 3.0e-17;     // grain H2 formation [cm^3 s^-1]
  double xi_cr = 2.0e-16;   // primary CR ionization rate [s^-1]
  double kcr = 3.0 * 2.0e-16; // CR H2 destruction [s^-1] (= 3 * xi_cr)
  // unit / composition
  double rho_unit = 5.467e-19;  // g/cm^3 per code density
  double t_unit = 1.476998822551e12; // s per code time
  double mu_n = 2.33;           // mean molecular weight of neutral gas
  double m_H = 1.6726219e-24;   // g
  double x_floor = 1.0e-20;     // abundance floor
  int    nsub_max = 200;        // cap on chemical sub-steps
  double cfl = 0.1;             // chemical sub-step CFL

  // number density of H nuclei [cm^-3] from code-unit density
  CHEM_FN double nH(double rho_code) const {
    return rho_code * rho_unit / (mu_n * m_H);
  }

  // RHS in CODE time units: ydot[i] = d x_i / d t_code
  CHEM_FN void rhs(const double *y, double n_H, double *ydot) const {
    const double rate_gr = kgr * n_H * y[iH];   // s^-1
    const double rate_cr = kcr * y[iH2];        // s^-1
    ydot[iH2] = (rate_gr - rate_cr) * t_unit;
    ydot[iH]  = (-2.0 * rate_gr + 2.0 * rate_cr) * t_unit;
  }

  // Analytic equilibrium x_H2 (for validation): solve kgr n_H x_H = kcr x_H2 with
  // x_H = 1 - 2 x_H2  =>  a*x_H2 = kgr n_H (1 - 2 x_H2)  =>  x_H2 = kgr n_H / (kcr + 2 kgr n_H)
  CHEM_FN double x_H2_eq(double rho_code) const {
    const double a = kgr * nH(rho_code);
    return a / (kcr + 2.0 * a);
  }

  // Operator-split per-cell integration over a code-time step dt_code via explicit
  // forward Euler with a chemical-time-step limiter (GPU-tractable stand-in for CVODE).
  // y[] are abundances (modified in place). Returns 1 if the integration was TRUNCATED
  // (nsub_max sub-steps exhausted before reaching dt_code), else 0.
  CHEM_FN int integrate_cell(double *y, double rho_code, double dt_code) const {
    const double n_H = nH(rho_code);
    double t = 0.0;
    int nsub = 0;
    while (t < dt_code && nsub < nsub_max) {
      double ydot[NSPEC_H2];
      rhs(y, n_H, ydot);
      // chemical timestep: cfl * min over species of |y/ydot|
      double dt_chem = dt_code - t;
      for (int s = 0; s < NSPEC_H2; ++s) {
        double yy = (y[s] > x_floor) ? y[s] : x_floor;
        double ad = (ydot[s] < 0 ? -ydot[s] : ydot[s]);
        if (ad > 0.0) {
          double dtc = cfl * yy / ad;
          if (dtc < dt_chem) dt_chem = dtc;
        }
      }
      if (t + dt_chem > dt_code) dt_chem = dt_code - t;
      for (int s = 0; s < NSPEC_H2; ++s) {
        y[s] += ydot[s] * dt_chem;
        if (y[s] < x_floor) y[s] = x_floor;
        if (y[s] > 1.0) y[s] = 1.0;
      }
      t += dt_chem;
      ++nsub;
    }
    return (t < dt_code) ? 1 : 0;
  }
};

} // namespace Chemistry

#endif // CHEMISTRY_NETWORK_H2_HPP_
