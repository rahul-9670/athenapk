#ifndef CHEMISTRY_NETWORK_GOW17_REDUCED_HPP_
#define CHEMISTRY_NETWORK_GOW17_REDUCED_HPP_
//========================================================================================
// AthenaPK GPU chemistry -- reduced H-C-O ionization network (a compact, device-callable
// reduction of Gong, Ostriker & Wolfire 2017, ApJ 843, 38, the network carried by
// athena++/src/chemistry/network/gow17.cpp). The purpose here is NOT the full thermal-
// chemistry of gow17 but a minimal set that evolves the FREE-ELECTRON ABUNDANCE x_e,
// which sets the non-ideal MHD diffusivities (ambipolar / Hall / Ohmic) self-consistently
// and time-dependently -- replacing the equilibrium x_e of ionization.hpp when chemistry
// is active.
//
// Species (NSPEC = 5, abundances x_i = n_i / n_H):
//     0 = H2     1 = H+ (Hp)    2 = C+ (Cp)    3 = CO     4 = e-
// Derived (not stored): x_H  = 1 - 2 x_H2 - x_Hp        (H-nucleus conservation)
//                       x_C0 = x_Ctot - x_Cp - x_CO     (neutral C; >= 0)
//                       x_e  = x_Hp + x_Cp              (charge neutrality; slaved)
//
// Reactions (cgs rate coefficients; UMIST/Gong+2017 representative values):
//   R1 H2 formation on grains      kgr  n_H x_H            -> +x_H2 , -2 x_H
//   R2 H2 CR dissociation/ioniz.   k_cr,H2 x_H2            -> -x_H2 , +2 x_H
//   R3 H  CR ionization            zeta_H x_H             -> +x_Hp , +x_e
//   R4 H+ recombination (gas+grain)(a_rr(T) n_e + k_gr n_H) x_Hp  -> -x_Hp, -x_e
//   R5 C  ionization (CR/photo)    zeta_C x_C0            -> +x_Cp , +x_e
//   R6 C+ recombination            a_C(T) n_e x_Cp        -> -x_Cp , -x_e
//   R7 CO formation (C+ chain)     k_CO n_H x_Cp          -> +x_CO , -x_Cp
//   R8 CO CR destruction           zeta_CO x_CO           -> -x_CO , +x_Cp
//
// Conservation (preserved by the RHS to round-off): H nuclei x_H + 2 x_H2 + x_Hp = 1;
// C nuclei x_Cp + x_CO + x_C0 = x_Ctot. Low-density gas-phase limit recovers the
// Spitzer/Oppenheimer-Dalgarno x_e ~ sqrt(zeta/(a_C n_H x_Ctot)).
//========================================================================================

// KOKKOS_INLINE_FUNCTION on device; plain inline for the standalone unit test.
#ifdef KOKKOS_INLINE_FUNCTION
#define CHEMG_FN KOKKOS_INLINE_FUNCTION
#else
#define CHEMG_FN inline
#endif

#include <cmath>

namespace Chemistry {

constexpr int NSPEC_GOW = 5;
constexpr int gH2 = 0;
constexpr int gHp = 1;
constexpr int gCp = 2;
constexpr int gCO = 3;
constexpr int ge  = 4;

struct Gow17ReducedNetwork {
  // --- cgs rate coefficients (representative; Gong+2017 / UMIST RATE12) ---
  double kgr = 3.0e-17;       // H2 grain formation        [cm^3 s^-1]
  double zeta = 1.0e-17;      // primary CR ionization rate [s^-1] (dense-core standard,
                              // Umebayashi-Nakano; matches ionization.hpp so the chem x_e
                              // and the equilibrium-eta cap use a consistent CR rate)
  double f_cr_H2 = 2.0;       // H2 CR destruction = f_cr_H2 * zeta
  double f_cr_H = 1.0;        // H  CR ionization  = f_cr_H  * zeta
  double f_cr_C = 3.0;        // C  ionization     = f_cr_C  * zeta (CR + CR-induced UV)
  double f_cr_CO = 10.0;      // CO CR destruction = f_cr_CO * zeta
  double a_rr0 = 2.7e-13;     // H+ radiative recomb a_rr = a_rr0 (T/1e4)^arr_exp [cm^3/s]
  double arr_exp = -0.70;
  double k_gr_Hp = 3.0e-17;   // H+ recomb on grains (per n_H)              [cm^3 s^-1]
  double a_C0 = 4.7e-12;      // C+ radiative recomb a_C = a_C0 (T/300)^aC_exp [cm^3/s]
  double aC_exp = -0.60;
  double k_CO = 5.0e-16;      // effective C+ -> CO formation (per n_H)     [cm^3 s^-1]

  // --- composition / units ---
  double x_Ctot = 1.6e-4;     // total gas-phase carbon abundance / n_H
  double rho_unit = 5.467e-19;
  double t_unit = 1.476998822551e12;
  double mu_n = 2.33;
  double m_H = 1.6726219e-24;
  double x_floor = 1.0e-20;
  double xe_floor = 1.0e-15;  // residual ionization (metals/CRs) floor on x_e
  int    nsub_max = 400;
  double cfl = 0.1;

  CHEMG_FN double nH(double rho_code) const {
    return rho_code * rho_unit / (mu_n * m_H);
  }

  // RHS in CODE time units for the 4 evolved species (e- is slaved). T in Kelvin.
  CHEMG_FN void rhs(const double *y, double n_H, double T, double *ydot) const {
    const double xH2 = y[gH2], xHp = y[gHp], xCp = y[gCp], xCO = y[gCO];
    double xH = 1.0 - 2.0 * xH2 - xHp;          if (xH  < 0.0) xH  = 0.0;
    double xC0 = x_Ctot - xCp - xCO;            if (xC0 < 0.0) xC0 = 0.0;
    double xe = xHp + xCp;                      if (xe  < xe_floor) xe = xe_floor;

    const double a_rr = a_rr0 * std::pow(T / 1.0e4, arr_exp);
    const double a_C  = a_C0  * std::pow(T / 300.0, aC_exp);

    const double R1 = kgr * n_H * xH;                          // H2 formation
    const double R2 = f_cr_H2 * zeta * xH2;                    // H2 CR dissoc
    const double R3 = f_cr_H  * zeta * xH;                     // H ionization
    const double R4 = (a_rr * n_H * xe + k_gr_Hp * n_H) * xHp; // H+ recomb (gas+grain)
    const double R5 = f_cr_C  * zeta * xC0;                    // C ionization
    const double R6 = a_C * n_H * xe * xCp;                    // C+ recomb
    const double R7 = k_CO * n_H * xCp;                        // CO formation
    const double R8 = f_cr_CO * zeta * xCO;                    // CO CR destruction

    ydot[gH2] = (R1 - R2) * t_unit;
    ydot[gHp] = (R3 - R4) * t_unit;
    ydot[gCp] = (R5 - R6 - R7 + R8) * t_unit;
    ydot[gCO] = (R7 - R8) * t_unit;
    ydot[ge]  = 0.0; // slaved: set in integrate_cell
  }

  // Production/loss split of the RHS in CODE time units: for each evolved species
  //   dy_s/dt = P[s] - L[s] * y[s],
  // with P[s] >= 0 (production) and L[s] >= 0 (linear loss-rate coefficient). Every
  // reaction that destroys species s is proportional to y[s], so the loss is linear in
  // y[s] at fixed coefficients. Within a sub-step P[s] and L[s] are FROZEN (L carries the
  // coupled x_e = x_Hp + x_Cp, so recombination is really quadratic in the destroyed
  // species; P[H2] hides a linear x_H2 loss via x_H) -- this is the standard Mott-Smith /
  // "asymptotic" linearization, first-order accurate, NOT an exact linear split. It enables
  // the unconditionally stable semi-implicit update in integrate_cell(), so the stiff (fast)
  // loss timescale no longer caps the sub-step. T in Kelvin.
  CHEMG_FN void rhs_PL(const double *y, double n_H, double T,
                       double *P, double *L) const {
    const double xH2 = y[gH2], xHp = y[gHp], xCp = y[gCp], xCO = y[gCO];
    double xH = 1.0 - 2.0 * xH2 - xHp;          if (xH  < 0.0) xH  = 0.0;
    double xC0 = x_Ctot - xCp - xCO;            if (xC0 < 0.0) xC0 = 0.0;
    double xe = xHp + xCp;                      if (xe  < xe_floor) xe = xe_floor;

    const double a_rr = a_rr0 * std::pow(T / 1.0e4, arr_exp);
    const double a_C  = a_C0  * std::pow(T / 300.0, aC_exp);

    // H2 : grain formation (P) ; CR dissociation (L ~ x_H2)
    P[gH2] = (kgr * n_H * xH) * t_unit;
    L[gH2] = (f_cr_H2 * zeta) * t_unit;
    // H+ : CR ionization of H (P) ; recombination gas+grain (L ~ x_Hp)
    P[gHp] = (f_cr_H * zeta * xH) * t_unit;
    L[gHp] = (a_rr * n_H * xe + k_gr_Hp * n_H) * t_unit;
    // C+ : C ionization + CO CR-destruction feed (P) ; recomb + CO formation (L ~ x_Cp)
    P[gCp] = (f_cr_C * zeta * xC0 + f_cr_CO * zeta * xCO) * t_unit;
    L[gCp] = (a_C * n_H * xe + k_CO * n_H) * t_unit;
    // CO : formation from C+ (P) ; CR destruction (L ~ x_CO)
    P[gCO] = (k_CO * n_H * xCp) * t_unit;
    L[gCO] = (f_cr_CO * zeta) * t_unit;
  }

  // Operator-split per-cell integration over code-time dt_code at fixed rho & T.
  // y[] abundances modified in place; T_code is the code-unit temperature, T_unit K/code.
  // Returns 1 if the integration was TRUNCATED (nsub_max sub-steps exhausted before
  // reaching dt_code, leaving the abundances under-integrated), else 0.
  CHEMG_FN int integrate_cell(double *y, double rho_code, double T_K, double dt_code) const {
    const double n_H = nH(rho_code);
    double t = 0.0;
    int nsub = 0;
    // Semi-implicit (production/loss-split) sub-cycler. Each species is advanced by
    //     y_new = (y + dt*P) / (1 + dt*L),
    // the linearized-implicit (asymptotic) update, which is unconditionally stable and
    // positivity-preserving for P,L >= 0. The stiff loss timescale therefore no longer
    // caps the step for STABILITY; the sub-step is limited only by an ACCURACY criterion
    // on the net change, and floored at dt_code/nsub_max so the full dt_code is always
    // covered in <= nsub_max steps (no truncation -- the old explicit Euler saturated
    // nsub_max because it had to resolve the fast loss timescale for stability).
    const double dt_floor = dt_code / (double)nsub_max;
    while (t < dt_code && nsub < nsub_max) {
      double P[NSPEC_GOW], L[NSPEC_GOW];
      rhs_PL(y, n_H, T_K, P, L);
      double dt_chem = dt_code - t;
      // accuracy limit from the NET tendency (r = P - L*y); stable at any magnitude.
      // Species already pinned at the abundance floor carry no accuracy information
      // (they stay at the floor), so they must NOT constrain the step -- otherwise their
      // O(1) relative residual would peg every sub-step at dt_floor forever.
      for (int s = 0; s < gCO + 1; ++s) {
        if (y[s] <= x_floor * 1.0000001) continue;
        double r  = P[s] - L[s] * y[s];
        double ad = (r < 0 ? -r : r);
        if (ad > 0.0) {
          double dtc = cfl * y[s] / ad;
          if (dtc < dt_chem) dt_chem = dtc;
        }
      }
      // Floor the step so stiffness can't stall progress: semi-implicit stability makes
      // these larger steps safe (they relax toward the local P/L equilibrium).
      if (dt_chem < dt_floor) dt_chem = dt_floor;
      if (t + dt_chem > dt_code) dt_chem = dt_code - t;
      for (int s = 0; s < gCO + 1; ++s) {
        y[s] = (y[s] + dt_chem * P[s]) / (1.0 + dt_chem * L[s]);
        if (y[s] < x_floor) y[s] = x_floor;
      }
      // Clamp to conservation limits. The per-species caps alone would still let the
      // SUMS 2 x_H2 + x_Hp and x_Cp + x_CO exceed their nucleus budgets (each species
      // at its individual cap), so joint renormalization enforces the budgets exactly.
      if (y[gH2] > 0.5) y[gH2] = 0.5;
      if (y[gHp] > 1.0) y[gHp] = 1.0;
      {
        const double hsum = 2.0 * y[gH2] + y[gHp];
        if (hsum > 1.0) {
          const double s = 1.0 / hsum;
          y[gH2] *= s;
          y[gHp] *= s;
        }
      }
      if (y[gCp] > x_Ctot) y[gCp] = x_Ctot;
      if (y[gCO] > x_Ctot) y[gCO] = x_Ctot;
      {
        const double csum = y[gCp] + y[gCO];
        if (csum > x_Ctot) {
          const double s = x_Ctot / csum;
          y[gCp] *= s;
          y[gCO] *= s;
        }
      }
      t += dt_chem;
      ++nsub;
    }
    // electron abundance slaved to charge neutrality (+ residual floor)
    double xe = y[gHp] + y[gCp];
    if (xe < xe_floor) xe = xe_floor;
    y[ge] = xe;
    // Truncated only if the sub-cycler genuinely failed to cover dt_code; a small
    // relative tolerance absorbs the floating-point undershoot from summing nsub_max
    // floored steps (which still fully integrate the abundances).
    return (t < dt_code * (1.0 - 1.0e-9)) ? 1 : 0;
  }
};

} // namespace Chemistry

#endif // CHEMISTRY_NETWORK_GOW17_REDUCED_HPP_
