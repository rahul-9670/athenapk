//========================================================================================
// AthenaPK - reduced ionization model for non-ideal MHD diffusivities (NICIL-class).
//
// Computes the gas/grain ionization state self-consistently (instead of a constant Q_A)
// and from it the three single-fluid magnetic diffusivities that set the internal B-field
// structure of the collapsing core:
//     ambipolar  eta_A   (ambipolar diffusion)
//     Hall       eta_H   (Hall drift; SIGNED -- sets the asymmetric field transport)
//     Ohmic      eta_O   (resistive dissipation)
//
// Charge carriers: electrons (e), one representative ion (i), and N_BIN MRN grain bins
// (dn/da ~ a^-3.5).  Each grain bin carries a mean charge <Z_k> set by electron/ion
// OML capture balance.  Ohmic and Hall come from inverting the Wardle conductivity tensor
// (sigma_O, sigma_H, sigma_P) built from e + i + all charged grain bins
// (Wardle 2007; Pandey & Wardle 2008; Marchand+2016; Wurster 2016 NICIL).
//
// AMBIPOLAR CLOSURE (ad_closure):
//   single_fluid (DEFAULT) : eta_A = B^2 / (4 pi gamma_AD rho_i rho)  -- the validated,
//       Athena++-matched single-fluid form (maps 1:1 to Athena++ eta_ad); grains enter
//       only through the x_e-dependent rho_i.  KEEP this for matched code comparisons.
//   tensor : eta_A = c^2/(4pi) sigma_P/sigma_perp^2 - eta_O  (grain-modified, ~3-5x lower;
//       floored at 0 where Hall/Ohmic dominate).  AthenaPK-only; NOT matched to Athena++.
//
// Ohmic and Hall are ALWAYS from the conductivity tensor (the single-fluid electron-only
// Hall over-predicts catastrophically at the ionization minimum -- the bug this fixes).
//
// Ionization sources/sinks ("CR + thermal + grains"):
//   (1) cosmic-ray ionization at rate zeta, balanced by gas-phase dissociative
//       recombination alpha(T) AND capture onto charged dust grains;
//   (2) thermal (Saha) ionization of potassium at T >~ 10^3 K (uses RT/EOS temperature);
//       grains sublimate above T_subl and drop out of the charge balance.
//
// Physics is evaluated in CGS internally (converting from the shared FHC code units via the
// unit scales); diffusivities are returned in code units (eta_code = eta_cgs t_unit/l_unit^2).
// All device-callable, no allocations; the MRN bin radii/mass-weights are filled host-side.
//========================================================================================
#ifndef HYDRO_DIFFUSION_IONIZATION_HPP_
#define HYDRO_DIFFUSION_IONIZATION_HPP_

#include <cmath>

#include <Kokkos_Core.hpp>

#include <basic_types.hpp>

namespace Ionization {

using parthenon::Real;

// Number of MRN grain bins (compile-time so the model is a POD captured by value).
constexpr int N_BIN = 5;

// Ambipolar-diffusion closure selector.
enum class ADClosure { single_fluid, tensor };

// CGS physical constants.
namespace cgs {
constexpr Real e_chg = 4.80320425e-10;  // elementary charge [esu]
constexpr Real m_e = 9.10938370e-28;    // electron mass [g]
constexpr Real m_H = 1.67262192e-24;    // hydrogen mass [g]
constexpr Real c_light = 2.99792458e10; // speed of light [cm/s]
constexpr Real k_B = 1.38064900e-16;    // Boltzmann [erg/K]
constexpr Real h_pl = 6.62607015e-27;   // Planck [erg s]
constexpr Real eV = 1.602176634e-12;    // [erg]
constexpr Real pi = 3.14159265358979324;
} // namespace cgs

//----------------------------------------------------------------------------------------
//! Reduced-ionization parameters. POD struct captured by value into device kernels.
struct IonizationModel {
  // --- unit conversions (code <-> cgs) ---
  // FLAGSHIP PHASE 1: these are set by BuildIonizationModel() from the single authoritative
  // PhysicalUnits (src/units/physical_units.hpp), NOT from these literals -- the literals are
  // only fallback defaults for standalone/test instantiation. T_unit uses mu_thermal = 2.29
  // (the c_s = 1.9e4 cm/s @ 10 K calibration), which is DISTINCT from the neutral mu_n = 2.33
  // below (n_n = rho/(mu_n m_H)); the two live in different objects so they can never be
  // conflated (that conflation was audit finding #4).
  Real rho_unit = 5.467e-19; // g/cm^3 per code density
  Real T_unit = 10.015;      // K per code temperature  (p = rho*T => T_code=1 <-> ~10 K)
  Real B_unit = 4.98e-5;     // G per code B (HL: B_unit = sqrt(4 pi rho_unit) v_unit)
  Real eta_unit = 1.874e-21; // eta_code = eta_cgs * (t_unit / l_unit^2)

  // --- gas / ionization parameters ---
  // zeta and mu_n are set by BuildIonizationModel() from the shared IonizationEnvironment
  // (src/units/ionization_environment.hpp), so Ohm/Hall/AD and chemistry share one CR rate.
  Real zeta = 1.0e-17;    // cosmic-ray ionization rate [s^-1]  (fallback default)
  Real mu_n = 2.33;       // mean molecular weight of neutrals
  Real m_ion = 24.3;      // representative ion mass [m_H] (metal/molecular, ~Mg+/HCO+)
  Real alpha0 = 2.4e-7;   // dissociative recomb alpha(T)=alpha0 (T/300)^alpha_Texp [cm^3/s]
  Real alpha_Texp = -0.69;
  Real sigma_en = 1.0e-15; // e-neutral momentum-transfer cross section [cm^2]
  Real sigv_in = 1.9e-9;   // ion-neutral (Langevin) rate coefficient [cm^3/s]

  // --- MRN dust grains (dn/da ~ a^-3.5, N_BIN bins; radii/weights filled host-side) ---
  Real f_dg = 0.01;       // dust-to-gas mass ratio
  Real rho_grain = 3.0;   // grain material density [g/cm^3]
  Real stick_e = 1.0;     // electron sticking coefficient
  Real stick_i = 1.0;     // ion sticking coefficient
  Real T_subl = 1.5e3;    // grain sublimation temperature [K]
  Real a_k[N_BIN] = {0.0}; // bin geometric-mean radius [cm]  (host-filled)
  Real mw_k[N_BIN] = {0.0}; // bin dust-mass fraction          (host-filled)

  // --- thermal (Saha) ionization of potassium ---
  Real x_K = 1.0e-7;     // potassium abundance relative to H nuclei
  Real chi_K = 4.34;     // first ionization potential [eV]
  Real x_H = 0.716;      // hydrogen MASS fraction (thermal H Saha: n_H = n_n mu_n x_H)

  // floors / options
  Real xe_floor = 1.0e-20;                       // minimum ionization fraction
  ADClosure ad_closure = ADClosure::single_fluid; // ambipolar closure (see header)
  Real gamma_AD = 3.5e13; // ion-neutral drag coeff for single_fluid AD [cm^3 g^-1 s^-1]
  // Charge-state solvers (SolveCharges + SahaThermal). false = the robust solves (default);
  // true = the pre-2026-07-25 pair, BOTH of which are defective: the grain fixed point does
  // not converge once grains dominate (neutrality violated by up to ~1e7), and the thermal
  // Saha bisection has an absolute resolution floor that injects spurious electrons in cold
  // gas. Kept only to bit-reproduce pre-fix results. See the Phase-3 conductivity gate.
  bool legacy_charge_solver = false;
};

//----------------------------------------------------------------------------------------
//! Fill the MRN bin radii (geometric-mean) and dust-mass fractions on the HOST.
//! dn/da ~ a^-p over [a_min, a_max]; mass per bin ~ integral a^(3-p) da.
inline void SetupGrainBins(IonizationModel &m, const Real a_min, const Real a_max,
                           const Real p = 3.5) {
  Real edges[N_BIN + 1];
  const Real llo = std::log10(a_min), lhi = std::log10(a_max);
  for (int e = 0; e <= N_BIN; ++e) {
    edges[e] = std::pow(10.0, llo + (lhi - llo) * static_cast<Real>(e) / N_BIN);
  }
  auto integ = [](Real lo, Real hi, Real q) {
    return (std::abs(q + 1.0) < 1e-9) ? std::log(hi / lo)
                                      : (std::pow(hi, q + 1.0) - std::pow(lo, q + 1.0)) /
                                            (q + 1.0);
  };
  Real msum = 0.0;
  for (int k = 0; k < N_BIN; ++k) {
    m.a_k[k] = std::sqrt(edges[k] * edges[k + 1]);
    m.mw_k[k] = integ(edges[k], edges[k + 1], 3.0 - p);
    msum += m.mw_k[k];
  }
  for (int k = 0; k < N_BIN; ++k)
    m.mw_k[k] /= msum;
}

//----------------------------------------------------------------------------------------
//! Mean thermal speed sqrt(8 k T / pi m).
KOKKOS_INLINE_FUNCTION
Real VBar(const Real T, const Real mass) {
  return std::sqrt(8.0 * cgs::k_B * T / (cgs::pi * mass));
}

//----------------------------------------------------------------------------------------
//! Thermal (Saha) free-electron density [cm^-3] from potassium AND hydrogen ionization.
//! K (chi=4.34 eV) provides the ~10^3 K electrons that re-couple the field once grains
//! sublimate; H (chi=13.6 eV) dominates from ~5000 K into the second core (K alone
//! saturates at x_e ~ x_K ~ 1e-7 and badly under-ionizes there). Both donors share n_e,
//! so solve neutrality n_e = sum_i n_i,tot f_i/(f_i + n_e) by bisection. The H term uses
//! the same Saha as the tabulated EOS, keeping the conductivity ionization consistent with
//! the thermodynamic ionization at high T.
KOKKOS_INLINE_FUNCTION
Real SahaThermal(const IonizationModel &m, const Real n_n, const Real T) {
  const Real le3 = std::pow(2.0 * cgs::pi * cgs::m_e * cgs::k_B * T / (cgs::h_pl * cgs::h_pl),
                            1.5);
  const Real fK = le3 * std::exp(-m.chi_K * cgs::eV / (cgs::k_B * T));  // K: g-ratio ~1
  const Real fH = le3 * std::exp(-13.598 * cgs::eV / (cgs::k_B * T));   // H: gII/gI net 1
  const Real n_K = m.x_K * n_n;
  const Real n_H = n_n * m.mu_n * m.x_H; // H-nucleus density = rho X_H/m_H = n_n mu_n X_H

  if (m.legacy_charge_solver) {
    // LEGACY (pre-2026-07-25): absolute bisection on [0, n_K+n_H] with a FIXED 64 steps.
    // Its resolution floor is (n_K+n_H)/2^64, so in cold gas -- where the true thermal
    // ionization underflows -- it returns ~(n_K+n_H)*2^-65 instead of ~0: a spurious
    // electron floor of x_e ~ 9e-20 (about 9x the intended xe_floor) that grows with
    // density and can dominate the real n_e in the dense core. Kept for bit-reproduction.
    Real lo = 0.0, hi = n_K + n_H;
    for (int it = 0; it < 64; ++it) {
      const Real ne = 0.5 * (lo + hi);
      const Real res = n_K * fK / (fK + ne) + n_H * fH / (fH + ne) - ne;
      if (res > 0.0) lo = ne; else hi = ne;
    }
    return 0.5 * (lo + hi);
  }

  // Relative-precision solve of  n_e = n_K fK/(fK+n_e) + n_H fH/(fH+n_e).
  // Newton from the weak-ionization limit n_e ~ sqrt(n_K fK + n_H fH) (exact when n_e >> f);
  // capped by full ionization n_K + n_H. Accurate down to underflow, unlike the fixed-step
  // absolute bisection above.
  const Real sum_nf = n_K * fK + n_H * fH;
  if (!(sum_nf > 0.0)) return 0.0; // both Boltzmann factors underflowed: no thermal electrons
  const Real ne_max = n_K + n_H;
  Real ne = std::sqrt(sum_nf);
  if (ne > ne_max) ne = ne_max;
  for (int it = 0; it < 100; ++it) {
    const Real gK = n_K * fK / (fK + ne);
    const Real gH = n_H * fH / (fH + ne);
    const Real g = ne - gK - gH;
    const Real dg = 1.0 + gK / (fK + ne) + gH / (fH + ne);
    Real dne = -g / dg;
    if (ne + dne <= 0.0) dne = -0.5 * ne; // keep the iterate positive
    ne += dne;
    if (std::abs(dne) <= 1e-15 * ne) break;
  }
  if (ne > ne_max) ne = ne_max;
  return ne;
}

//----------------------------------------------------------------------------------------
//! Self-consistent charge state: electrons n_e, ions n_i [cm^-3], and per-bin grain mean
//! charge Z_k (signed) and number density ng_k [cm^-3], given neutral density n_n and T.
//! CR ionization + grain capture solved by fixed-point; thermal Saha added on top.
KOKKOS_INLINE_FUNCTION
void SolveCharges(const IonizationModel &m, const Real n_n, const Real T, Real &n_e,
                  Real &n_i, Real Zk[N_BIN], Real ng[N_BIN]) {
  // grain number densities per bin (zero if sublimated)
  const Real rho_d = m.f_dg * (n_n * m.mu_n * cgs::m_H);
  Real tau[N_BIN], ke0[N_BIN], ki0[N_BIN];
  const Real ve = VBar(T, cgs::m_e);
  const Real vi = VBar(T, m.m_ion * cgs::m_H);
  for (int k = 0; k < N_BIN; ++k) {
    const Real a = m.a_k[k];
    const Real m_g = (4.0 / 3.0) * cgs::pi * a * a * a * m.rho_grain;
    ng[k] = (T >= m.T_subl) ? 0.0 : m.mw_k[k] * rho_d / m_g;
    tau[k] = a * cgs::k_B * T / (cgs::e_chg * cgs::e_chg); // reduced temperature
    ke0[k] = cgs::pi * a * a * ve * m.stick_e;             // neutral-grain e capture
    ki0[k] = cgs::pi * a * a * vi * m.stick_i;             // neutral-grain i capture
    Zk[k] = -0.5;
  }

  const Real alpha = m.alpha0 * std::pow(T / 300.0, m.alpha_Texp);
  const Real P = m.zeta * n_n; // CR ion-pair production rate

  // ---------------------------------------------------------------------------------
  // ROBUST grain-charge solve (default; Phase-3 fix 2026-07-25).
  //
  // The legacy path below is a relaxed fixed point on n_i with an inner Newton on psi; it
  // DIVERGES once grains dominate the charge budget (rho >~ 1e-12 g/cm^3 on the collapse
  // track), returning states that violate its own neutrality constraint by up to ~1e7 in
  // relative terms -- and even n_e > n_i, impossible with negatively charged grains.
  //
  // Reduction used here: the per-bin capture balance n_e v_e e^psi = n_i v_i (1-psi) has NO
  // bin dependence (the pi a^2 cross-sections cancel between electrons and ions), so psi is a
  // single global unknown with Z_k = psi tau_k. That makes r = n_e/n_i an EXPLICIT function
  //   r(psi) = (v_i/v_e) (1 - psi) e^{-psi},
  // and the remaining constraints close on psi alone:
  //   neutrality : n_i (1 - r) = Q = -psi S,      S  = sum_k tau_k ng_k
  //   ionization : alpha r n_i^2 + (1-psi) K0 n_i = P,   K0 = sum_k ng_k ki0_k
  // R(psi) is +inf at psi_lo (where r -> 1, n_i -> inf) and -P at psi -> 0^-, so a bracketed
  // bisection converges unconditionally -- no relaxation, no iteration-count failure mode.
  // ---------------------------------------------------------------------------------
  Real Ssum = 0.0, K0 = 0.0;
  for (int k = 0; k < N_BIN; ++k) {
    Ssum += tau[k] * ng[k];
    K0 += ng[k] * ki0[k];
  }

  if (!m.legacy_charge_solver) {
    if (Ssum <= 0.0) { // no grains (sublimated or f_dg = 0): pure gas-phase balance
      n_i = std::sqrt(P / alpha);
      n_e = n_i;
      for (int k = 0; k < N_BIN; ++k)
        Zk[k] = 0.0;
    } else {
      const Real vr = vi / ve; // = r at psi = 0
      // psi_lo: the psi where r(psi) = 1 (n_i -> inf). r is strictly decreasing in psi.
      Real lo = -200.0, hi = 0.0;
      for (int it = 0; it < 60; ++it) {
        const Real mid = 0.5 * (lo + hi);
        const Real ex = std::exp((-mid > 200.0) ? 200.0 : -mid);
        if (vr * (1.0 - mid) * ex > 1.0) lo = mid;
        else hi = mid;
      }
      // bracket [a,b] with R(a) > 0 > R(b); b = 0^- gives R = -P.
      Real a = hi, b = 0.0;
      for (int it = 0; it < 80; ++it) {
        const Real psi = 0.5 * (a + b);
        const Real ex = std::exp((-psi > 200.0) ? 200.0 : -psi);
        const Real r = vr * (1.0 - psi) * ex;
        const Real ni = -psi * Ssum / (1.0 - r);
        const Real R = alpha * r * ni * ni + (1.0 - psi) * K0 * ni - P;
        if (R > 0.0) a = psi;
        else b = psi;
      }
      const Real psi = 0.5 * (a + b);
      const Real ex = std::exp((-psi > 200.0) ? 200.0 : -psi);
      const Real r = vr * (1.0 - psi) * ex;
      n_i = -psi * Ssum / (1.0 - r);
      n_e = r * n_i;
      if (n_e < 1e-30) n_e = 1e-30;
      for (int k = 0; k < N_BIN; ++k)
        Zk[k] = (ng[k] > 0.0) ? psi * tau[k] : 0.0;
    }
    // thermal channel + floor (shared tail with the legacy path)
    const Real n_th_r = SahaThermal(m, n_n, T);
    n_e += n_th_r;
    n_i += n_th_r;
    const Real n_e_min_r = m.xe_floor * n_n;
    if (n_e < n_e_min_r) n_e = n_e_min_r;
    return;
  }

  // ---------------- LEGACY relaxed fixed point (non-convergent; see above) -------------
  // initial guess: pure gas-phase recombination
  n_i = std::sqrt(P / alpha);
  n_e = n_i;

  const Real relax = 0.5;
  for (int it = 0; it < 300; ++it) {
    // per-bin mean charge from capture balance n_e ve exp(psi) = n_i vi (1 - psi)
    Real Ki = 0.0; // total ion->grain capture coefficient [s^-1 per n_i]
    Real Q = 0.0;  // grain negative-charge density [cm^-3]
    Real maxdZ = 0.0;
    for (int k = 0; k < N_BIN; ++k) {
      if (ng[k] <= 0.0) {
        Zk[k] = 0.0;
        continue;
      }
      Real psi = Zk[k] / tau[k];
      for (int j = 0; j < 60; ++j) {
        const Real ex = std::exp((psi < -200.0) ? -200.0 : (psi > 50.0 ? 50.0 : psi));
        const Real f = n_e * ve * ex - n_i * vi * (1.0 - psi);
        const Real df = n_e * ve * ex + n_i * vi;
        Real dp = -f / df;
        if (dp > 5.0) dp = 5.0;
        if (dp < -5.0) dp = -5.0;
        psi += dp;
        if (std::abs(dp) < 1e-10) break;
      }
      const Real Zk_new = psi * tau[k];
      maxdZ = std::max(maxdZ, std::abs(Zk_new - Zk[k]));
      Zk[k] = Zk_new;
      Ki += ng[k] * ki0[k] * (1.0 - psi);
      Q += -Zk_new * ng[k];
    }

    // ionization balance with neutrality n_e = n_i - Q:
    //   P = alpha (n_i - Q) n_i + Ki n_i  =>  alpha n_i^2 + (Ki - alpha Q) n_i - P = 0
    const Real b = Ki - alpha * Q;
    const Real disc = b * b + 4.0 * alpha * P;
    const Real n_i_new = (-b + std::sqrt(disc)) / (2.0 * alpha);
    Real n_e_new = n_i_new - Q;
    if (n_e_new < 1e-30) n_e_new = 1e-30;

    const Real dni = std::abs(n_i_new - n_i);
    n_i = (1.0 - relax) * n_i + relax * n_i_new;
    n_e = (1.0 - relax) * n_e + relax * n_e_new;
    if (dni < 1e-10 * n_i && maxdZ < 1e-8) break;
  }

  // thermal channel (electrons & ions; grains gone at these temperatures)
  const Real n_th = SahaThermal(m, n_n, T);
  n_e += n_th;
  n_i += n_th;
  const Real n_e_min = m.xe_floor * n_n;
  if (n_e < n_e_min) n_e = n_e_min;
}

//----------------------------------------------------------------------------------------
//! All three magnetic diffusivities in CODE units, given code-unit (rho, T, |B|).
//! eta_O (Ohmic) and eta_H (Hall, signed) from the Wardle conductivity tensor over
//! e + i + grain bins; eta_A per ad_closure (single_fluid default, or tensor).
KOKKOS_INLINE_FUNCTION
void Diffusivities(const IonizationModel &m, const Real rho_code, const Real T_code,
                   const Real B_code, Real &eta_O, Real &eta_H, Real &eta_A) {
  const Real rho = rho_code * m.rho_unit; // g/cm^3
  const Real T = T_code * m.T_unit;       // K
  // Floor |B|: the tensor assembly divides by B (ecB = e*c/B), so B = 0 exactly (e.g. a
  // field null in a turbulent run) would give 0*inf = NaN. 1e-20 G is far below any
  // physical field; all eta's go smoothly to their B->0 limits.
  const Real B = std::max(std::abs(B_code) * m.B_unit, 1.0e-20); // G
  const Real m_n = m.mu_n * cgs::m_H;
  const Real n_n = rho / m_n; // cm^-3

  Real n_e, n_i, Zk[N_BIN], ng[N_BIN];
  SolveCharges(m, n_n, T, n_e, n_i, Zk, ng);

  // collision frequencies with neutrals  nu_jn = (m_n/(m_j+m_n)) n_n <sigv>_jn
  const Real ve = VBar(T, cgs::m_e);
  const Real m_ion = m.m_ion * cgs::m_H;
  const Real nu_e = (m_n / (cgs::m_e + m_n)) * n_n * (m.sigma_en * ve);
  const Real nu_i = (m_n / (m_ion + m_n)) * n_n * m.sigv_in;
  const Real vn = VBar(T, m_n);

  // Conductivity tensor (Pandey & Wardle 2008): sigma_{O,H,P} = (e c / B) sum_j n_j Z_j f_j
  //   with signed Hall parameter beta_j = Z_j e B / (m_j c) / nu_j.
  const Real ecB = cgs::e_chg * cgs::c_light / B;
  auto beta = [&](Real Zsigned, Real mass, Real nu) {
    return (Zsigned * cgs::e_chg * B) / (mass * cgs::c_light) / nu;
  };
  // electrons (Z=-1) and ions (Z=+1)
  Real sO = 0.0, sH = 0.0, sP = 0.0;
  {
    const Real be = beta(-1.0, cgs::m_e, nu_e);
    sO += n_e * (-1.0) * be;
    sH += n_e * (-1.0) / (1.0 + be * be);
    sP += n_e * (-1.0) * be / (1.0 + be * be);
    const Real bi = beta(1.0, m_ion, nu_i);
    sO += n_i * (1.0) * bi;
    sH += n_i * (1.0) / (1.0 + bi * bi);
    sP += n_i * (1.0) * bi / (1.0 + bi * bi);
  }
  // charged grain bins (mean charge Z_k, mass m_gk)
  for (int k = 0; k < N_BIN; ++k) {
    if (ng[k] <= 0.0 || std::abs(Zk[k]) < 1e-30) continue;
    const Real a = m.a_k[k];
    const Real m_gk = (4.0 / 3.0) * cgs::pi * a * a * a * m.rho_grain;
    const Real nu_g = (m_n / (m_gk + m_n)) * n_n * (cgs::pi * a * a * vn);
    const Real bg = beta(Zk[k], m_gk, nu_g);
    sO += ng[k] * Zk[k] * bg;
    sH += ng[k] * Zk[k] / (1.0 + bg * bg);
    sP += ng[k] * Zk[k] * bg / (1.0 + bg * bg);
  }
  sO *= ecB;
  sH *= ecB;
  sP *= ecB;

  const Real pref = cgs::c_light * cgs::c_light / (4.0 * cgs::pi);
  const Real sperp2 = sH * sH + sP * sP;
  const Real eta_O_cgs = pref / sO;          // Ohmic = c^2/(4 pi sigma_O)
  const Real eta_H_cgs = pref * sH / sperp2; // Hall (signed)

  // Ambipolar closure.
  Real eta_A_cgs;
  if (m.ad_closure == ADClosure::tensor) {
    eta_A_cgs = pref * sP / sperp2 - eta_O_cgs; // grain-modified; floor at 0
    if (eta_A_cgs < 0.0) eta_A_cgs = 0.0;
  } else {
    // single-fluid (Athena++-matched): eta_A = B^2 / (4 pi gamma_AD rho_i rho)
    const Real rho_i = n_i * m_ion;
    eta_A_cgs = B * B / (4.0 * cgs::pi * m.gamma_AD * rho_i * rho);
  }

  // Convert to code units.
  eta_O = eta_O_cgs * m.eta_unit;
  eta_H = eta_H_cgs * m.eta_unit;
  eta_A = eta_A_cgs * m.eta_unit;
}

//----------------------------------------------------------------------------------------
//! Grain charge state with the FREE-ELECTRON density FIXED by chemistry (n_e = x_e n_n),
//! instead of solved from CR ionization equilibrium. The grains equilibrate to the ambient
//! (n_e, n_i) on times << t_dyn, and charge neutrality closes the ion density:
//!   n_i = n_e + Q,   Q = sum_k (-Z_k) ng_k  (grain net negative charge density).
//! This is the inner grain loop of SolveCharges with n_e held fixed.
KOKKOS_INLINE_FUNCTION
void SolveGrainsFixedNe(const IonizationModel &m, const Real n_n, const Real T,
                        const Real n_e_chem, Real &n_e, Real &n_i, Real Zk[N_BIN],
                        Real ng[N_BIN]) {
  const Real rho_d = m.f_dg * (n_n * m.mu_n * cgs::m_H);
  Real tau[N_BIN];
  const Real ve = VBar(T, cgs::m_e);
  const Real vi = VBar(T, m.m_ion * cgs::m_H);
  for (int k = 0; k < N_BIN; ++k) {
    const Real a = m.a_k[k];
    const Real m_g = (4.0 / 3.0) * cgs::pi * a * a * a * m.rho_grain;
    ng[k] = (T >= m.T_subl) ? 0.0 : m.mw_k[k] * rho_d / m_g;
    tau[k] = a * cgs::k_B * T / (cgs::e_chg * cgs::e_chg);
    Zk[k] = -0.5;
  }
  n_e = n_e_chem; // fixed by the evolved chemistry x_e
  n_i = n_e_chem; // initial guess; corrected each sweep by neutrality with grains
  const Real relax = 0.5;
  for (int it = 0; it < 300; ++it) {
    Real Q = 0.0, maxdZ = 0.0;
    for (int k = 0; k < N_BIN; ++k) {
      if (ng[k] <= 0.0) {
        Zk[k] = 0.0;
        continue;
      }
      Real psi = Zk[k] / tau[k];
      for (int j = 0; j < 60; ++j) {
        const Real ex = std::exp((psi < -200.0) ? -200.0 : (psi > 50.0 ? 50.0 : psi));
        const Real f = n_e * ve * ex - n_i * vi * (1.0 - psi);
        const Real df = n_e * ve * ex + n_i * vi;
        Real dp = -f / df;
        if (dp > 5.0) dp = 5.0;
        if (dp < -5.0) dp = -5.0;
        psi += dp;
        if (std::abs(dp) < 1e-10) break;
      }
      const Real Zk_new = psi * tau[k];
      maxdZ = std::max(maxdZ, std::abs(Zk_new - Zk[k]));
      Zk[k] = Zk_new;
      Q += -Zk_new * ng[k];
    }
    Real n_i_new = n_e + Q; // neutrality: n_i(+) = n_e(-) + grain net negative charge
    if (n_i_new < n_e) n_i_new = n_e;
    const Real dni = std::abs(n_i_new - n_i);
    n_i = (1.0 - relax) * n_i + relax * n_i_new;
    if (dni < 1e-10 * n_i && maxdZ < 1e-8) break;
  }
}

//----------------------------------------------------------------------------------------
//! All three magnetic diffusivities in CODE units from an EXTERNALLY-supplied electron
//! abundance x_e (evolved by the chemistry package), via the SAME Wardle conductivity
//! tensor (e + i + charged grains) as Diffusivities() -- only the electron density comes
//! from chemistry instead of the CR-equilibrium charge solve. eta_O, eta_H (signed) from
//! the tensor; eta_A single-fluid (Athena++-matched) from the ion density. Using the smooth,
//! monotonic chemistry x_e here also removes the high-density non-convergence noise of the
//! equilibrium solve at the ionization minimum.
KOKKOS_INLINE_FUNCTION
void DiffusivitiesFromXe(const IonizationModel &m, const Real rho_code, const Real T_code,
                         const Real B_code, const Real xe, Real &eta_O, Real &eta_H,
                         Real &eta_A) {
  const Real rho = rho_code * m.rho_unit;
  const Real T = T_code * m.T_unit;
  // Floor |B| against the 1/B in the tensor assembly (see Diffusivities).
  const Real B = std::max(std::abs(B_code) * m.B_unit, 1.0e-20);
  const Real m_n = m.mu_n * cgs::m_H;
  const Real n_n = rho / m_n;
  const Real xe_use = (xe > m.xe_floor) ? xe : m.xe_floor;

  Real n_e, n_i, Zk[N_BIN], ng[N_BIN];
  SolveGrainsFixedNe(m, n_n, T, xe_use * n_n, n_e, n_i, Zk, ng);

  const Real ve = VBar(T, cgs::m_e);
  const Real m_ion = m.m_ion * cgs::m_H;
  const Real nu_e = (m_n / (cgs::m_e + m_n)) * n_n * (m.sigma_en * ve);
  const Real nu_i = (m_n / (m_ion + m_n)) * n_n * m.sigv_in;
  const Real vn = VBar(T, m_n);

  const Real ecB = cgs::e_chg * cgs::c_light / B;
  auto beta = [&](Real Zsigned, Real mass, Real nu) {
    return (Zsigned * cgs::e_chg * B) / (mass * cgs::c_light) / nu;
  };
  Real sO = 0.0, sH = 0.0, sP = 0.0;
  {
    const Real be = beta(-1.0, cgs::m_e, nu_e);
    sO += n_e * (-1.0) * be;
    sH += n_e * (-1.0) / (1.0 + be * be);
    sP += n_e * (-1.0) * be / (1.0 + be * be);
    const Real bi = beta(1.0, m_ion, nu_i);
    sO += n_i * (1.0) * bi;
    sH += n_i * (1.0) / (1.0 + bi * bi);
    sP += n_i * (1.0) * bi / (1.0 + bi * bi);
  }
  for (int k = 0; k < N_BIN; ++k) {
    if (ng[k] <= 0.0 || std::abs(Zk[k]) < 1e-30) continue;
    const Real a = m.a_k[k];
    const Real m_gk = (4.0 / 3.0) * cgs::pi * a * a * a * m.rho_grain;
    const Real nu_g = (m_n / (m_gk + m_n)) * n_n * (cgs::pi * a * a * vn);
    const Real bg = beta(Zk[k], m_gk, nu_g);
    sO += ng[k] * Zk[k] * bg;
    sH += ng[k] * Zk[k] / (1.0 + bg * bg);
    sP += ng[k] * Zk[k] * bg / (1.0 + bg * bg);
  }
  sO *= ecB;
  sH *= ecB;
  sP *= ecB;

  const Real pref = cgs::c_light * cgs::c_light / (4.0 * cgs::pi);
  const Real sperp2 = sH * sH + sP * sP;
  const Real eta_O_cgs = pref / sO;
  const Real eta_H_cgs = pref * sH / sperp2;
  // Ambipolar: single-fluid (Athena++-matched) from the chemistry-set ion density.
  const Real rho_i = n_i * m_ion;
  const Real eta_A_cgs = B * B / (4.0 * cgs::pi * m.gamma_AD * rho_i * rho);

  eta_O = eta_O_cgs * m.eta_unit;
  eta_H = eta_H_cgs * m.eta_unit;
  eta_A = eta_A_cgs * m.eta_unit;
}

//----------------------------------------------------------------------------------------
//! Single-fluid ambipolar diffusivity (CODE units) from an EXTERNALLY-supplied electron
//! abundance x_e (e.g. evolved time-dependently by the chemistry package), bypassing the
//! equilibrium charge solve:  eta_A = B^2 / (4 pi gamma_AD rho_i rho),  rho_i = x_e n_n m_ion
//! (ions track electrons by charge balance). This is the Athena++-matched single-fluid form
//! with the only change being that x_e comes from chemistry instead of CR equilibrium.
KOKKOS_INLINE_FUNCTION
Real AmbipolarEtaFromXe(const IonizationModel &m, const Real rho_code, const Real B_code,
                        const Real xe) {
  const Real rho = rho_code * m.rho_unit;       // g/cm^3
  const Real B = B_code * m.B_unit;             // G
  const Real m_n = m.mu_n * cgs::m_H;
  const Real n_n = rho / m_n;                   // cm^-3
  Real xe_use = (xe > m.xe_floor) ? xe : m.xe_floor;
  const Real n_i = xe_use * n_n;                // ion density ~ electron density
  const Real rho_i = n_i * (m.m_ion * cgs::m_H);
  const Real eta_A_cgs = B * B / (4.0 * cgs::pi * m.gamma_AD * rho_i * rho);
  return eta_A_cgs * m.eta_unit;
}

} // namespace Ionization

#endif // HYDRO_DIFFUSION_IONIZATION_HPP_
