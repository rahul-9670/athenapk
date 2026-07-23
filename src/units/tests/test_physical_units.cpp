// Flagship Phase 1 acceptance test for the authoritative unit system.
//
// Host g++ standalone (same convention as src/chemistry/tests/test_thermo.cpp): this test
// INDEPENDENTLY re-derives the FHC unit scales from the base (rho, v, length, mu) and checks
// the physical identities the PhysicalUnits header must satisfy. Because physical_units.hpp
// pulls in Parthenon (parameter_input / package), it cannot be compiled host-only; the
// derivation formulas below are a deliberate independent mirror -- keep them in sync with
// src/units/physical_units.hpp (BuildPhysicalUnits) and src/units.hpp (CodeMagneticCgs).
// Phase 1 follow-up (pending the units-authority decision): collapse this mirror onto a
// dependency-free units_core.hpp shared with the header so there is a single derivation.
//
// Gates:
//   (1) cgs round-trip: physical quantity -> code -> cgs recovers the input (1000 samples).
//   (2) Heaviside-Lorentz <-> Gaussian bridge: P_mag = B_code^2/2 (HL, code) equals
//       B_cgs^2/(8*pi) (Gaussian, cgs) converted back to code. This is the #1 analysis
//       landmine (the sqrt(4*pi)); it must hold to round-off. It reduces to the exact
//       identity  energy_density_unit == magnetic_unit^2 / (4*pi).
//   (3) single temperature calibration (audit #4): temperature_unit depends ONLY on (mu, v),
//       independent of rho/length; mu=2.29 -> T0 = 10.015 K, the one value every package uses.
//   (4) RSLA / radiation sanity: c_code and arad_code finite, positive, correct magnitude.
//   (5) DRIFT SENTINEL (informational): PhysicalUnits <units>-block defaults
//       (5.467e-19, 1.9e4, 2.81e16) vs the EXACT BE normalization derived from the fiducial
//       IC (mass=6, temp=10, f=5). diffusivity_unit and opacity_unit differ by ~1.3e-3 --
//       the Phase-1 finding. A consistent config drives every column to 0.
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace {
constexpr double PI = 3.14159265358979323846;
// cgs constants -- MUST match src/units/physical_units.hpp namespace cgs.
constexpr double m_H = 1.6726219e-24;     // g
constexpr double k_B = 1.380649e-16;      // erg/K
constexpr double a_rad = 7.5657e-15;      // erg cm^-3 K^-4
constexpr double c_light = 2.99792458e10; // cm/s

// Independent mirror of BuildPhysicalUnits' derivation (see file header).
struct Scales {
  double rho, v, length, mu; // base
  double time, mass, energy_density, magnetic, temperature, diffusivity, opacity;
};
Scales Derive(double rho, double v, double length, double mu) {
  Scales s;
  s.rho = rho; s.v = v; s.length = length; s.mu = mu;
  s.time = length / v;
  s.mass = rho * length * length * length;
  s.energy_density = rho * v * v;
  // CodeMagneticCgs(mass, length, time) = sqrt(4pi)*sqrt(m)/sqrt(l)/t  (Heaviside-Lorentz).
  s.magnetic = std::sqrt(4.0 * PI) * std::sqrt(s.mass) / std::sqrt(length) / s.time;
  s.temperature = mu * m_H * v * v / k_B;
  s.diffusivity = s.time / (length * length);
  s.opacity = rho * length;
  return s;
}
double arad_code(const Scales &s) {
  const double T2 = s.temperature * s.temperature;
  return a_rad * T2 * T2 / s.energy_density;
}
double relerr(double a, double b) { return std::fabs(a - b) / std::fabs(b == 0.0 ? 1.0 : b); }

// Deterministic LCG so the 1000-sample sweep is reproducible.
struct Rng {
  unsigned long s = 88172645463325252ULL;
  double next() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return (s >> 11) * (1.0 / 9007199254740992.0); }
  double loguniform(double lo, double hi) { return std::pow(10.0, std::log10(lo) + next() * (std::log10(hi / lo))); }
};
} // namespace

int main() {
  int fail = 0;
  const double tol = 1.0e-12;

  // Fiducial FHC normalization (the <units>-block defaults BuildPhysicalUnits ships).
  const Scales U = Derive(5.467e-19, 1.9e4, 2.81e16, 2.29);

  // (1) cgs round-trip on 1000 samples across many decades.
  printf("--- GATE 1: cgs round-trip (1000 samples) ---\n");
  Rng rng;
  double worst_rt = 0.0;
  for (int i = 0; i < 1000; ++i) {
    const double rho_cgs = rng.loguniform(1e-20, 1e-8);   // g/cm^3
    const double B_cgs = rng.loguniform(1e-7, 1.0);       // G
    const double kap_cgs = rng.loguniform(1e-3, 10.0);    // cm^2/g
    const double eta_cgs = rng.loguniform(1e15, 1e22);    // cm^2/s
    const double rho_code = rho_cgs / U.rho; // divide by unit to go cgs->code
    const double B_code = B_cgs / U.magnetic;
    const double kap_code = kap_cgs * U.opacity;   // kappa_code = kappa_cgs * opacity_unit
    const double eta_code = eta_cgs * U.diffusivity;
    worst_rt = std::max(worst_rt, relerr(rho_code * U.rho, rho_cgs));
    worst_rt = std::max(worst_rt, relerr(B_code * U.magnetic, B_cgs));
    worst_rt = std::max(worst_rt, relerr(kap_code / U.opacity, kap_cgs));
    worst_rt = std::max(worst_rt, relerr(eta_code / U.diffusivity, eta_cgs));
  }
  const bool g1 = worst_rt <= tol;
  printf("  worst round-trip rel err = %.2e  %s\n", worst_rt, g1 ? "PASS" : "FAIL");
  fail += !g1;

  // (2) Heaviside-Lorentz <-> Gaussian magnetic-pressure bridge.
  printf("--- GATE 2: HL P_mag=B^2/2 == Gaussian B_cgs^2/8pi ---\n");
  double worst_hl = 0.0;
  for (int i = 0; i < 1000; ++i) {
    const double B_code = rng.loguniform(1e-4, 1e2);
    const double pmag_hl_code = 0.5 * B_code * B_code;                 // HL, code units
    const double B_cgs = B_code * U.magnetic;
    const double pmag_gauss_code = (B_cgs * B_cgs / (8.0 * PI)) / U.energy_density; // Gaussian->code
    worst_hl = std::max(worst_hl, relerr(pmag_gauss_code, pmag_hl_code));
  }
  // Underlying exact identity: energy_density_unit == magnetic_unit^2/(4pi).
  const double id_err = relerr(U.energy_density, U.magnetic * U.magnetic / (4.0 * PI));
  const bool g2 = (worst_hl <= tol) && (id_err <= tol);
  printf("  worst P_mag rel err = %.2e ; energy_density==B_unit^2/4pi rel err = %.2e  %s\n",
         worst_hl, id_err, g2 ? "PASS" : "FAIL");
  fail += !g2;

  // (3) single temperature calibration: T_unit(mu,v) independent of rho,length; T0=10.015 K.
  printf("--- GATE 3: single T0 (audit #4) ---\n");
  const Scales A = Derive(1e-19, 1.9e4, 1e16, 2.29); // "chemistry"-like base
  const Scales B = Derive(9e-18, 1.9e4, 5e17, 2.29); // "radiation"-like base
  const double t_consistency = relerr(A.temperature, B.temperature);
  const double T0 = U.temperature; // K per code T; c_s,iso=1 at the reference T
  const bool g3 = (t_consistency <= tol) && (relerr(T0, 10.015) <= 5e-4);
  printf("  T_unit(chem) vs T_unit(rad) rel err = %.2e ; T0 = %.4f K  %s\n",
         t_consistency, T0, g3 ? "PASS" : "FAIL");
  fail += !g3;

  // (4) RSLA / radiation-constant sanity.
  printf("--- GATE 4: c_code, arad_code magnitudes ---\n");
  const double c_code = c_light / U.v;
  const double ar = arad_code(U);
  const bool g4 = std::isfinite(c_code) && c_code > 1.5e6 && c_code < 1.7e6 &&
                  std::isfinite(ar) && ar > 0.0;
  printf("  c_code = %.4e (expect ~1.578e6) ; arad_code = %.4e  %s\n",
         c_code, ar, g4 ? "PASS" : "FAIL");
  fail += !g4;

  // (5) DRIFT SENTINEL (informational): defaults vs exact BE (fiducial IC mass=6,temp=10,f=5).
  printf("--- SENTINEL: PhysicalUnits defaults vs exact BE normalization ---\n");
  const double msun = 1.9891e33, bemass = 197.561, cs10 = 1.9e4, G = 6.67259e-8;
  const double mass = 6.0, temp = 10.0, f = 5.0;
  const double m0 = mass * msun / (bemass * f);
  const double v0 = cs10 * std::sqrt(temp / 10.0);
  const double rho0 = std::pow(v0, 6) / (m0 * m0) / (64.0 * PI * PI * PI * G * G * G);
  const double t0 = 1.0 / std::sqrt(4.0 * PI * G * rho0);
  const double l0 = v0 * t0;
  const Scales be = Derive(rho0, v0, l0, 2.29);
  printf("  %-14s %14s %14s %10s\n", "unit", "exact BE", "defaults", "rel diff");
  auto row = [&](const char *n, double e, double d) {
    printf("  %-14s %14.6e %14.6e %10.2e\n", n, e, d, relerr(d, e));
  };
  row("diffusivity", be.diffusivity, U.diffusivity);
  row("opacity", be.opacity, U.opacity);
  row("magnetic", be.magnetic, U.magnetic);
  row("temperature", be.temperature, U.temperature);
  printf("  NOTE: this column is the HISTORICAL drift (rounded <units> defaults vs exact BE).\n");
  printf("        Phase 1 Option A: BuildPhysicalUnits now DERIVES the base scales from the\n");
  printf("        BE IC for collapse_be, so the live code path uses the 'exact BE' column and\n");
  printf("        the drift is 0. The defaults remain only as the non-collapse fallback.\n");

  printf("\n%s (%d gate failures)\n", fail == 0 ? "ALL GATES PASS" : "GATE FAILURES", fail);
  return fail == 0 ? 0 : 1;
}
