//========================================================================================
// AthenaPK - single authoritative unit system for the shared FHC normalization.
//
// FLAGSHIP PHASE 1 (audit Workstream A infra): one immutable set of base + derived unit
// scales, built ONCE before any physics package and consumed by ALL of them (radiation,
// chemistry, non-ideal-MHD ionization, and the collapse problem generator). It replaces
// the per-package hardcoded calibrations that were the root of audit findings #3 and #4:
//   #3 magnetic unit dropped the Heaviside-Lorentz sqrt(4*pi) in the collapse sidecar;
//   #4 radiation used a temperature unit (mu=2.33 -> T0~10.19 K) inconsistent with the
//      chemistry/diffusion T0~10.015 K.
//
// The FHC codes are normalized by (rho_unit, v_unit, length_unit) with four_pi_G = 1
// (see CLAUDE.md "Code units are shared"). Everything else is DERIVED here so there is no
// second, drifting copy:
//   time_unit          = length/v
//   mass_unit          = rho * length^3
//   energy_density_unit= rho * v^2                       (== pressure unit)
//   magnetic_unit      = sqrt(4*pi*rho) * v              (Heaviside-Lorentz; via Units)
//   temperature_unit   = mu_thermal * m_H * v^2 / k_B    (c_s,iso = 1 at the reference T)
//   diffusivity_unit   = time/length^2                   (eta_code = eta_cgs*diffusivity_unit)
//   opacity_unit       = rho * length                    (kappa_code = kappa_cgs*opacity_unit)
//
// TEMPERATURE / audit #4:  mu_thermal (default 2.29) is the mean molecular weight IMPLIED
// by the shared normalization (c_s = 1.9e4 cm/s at 10 K); it gives T0 = 10.015 K and is
// the SINGLE temperature calibration. It is deliberately DISTINCT from the neutral mean
// molecular weight mu_n = 2.33 (n = rho/(mu_n m_H)), which is composition and now lives in
// IonizationEnvironment, not here -- conflating the two was audit finding #4.
//========================================================================================
#ifndef UNITS_PHYSICAL_UNITS_HPP_
#define UNITS_PHYSICAL_UNITS_HPP_

#include <cmath>

#include <parameter_input.hpp>
#include <parthenon/package.hpp>

#include "basic_types.hpp"
#include "../units.hpp"        // canonical Units: the sqrt(4*pi) magnetic-unit definition (#3)
#include "be_normalization.hpp" // Phase 1 Option A: shared BE base-scale derivation

namespace PhysUnits {

using parthenon::Real;

// CGS physical constants used to derive the temperature/radiation units. Single copy
// consumed by every package (radiation/chemistry/diffusion previously each carried their
// own; the values here match radiation_opacity.hpp so the migration is bit-identical).
namespace cgs {
constexpr Real m_H = 1.6726219e-24;     // hydrogen mass [g]
constexpr Real k_B = 1.380649e-16;      // Boltzmann [erg/K]
constexpr Real a_rad = 7.5657e-15;      // radiation constant a_R [erg cm^-3 K^-4]
constexpr Real c_light = 2.99792458e10; // speed of light [cm/s]
} // namespace cgs

//----------------------------------------------------------------------------------------
//! Immutable unit system. POD (trivially copyable) so it can be captured by value into
//! device kernels or stashed as a package Param and read identically everywhere.
struct PhysicalUnits {
  // --- base scales (cgs per code unit) ---
  Real rho_unit = 5.467e-19;    // g/cm^3
  Real v_unit = 1.9e4;          // cm/s
  Real length_unit = 2.81e16;   // cm
  Real time_unit = 0.0;         // s      (= length/v)
  Real mass_unit = 0.0;         // g      (= rho*length^3)
  // --- derived scales ---
  Real energy_density_unit = 0.0; // erg/cm^3 (= rho*v^2, the pressure unit)
  Real magnetic_unit = 0.0;       // G       (= sqrt(4*pi*rho)*v, Heaviside-Lorentz)
  Real temperature_unit = 0.0;    // K per code T
  Real diffusivity_unit = 0.0;    // eta_code = eta_cgs*diffusivity_unit (= time/length^2)
  Real opacity_unit = 0.0;        // kappa_code = kappa_cgs*opacity_unit (= rho*length)
  // --- calibration ---
  Real mu_thermal = 2.29;         // mean mol. wt setting temperature_unit (NOT mu_n)

  // Radiation constant in code units: E_eq = arad_code * T_code^4.
  KOKKOS_INLINE_FUNCTION Real arad_code() const {
    const Real T2 = temperature_unit * temperature_unit;
    return cgs::a_rad * T2 * T2 / energy_density_unit;
  }
  // Speed of light in code units (v_unit); RSLA is mandatory for the FHC (~1.58e6).
  KOKKOS_INLINE_FUNCTION Real c_code() const { return cgs::c_light / v_unit; }
};

//----------------------------------------------------------------------------------------
//! Build the authoritative unit system from the input deck. Deterministic in `pin`, so
//! calling it from several packages yields the SAME object (a shared definition even when
//! not a shared instance -- the packages are Initialize()d independently). Base scales are
//! read from the <units> block with the shared-FHC defaults; every other scale is derived.
inline PhysicalUnits BuildPhysicalUnits(parthenon::ParameterInput *pin) {
  const std::string bn = "units";
  PhysicalUnits U;
  // Phase 1 Option A: for the BE collapse problem the base scales are DERIVED from the
  // physical IC via the SAME normalization the pgen uses (single source of truth), so the
  // microphysics units match the dynamics exactly (kills the ~0.13% rounded-defaults drift).
  // Any other problem uses the <units> block (defaults = the historical FHC scales). An
  // explicit <units> base-scale override for a BE run is REJECTED, not silently honored: it
  // would desync the microphysics cgs mapping from the IC-fixed dynamics -- exactly the drift
  // this consolidation removes. (mu_thermal is still read from <units>; it is a calibration,
  // not a base scale.)
  const std::string pb = "problem/collapse_be";
  const bool have_be = pin->DoesParameterExist(pb, "mass") &&
                       pin->DoesParameterExist(pb, "temperature") &&
                       pin->DoesParameterExist(pb, "f");
  if (have_be) {
    PARTHENON_REQUIRE(
        !(pin->DoesParameterExist(bn, "rho_unit_cgs") ||
          pin->DoesParameterExist(bn, "v_unit_cgs") ||
          pin->DoesParameterExist(bn, "length_unit_cgs")),
        "<units> rho_unit_cgs/v_unit_cgs/length_unit_cgs override is not allowed for the "
        "collapse_be problem: the base scales are fixed by the BE IC (mass, temperature, f) "
        "so the microphysics units stay in sync with the dynamics. Remove the override.");
    const auto be = DeriveBENormalization(pin->GetReal(pb, "mass"),
                                          pin->GetReal(pb, "temperature"),
                                          pin->GetReal(pb, "f"));
    U.rho_unit = be.rho0;
    U.v_unit = be.v0;
    U.length_unit = be.l0;
  } else {
    U.rho_unit = pin->GetOrAddReal(bn, "rho_unit_cgs", 5.467e-19);
    U.v_unit = pin->GetOrAddReal(bn, "v_unit_cgs", 1.9e4);
    U.length_unit = pin->GetOrAddReal(bn, "length_unit_cgs", 2.81e16);
  }
  U.mu_thermal = pin->GetOrAddReal(bn, "mu_thermal", 2.29);

  U.time_unit = U.length_unit / U.v_unit;
  U.mass_unit = U.rho_unit * U.length_unit * U.length_unit * U.length_unit;
  U.energy_density_unit = U.rho_unit * U.v_unit * U.v_unit;
  // Route the magnetic unit through the canonical Units definition (audit #3 single source).
  U.magnetic_unit = Units::CodeMagneticCgs(U.mass_unit, U.length_unit, U.time_unit);
  U.temperature_unit = U.mu_thermal * cgs::m_H * U.v_unit * U.v_unit / cgs::k_B;
  U.diffusivity_unit = U.time_unit / (U.length_unit * U.length_unit);
  U.opacity_unit = U.rho_unit * U.length_unit;
  return U;
}

} // namespace PhysUnits

#endif // UNITS_PHYSICAL_UNITS_HPP_
