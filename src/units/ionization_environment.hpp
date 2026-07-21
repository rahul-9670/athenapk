//========================================================================================
// AthenaPK - shared ionization environment (cosmic-ray rate + composition).
//
// FLAGSHIP PHASE 1: a single owner of the ionization parameters that MUST agree between
// the chemistry network (which evolves x_e that ambipolar_coeff=ionization_chem consumes)
// and the equilibrium ionization model used by Ohmic/Hall/ambipolar_coeff=ionization.
// This GENERALIZES audit fix #1 from an init-time consistency ASSERTION (buried in
// chemistry.cpp) into a shared object: both consumers build it from the same input deck
// and read the same cosmic-ray rate and neutral composition, so a silent sqrt(10)-class
// mismatch in the charge population is structurally impossible.
//
// Ownership split vs PhysicalUnits: PhysicalUnits owns *units* (incl. mu_thermal, the
// temperature calibration). IonizationEnvironment owns the *charge microphysics inputs*:
// the cosmic-ray ionization rate and the neutral mean molecular weight mu_n used for
// number densities n = rho/(mu_n m_H). (Keeping mu_thermal and mu_n in separate objects
// is the structural fix for audit finding #4's conflation of the two.)
//
// Room to grow (audit Workstream A.2 / Phase 3): CR attenuation with column, radionuclide
// floor, grain population, thermal-ionization params can move here so every non-ideal
// coefficient and the chemistry network draw one composition.
//========================================================================================
#ifndef UNITS_IONIZATION_ENVIRONMENT_HPP_
#define UNITS_IONIZATION_ENVIRONMENT_HPP_

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include <parameter_input.hpp>
#include <parthenon/package.hpp>

#include "basic_types.hpp"

namespace PhysUnits {

using parthenon::Real;

//----------------------------------------------------------------------------------------
//! Shared charge-microphysics environment. POD, trivially copyable.
struct IonizationEnvironment {
  Real zeta_cr_cgs = 1.0e-16; // cosmic-ray ionization rate [s^-1] (single source of truth)
  Real mu_n = 2.33;           // neutral mean molecular weight (n = rho/(mu_n m_H))
};

//----------------------------------------------------------------------------------------
//! Build the shared environment from the input deck. Deterministic in `pin` so chemistry
//! and the diffusion ionization model obtain identical values from independent Initialize()
//! calls. The cosmic-ray rate is owned by chemistry when it is enabled (it evolves the x_e
//! that ambipolar diffusion consumes); otherwise by the diffusion equilibrium model.
inline IonizationEnvironment BuildIonizationEnvironment(parthenon::ParameterInput *pin) {
  IonizationEnvironment env;
  const bool have_chem = pin->GetOrAddBoolean("physics", "chemistry", false);

  // Canonical cosmic-ray rate.
  if (have_chem) {
    env.zeta_cr_cgs = pin->GetOrAddReal("chemistry", "zeta_cr_cgs", 1.0e-16);
    env.mu_n = pin->GetOrAddReal("chemistry", "mu_n", 2.33);
  } else {
    env.zeta_cr_cgs = pin->GetOrAddReal("diffusion", "ion_zeta", 1.0e-17);
    env.mu_n = pin->GetOrAddReal("diffusion", "ion_mu_n", 2.33);
  }

  // Consistency guard (audit #1): if BOTH the chemistry network AND a diffusion equilibrium
  // ionization coefficient are configured, the two cosmic-ray rates must describe the same
  // charge population. The diffusion rate defaults to the (authoritative) chemistry rate,
  // so an unset ion_zeta auto-agrees; a user-set, disagreeing ion_zeta is a hard error.
  auto uses_ion_model = [&](const std::string &key) {
    const std::string c = pin->GetOrAddString("diffusion", key, "none");
    return c == "ionization" || c == "ionization_chem";
  };
  const bool ion_nonideal = uses_ion_model("resistivity_coeff") ||
                            uses_ion_model("hall_coeff") ||
                            uses_ion_model("ambipolar_coeff");
  if (have_chem && ion_nonideal) {
    const Real ion_zeta = pin->GetOrAddReal("diffusion", "ion_zeta", env.zeta_cr_cgs);
    const bool ok = std::abs(ion_zeta - env.zeta_cr_cgs) <=
                    1.0e-6 * std::max(ion_zeta, env.zeta_cr_cgs);
    char msg[512];
    std::snprintf(msg, sizeof(msg),
                  "Inconsistent cosmic-ray ionization rate: <chemistry> zeta_cr_cgs=%.3e "
                  "s^-1 but <diffusion> ion_zeta=%.3e s^-1. The chemistry x_e "
                  "(ambipolar_coeff=ionization_chem) and the equilibrium ionization model "
                  "(resistivity/hall/ambipolar_coeff=ionization) must describe the same "
                  "charge population. Set <diffusion> ion_zeta = <chemistry> zeta_cr_cgs "
                  "(or leave it unset to inherit it).",
                  env.zeta_cr_cgs, ion_zeta);
    PARTHENON_REQUIRE(ok, std::string(msg));
  }
  return env;
}

} // namespace PhysUnits

#endif // UNITS_IONIZATION_ENVIRONMENT_HPP_
