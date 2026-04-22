//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Self-gravity package, ported from Artemis (LANL).
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
#ifndef SELF_GRAVITY_SELF_GRAVITY_HPP_
#define SELF_GRAVITY_SELF_GRAVITY_HPP_

#include <memory>
#include <string>

#include <parthenon/package.hpp>

using namespace parthenon::package::prelude;

#include <parthenon/driver.hpp>
using namespace parthenon::driver::prelude;

namespace SelfGravity {

// Field types — use the VARIABLE macro pattern from Parthenon's poisson_package.hpp
// "grav.phi" and "grav.rhs" are the HDF5 field names you'll see in outputs.
#define SG_VARIABLE(ns, varname)                                                         \
  struct varname : public parthenon::variable_names::base_t<false> {                     \
    template <class... Ts>                                                               \
    KOKKOS_INLINE_FUNCTION varname(Ts &&...args)                                         \
        : parthenon::variable_names::base_t<false>(std::forward<Ts>(args)...) {}         \
    static std::string name() { return #ns "." #varname; }                               \
  }

namespace grav {
SG_VARIABLE(grav, phi);
SG_VARIABLE(grav, rhs);
} // namespace grav

// Package registration — called from Hydro::ProcessPackages
std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin);

// FillDerived function: computes rhs = 4piG * (rho - rho_mean) on every cell
// including ghosts (ghosts so the solver's boundary logic sees the right values).
void FillPoissonRHS(MeshData<Real> *md);

// Apply gravitational acceleration to momentum and energy using phi.
// Flux-weighted energy update (Artemis style) for better AMR energy conservation.
TaskStatus ApplyGravitySource(MeshData<Real> *md, const parthenon::SimTime &tm,
                              const Real beta_dt);

// Build and submit the Poisson solve into the task collection.
// Called from HydroDriver::MakeTaskCollection on the final stage.
void SolvePoisson(TaskCollection &tc, Mesh *pmesh);

} // namespace SelfGravity

#endif // SELF_GRAVITY_SELF_GRAVITY_HPP_
