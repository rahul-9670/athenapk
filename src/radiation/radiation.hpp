//========================================================================================
// AthenaPK - M1 moment radiation transport package.
// Ported from Artemis (LANL, BSD) src/radiation/moments. Licensed BSD 3-Clause.
//
// INCREMENT 1 (scaffold): registers the M1 two-moment fields (Er, Fr1..3) and
// physical-constant/opacity params, gated behind <physics> radiation = false.
// Transport (CalculateFluxes/FluxSource), M1 closure, and matter coupling land in
// later increments. No tasks are submitted yet, so enabling this is a no-op on
// dynamics until increment 2.
//========================================================================================
#ifndef RADIATION_RADIATION_HPP_
#define RADIATION_RADIATION_HPP_

#include <memory>
#include <string>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

using namespace parthenon::package::prelude;
using namespace parthenon::driver::prelude;

namespace Radiation {

// Field-name macro (mirrors SelfGravity's SG_VARIABLE). HDF5 field names are
// "rad.Er", "rad.Fr1", ... Scalar (base_t<false>) for the scaffold; the flux may
// be consolidated into a single 3-vector field in a later increment.
#define RAD_VARIABLE(varname)                                                            \
  struct varname : public parthenon::variable_names::base_t<false> {                     \
    template <class... Ts>                                                               \
    KOKKOS_INLINE_FUNCTION varname(Ts &&...args)                                         \
        : parthenon::variable_names::base_t<false>(std::forward<Ts>(args)...) {}         \
    static std::string name() { return "rad." #varname; }                                \
  }

namespace rad {
RAD_VARIABLE(Er);   // radiation energy density (lab frame), code units
RAD_VARIABLE(Fr1);  // radiation flux, x1-component
RAD_VARIABLE(Fr2);  // radiation flux, x2-component
RAD_VARIABLE(Fr3);  // radiation flux, x3-component
} // namespace rad

// Package factory — called from Hydro::ProcessPackages when <physics> radiation=true.
std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin);

// --- Increment 2b: explicit M1 transport (free-streaming + diffusion) ---------
// Donor-cell M1 HLL fluxes for (Er, Fr1, Fr2, Fr3) into the WithFluxes flux arrays.
TaskStatus CalculateRadFluxes(MeshData<Real> *md);
// Explicit FV update U += dt * (-div F), then causal clamp (Er>=floor, |F|<=chat*Er).
TaskStatus ApplyRadUpdate(MeshData<Real> *md, const Real dt);
// Radiation CFL timestep over the whole mesh: cfl * dx_min / chat (reduced c).
Real RadDtMesh(Mesh *pmesh);

// --- Increment 3: matter coupling (RT owns the gas thermal energy) ------------
// Implicit (per-cell Newton) gray absorption/emission + flux--momentum exchange
// between the M1 moments and the gas (cons). Ported from Artemis "Simple" coupling
// (radiation/moments/matter_coupling.hpp). Replaces the barotropic EOS floor:
// when this package is active, collapse_be skips its e_th overwrite so RT sets T.
TaskStatus MatterCoupling(MeshData<Real> *md, const Real dt);

// Operator-split M1 update entry point. Mirrors SelfGravity::SolvePoisson; called
// from HydroDriver::MakeTaskCollection. Sub-cycles the transport at the radiation CFL.
void AddRadiationTasks(TaskCollection &tc, Mesh *pmesh, const Real dt);

} // namespace Radiation

#endif // RADIATION_RADIATION_HPP_
