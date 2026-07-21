//========================================================================================
// AthenaPK WS-4 dust package interface. Mirrors the Chemistry package: evolves the two
// dust passive scalars (f_dg, a_c) with an operator-split growth/sublimation source, gated
// behind <dust> evolve=false (default) => zero impact on existing runs.
//========================================================================================
#ifndef DUST_DUST_PKG_HPP_
#define DUST_DUST_PKG_HPP_

#include <memory>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

using namespace parthenon::package::prelude;
using namespace parthenon::driver::prelude;

namespace Dust {

// Package factory — called from Hydro::ProcessPackages when <physics> dust=true.
std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin);

// Operator-split per-cell dust update on the conserved passive scalars carrying
// (f_dg, a_c). Reads them, integrates growth+sublimation over dt at fixed gas state,
// writes them back. a_c is stored in code length units internally? No -- stored in cm
// directly (scalar just advects the value); f_dg is dimensionless.
TaskStatus ReactDust(MeshData<Real> *md, const Real dt);

// Adds the dust reaction task; called from HydroDriver on the final stage after chemistry.
void AddDustTasks(TaskCollection &tc, Mesh *pmesh, const Real dt);

} // namespace Dust

#endif // DUST_DUST_PKG_HPP_
