//========================================================================================
// AthenaPK GPU chemistry package.
//
// Greenfield, device-callable reaction network operating on AthenaPK's existing
// passive scalars (the species ARE the scalars: <hydro> nscalars carries the
// abundances x_i = n_i / n_H as scalar_density_i / rho). The advection of the
// scalars is handled by the hydro integrator; this package adds ONLY the
// operator-split reaction source term, run on the final stage after RT.
//
// Gated behind <physics> chemistry = false (default) => zero impact on existing
// runs. Mirrors the structure of the M1 Radiation package (Initialize +
// AddRadiationTasks).
//========================================================================================
#ifndef CHEMISTRY_CHEMISTRY_HPP_
#define CHEMISTRY_CHEMISTRY_HPP_

#include <memory>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

using namespace parthenon::package::prelude;
using namespace parthenon::driver::prelude;

namespace Chemistry {

// Package factory — called from Hydro::ProcessPackages when <physics> chemistry=true.
std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin);

// Operator-split per-cell reaction update on the conserved passive scalars
// (scalar_density_i). Reads x_i = scalar_density_i / rho, integrates the network
// over the hydro step dt at fixed rho, writes scalar_density_i = rho * x_i.
TaskStatus ReactScalars(MeshData<Real> *md, const Real dt);

// Operator-split chemistry entry point. Mirrors Radiation::AddRadiationTasks;
// called from HydroDriver::MakeTaskCollection on the final stage only.
void AddChemistryTasks(TaskCollection &tc, Mesh *pmesh, const Real dt);

} // namespace Chemistry

#endif // CHEMISTRY_CHEMISTRY_HPP_
