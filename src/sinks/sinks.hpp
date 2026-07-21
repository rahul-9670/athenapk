//========================================================================================
// AthenaPK - sink-particle package (WS-1 of PHYSICS_COMPLETION_PLAN.md)
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
#ifndef SINKS_SINKS_HPP_
#define SINKS_SINKS_HPP_

#include <memory>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

#include "../main.hpp"
#include "basic_types.hpp"

using namespace parthenon::driver::prelude;
using namespace parthenon::package::prelude;

namespace Sinks {

// Registers the "sinks" package. Increment 1 (swarm plumbing): a Parthenon swarm carrying
// {mass, v, L, t_created, id} that is advected ballistically, migrates across blocks/ranks,
// and survives restart. Creation/accretion/gravity coupling are later WS-1 increments.
// Default OFF (sinks/enabled=false) -> registers only the flag, no swarm, bit-identical.
std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin);

// UserWorkBeforeLoopMesh: seed the initial sink(s) from <sinks> input params (skipped on
// restart, where Parthenon restores the swarm from the restart file).
void SeedInitialSinks(Mesh *pmesh, ParameterInput *pin, parthenon::SimTime &tm);

// Mesh-level sink advance: gather all sinks, integrate the N-body system with a subcycled KDK
// leapfrog (WS-1 inc 3; reduces to exact drift for a single/force-free sink = inc 1), write
// each rank's sinks back, and update neighbor-block indices for SwarmContainer migration.
TaskStatus AdvanceSinksNBody(MeshData<Real> *md, const Real dt);

// WS-1 increment 2: gather every active sink's {x,y,z,mass} onto all ranks, then apply the
// sinks' softened point-mass gravity to the gas as an operator-split source.
TaskStatus GatherSinks(MeshData<Real> *md);
TaskStatus ApplySinkGravity(MeshData<Real> *md, const parthenon::SimTime &tm, const Real dt);

// WS-1 increment 4: spawn one sink per step at the densest cell meeting all formation criteria.
TaskStatus CreateSinks(MeshData<Real> *md, const parthenon::SimTime &tm);

// WS-1 increment 5: conservatively accrete bound gas within r_acc onto the sinks (mass +
// momentum conserving; caps the sink-region density at rho_sink/3 so the timestep recovers).
TaskStatus AccreteSinks(MeshData<Real> *md, const parthenon::SimTime &tm, const Real dt);

} // namespace Sinks

#endif // SINKS_SINKS_HPP_
