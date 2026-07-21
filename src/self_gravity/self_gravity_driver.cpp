//========================================================================================
// AthenaPK - Self-gravity solver driver
// Adapted from Artemis (LANL) and parthenon poisson_gmg example.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include <memory>
#include <string>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <solvers/solver_utils.hpp>

#include "self_gravity.hpp"
#include "poisson_equation.hpp"
#include <solvers/internal_prolongation.hpp>
#include <solvers/mg_solver.hpp>
#include <solvers/bicgstab_solver.hpp>
#include <solvers/solver_base.hpp>

using PoissEq = SelfGravity::PoissonEquation<SelfGravity::grav::phi>;
using prolongator_t = parthenon::solvers::ProlongationBlockInteriorZeroDirichlet;
using preconditioner_t = parthenon::solvers::MGSolver<PoissEq, prolongator_t>;
using SolverT = parthenon::solvers::BiCGSTABSolver<PoissEq, preconditioner_t>;

namespace SelfGravity {

// NOTE (resolved): Multi-rank GPU runs of the GMG-preconditioned BiCGSTAB solver
// once required CUDA_LAUNCH_BLOCKING=1 at runtime to avoid an asynchronous-execution
// race in Parthenon's boundary exchange (a sender posted its MPI send before the
// device-side buffer-pack kernel had completed, so a neighbor rank received garbage
// ghost data and the preconditioned solve diverged to grav.phi -> NaN).
// This is now FIXED at the source by a Kokkos::fence() before the MPI send in
// external/parthenon/src/bvals/comms/boundary_communication.cpp (~line 158, plus an
// unpack-side fence ~line 358). CUDA_LAUNCH_BLOCKING is therefore NO LONGER NEEDED
// and should be removed from submit scripts (it serializes every kernel and badly
// throttles the many small multigrid kernels). Verified: launch-blocking-off multi-
// rank GPU runs are bit-identical to the old launch-blocking baseline (no NaN).
void SolvePoisson(TaskCollection &tc, Mesh *pmesh) {
  using namespace parthenon;
  TaskID none(0);

  auto pkg = pmesh->packages.Get("self_gravity");
  auto psolver = pkg->Param<std::shared_ptr<parthenon::solvers::SolverBase>>("solver_pointer");

  auto partitions = pmesh->GetDefaultBlockPartitions();
  const int num_partitions = partitions.size();
  TaskRegion &region = tc.AddRegion(num_partitions);
  for (int i = 0; i < num_partitions; ++i) {
    TaskList &tl = region[i];

    auto &md = pmesh->mesh_data.Add("base", partitions[i]);
    auto &md_phi = pmesh->mesh_data.Add("phi", md, {grav::phi::name()});
    auto &md_rhs = pmesh->mesh_data.Add("rhs", md, {grav::phi::name()});

    // The solver expects both "phi" container and "rhs" container to hold
    // fields named grav::phi (it operates on IndependentVars = {grav::phi}).
    // rhs lives in field grav::rhs in md. Copy into grav::phi slot of md_rhs.
    auto copy_rhs = tl.AddTask(
        none, TF(parthenon::solvers::utils::between_fields::CopyData<grav::rhs, grav::phi>),
        md);
    copy_rhs = tl.AddTask(
        copy_rhs,
        TF(parthenon::solvers::utils::CopyData<parthenon::TypeList<grav::phi>>),
        md, md_rhs);

    // Solve
    auto setup = psolver->AddSetupTasks(tl, copy_rhs, i, pmesh);
    auto solve = psolver->AddTasks(tl, setup, i, pmesh);

    // Communicate phi ghost cells after solve (so ApplyGravitySource sees full halo).
    auto bcs = parthenon::AddBoundaryExchangeTasks(solve, tl, md_phi, pmesh->multilevel);

    // Copy solution from md_phi's grav::phi slot back into md (the base container)
    // so downstream tasks (ApplyGravitySource) and outputs can find it.
    tl.AddTask(
        bcs,
        TF(parthenon::solvers::utils::CopyData<parthenon::TypeList<grav::phi>>),
        md_phi, md);
  }
}

} // namespace SelfGravity
