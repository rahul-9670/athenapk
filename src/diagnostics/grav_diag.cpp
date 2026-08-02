//========================================================================================
// AthenaPK - self-gravity (Poisson) solver convergence diagnostics. See grav_diag.hpp.
//========================================================================================

#include <iostream>
#include <string>

#include <parthenon/package.hpp>
#include <solvers/solver_base.hpp>

#include "grav_diag.hpp"

namespace Diagnostics {

using parthenon::Real;
using parthenon::TaskStatus;

namespace {
// Param keys used to hand the recorded scalars from the post-solve task to the history
// reductions. Mutable, and deliberately NOT Restart mutability: they describe the most
// recent solve, so a restart legitimately starts them at "no solve yet" (-1).
constexpr char kIters[] = "grav_last_iters";
constexpr char kRes[] = "grav_last_res";
constexpr char kNonconv[] = "grav_last_nonconv";
constexpr char kNonconvCount[] = "grav_nonconv_count";
} // namespace

Real GravDiagReport(MeshData<Real> *md, GravDiag which) {
  auto pkg = md->GetParentPointer()->packages.Get("self_gravity");
  switch (which) {
  case GravDiag::iters:
    return pkg->Param<Real>(kIters);
  case GravDiag::res:
    return pkg->Param<Real>(kRes);
  case GravDiag::nonconv:
    return pkg->Param<Real>(kNonconv);
  }
  return -1.0;
}

TaskStatus GravDiagRecord(Mesh *pmesh) {
  auto pkg = pmesh->packages.Get("self_gravity");
  auto psolver =
      pkg->Param<std::shared_ptr<parthenon::solvers::SolverBase>>("solver_pointer");

  // GetFinalIterations() returns final_iteration + 1 (see solver_base.hpp), i.e. the
  // number of iterations actually taken.
  const int iters = psolver->GetFinalIterations();
  const Real res = psolver->GetFinalResidual();

  // B2 (2026-08-02): ask the solver directly instead of inferring from the iteration count.
  // The old test `iters >= max_iters` was a heuristic and had a real false-positive mode -- a
  // solve that converges exactly ON the last permitted iteration met its tolerance but was
  // reported as a bail-out. SolverBase::GetConverged() now records which branch actually fired.
  const bool bailed = !psolver->GetConverged();

  pkg->UpdateParam<Real>(kIters, static_cast<Real>(iters));
  pkg->UpdateParam<Real>(kRes, res);
  pkg->UpdateParam<Real>(kNonconv, bailed ? 1.0 : 0.0);
  return TaskStatus::complete;
}

//! Unconditional non-convergence warning -- runs whether or not `solver_diag` is on.
//! Reads solver state and a package Param only; touches no field data, so it cannot perturb
//! results. See the rationale block in self_gravity.cpp where the Params are registered.
TaskStatus GravConvergenceWarn(Mesh *pmesh) {
  auto pkg = pmesh->packages.Get("self_gravity");
  auto psolver =
      pkg->Param<std::shared_ptr<parthenon::solvers::SolverBase>>("solver_pointer");
  if (psolver->GetConverged()) return TaskStatus::complete;

  const int n = pkg->Param<int>(kNonconvCount) + 1;
  pkg->UpdateParam<int>(kNonconvCount, n);
  // Warn on the first occurrences and then geometrically, so a persistently non-converging
  // run neither hides nor drowns the log.
  if (parthenon::Globals::my_rank == 0 &&
      (n <= 4 || n == 10 || n == 100 || n == 1000 || n % 10000 == 0)) {
    std::cout << "## WARNING [self-gravity] Poisson solve did NOT converge: used all "
              << psolver->GetFinalIterations()
              << " of max_iterations=" << pkg->Param<int>("grav_max_iters")
              << " with rms residual " << psolver->GetFinalResidual()
              << " (tolerance not met). Occurrence #" << n
              << ". The gravitational potential this step is only as good as that "
                 "residual -- see WP-5/B2 in VALIDATION_PLAN.md."
              << std::endl;
  }
  return TaskStatus::complete;
}

} // namespace Diagnostics
