//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2020-2025, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Parthenon headers
#include "amr_criteria/refinement_package.hpp"
#include "basic_types.hpp"
#include "bvals/comms/bvals_in_one.hpp"
#include "prolong_restrict/prolong_restrict.hpp"
#include "utils/error_checking.hpp"
#include <parthenon/parthenon.hpp>
// AthenaPK headers
#include "../eos/adiabatic_hydro.hpp"
#include "../pgen/cluster/agn_triggering.hpp"
#include "../pgen/cluster/magnetic_tower.hpp"
#include "../tracers/tracers.hpp"
#include "../sinks/sinks.hpp"
#include "ct/ct.hpp"
#include "diffusion/diffusion.hpp"
#include "glmmhd/glmmhd.hpp"
#include "hydro.hpp"
#include "hydro_driver.hpp"
#include "../self_gravity/self_gravity.hpp"
#include "../radiation/radiation.hpp"
#include "../chemistry/chemistry.hpp"
#include "../dust/dust_pkg.hpp"
#include "../pgen/pgen.hpp"

using namespace parthenon::driver::prelude;

namespace Hydro {

HydroDriver::HydroDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pm)
    : MultiStageDriver(pin, app_in, pm) {
  // fail if these are not specified in the input file
  pin->CheckRequired("hydro", "eos");

  // warn if these fields aren't specified in the input file
  pin->CheckDesired("parthenon/time", "cfl");

  // WP-13 (2026-08-02): fail LOUDLY if the collapse_be Params are missing.
  //
  // MakeTaskCollection adds ApplyBarotropicCooling only `if hasKey("collapse_be_rhocrit")`.
  // That key used to be registered exclusively in ProblemGenerator -- which Parthenon does
  // NOT call on restart -- so on every resume the key was absent, the task was quietly never
  // added, and the run silently lost the outside-sphere momentum BC (and, on a non-RT deck,
  // the barotropic cooling too). No error, no warning; it was found only by diffing a
  // straight run against a restarted one. Registration now lives in
  // collapse_be::ProblemInitPackageData, which runs on both paths. This assert exists so that
  // if that wiring is ever broken again the run dies at startup instead of silently changing
  // physics mid-chain -- which matters because every production submit script self-resumes.
  // Read WITHOUT registering a default: another call site already registers job/problem_id
  // with its own default, and Parthenon aborts if two defaults for one key disagree.
  const std::string problem_id = pin->DoesParameterExist("job", "problem_id")
                                     ? pin->GetString("job", "problem_id")
                                     : std::string("");
  if (pm != nullptr && problem_id == "collapse_be") {
    auto hydro_pkg = pm->packages.Get("Hydro");
    PARTHENON_REQUIRE(
        hydro_pkg->AllParams().hasKey("collapse_be_rhocrit"),
        "collapse_be: Param 'collapse_be_rhocrit' is missing, so ApplyBarotropicCooling -- "
        "and with it the outside-sphere momentum BC -- would be SILENTLY skipped. Is "
        "collapse_be::ProblemInitPackageData still wired up in main.cpp? See WP-13 "
        "(docs/validation/WP13_restart_reproducibility.md).");
  }
}

// Sets all fluxes to 0
TaskStatus ResetFluxes(MeshData<Real> *md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  // Pack "cons" by NAME rather than by the {Independent} flag: the STS resets the diffusive
  // fluxes, which are deposited ONLY onto cons. Under Constrained Transport the {Independent}
  // flag also matches the FACE field Bf (WithFluxes) whose EDGE fluxes (E1/E2/E3) have no
  // X1/X2/X3 cell-flux slots -- indexing them here segfaults. Bf's edge EMF is zeroed
  // separately by CT_ZeroEMF. cons-by-name is also the v6-STS-bug-hardened choice.
  const std::vector<std::string> cons_name{"cons"};
  auto cons_pack = md->PackVariablesAndFluxes(cons_name, cons_name);

  const int ndim = pmb->pmy_mesh->ndim;
  // Using separate loops for each dim as the launch overhead should be hidden
  // by enough work over the entire pack and it allows to not use any conditionals.
  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "ResetFluxes X1", parthenon::DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, 0, cons_pack.GetDim(4) - 1, kb.s, kb.e, jb.s, jb.e, ib.s,
      ib.e + 1,
      KOKKOS_LAMBDA(const int b, const int v, const int k, const int j, const int i) {
        auto &cons = cons_pack(b);
        cons.flux(X1DIR, v, k, j, i) = 0.0;
      });

  if (ndim < 2) {
    return TaskStatus::complete;
  }
  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "ResetFluxes X2", parthenon::DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, 0, cons_pack.GetDim(4) - 1, kb.s, kb.e, jb.s, jb.e + 1,
      ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int v, const int k, const int j, const int i) {
        auto &cons = cons_pack(b);
        cons.flux(X2DIR, v, k, j, i) = 0.0;
      });

  if (ndim < 3) {
    return TaskStatus::complete;
  }
  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "ResetFluxes X3", parthenon::DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, 0, cons_pack.GetDim(4) - 1, kb.s, kb.e + 1, jb.s, jb.e,
      ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int v, const int k, const int j, const int i) {
        auto &cons = cons_pack(b);
        cons.flux(X3DIR, v, k, j, i) = 0.0;
      });
  return TaskStatus::complete;
}

TaskStatus RKL2StepFirst(MeshData<Real> *md_Y0, MeshData<Real> *md_Yjm1,
                         MeshData<Real> *md_Yjm2, MeshData<Real> *md_MY0, const int s_rkl,
                         const Real tau) {
  auto pmb = md_Y0->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  // Compute coefficients. Meyer+2014 eq. (18)
  Real mu_tilde_1 = 4. / 3. /
                    (static_cast<Real>(s_rkl) * static_cast<Real>(s_rkl) +
                     static_cast<Real>(s_rkl) - 2.);

  // Pack by NAME ("cons"), not by {Metadata::Independent}: the diffusion terms update
  // only the gas conserved variables, but the Independent flag also matches rad.Er/Fr
  // (radiation) and grav.phi. Since Yjm1 aliases "base" and the Y0 snapshot ("u1") only
  // deep-copies cons/prim, a flag-based pack overwrites the radiation fields and phi in
  // base with the u1 container's stale/zero content on every STS call (this silently
  // zeroed rad.Er for the whole run past cycle 71000, see DEV_LOG 2026-07-15).
  const std::vector<std::string> sts_vars{"cons"};
  auto Y0 = md_Y0->PackVariablesAndFluxes(sts_vars, sts_vars);
  auto Yjm1 = md_Yjm1->PackVariablesAndFluxes(sts_vars, sts_vars);
  auto Yjm2 = md_Yjm2->PackVariablesAndFluxes(sts_vars, sts_vars);
  auto MY0 = md_MY0->PackVariablesAndFluxes(sts_vars, sts_vars);

  // Using separate loops for each dim as the launch overhead should be hidden
  // by enough work over the entire pack and it allows to not use any conditionals.
  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "RKL first step", parthenon::DevExecSpace(), 0,
      Y0.GetDim(5) - 1, 0, Y0.GetDim(4) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int v, const int k, const int j, const int i) {
        Yjm1(b, v, k, j, i) =
            Y0(b, v, k, j, i) + mu_tilde_1 * tau * MY0(b, v, k, j, i); // Y_1
        Yjm2(b, v, k, j, i) = Y0(b, v, k, j, i);                       // Y_0
      });

  return TaskStatus::complete;
}

TaskStatus RKL2StepOther(MeshData<Real> *md_Y0, MeshData<Real> *md_Yjm1,
                         MeshData<Real> *md_Yjm2, MeshData<Real> *md_MY0, const Real mu_j,
                         const Real nu_j, const Real mu_tilde_j, const Real gamma_tilde_j,
                         const Real tau) {
  auto pmb = md_Y0->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  // Pack by NAME ("cons") — same constraint as RKL2StepFirst above: a flag-based
  // {Independent} pack pulls in rad.Er/Fr and grav.phi and corrupts them via the
  // cons/prim-only u1 snapshot.
  const std::vector<std::string> sts_vars{"cons"};
  auto Y0 = md_Y0->PackVariablesAndFluxes(sts_vars, sts_vars);
  auto Yjm1 = md_Yjm1->PackVariablesAndFluxes(sts_vars, sts_vars);
  auto Yjm2 = md_Yjm2->PackVariablesAndFluxes(sts_vars, sts_vars);
  auto MY0 = md_MY0->PackVariablesAndFluxes(sts_vars, sts_vars);

  const int ndim = pmb->pmy_mesh->ndim;
  // Using separate loops for each dim as the launch overhead should be hidden
  // by enough work over the entire pack and it allows to not use any conditionals.
  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "RKL other step", parthenon::DevExecSpace(), 0,
      Y0.GetDim(5) - 1, 0, Y0.GetDim(4) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int v, const int k, const int j, const int i) {
        // First calc this step
        const auto &coords = Yjm1.GetCoords(b);
        const Real MYjm1 =
            parthenon::Update::FluxDivHelper(v, k, j, i, ndim, coords, Yjm1(b));
        const Real Yj = mu_j * Yjm1(b, v, k, j, i) + nu_j * Yjm2(b, v, k, j, i) +
                        (1.0 - mu_j - nu_j) * Y0(b, v, k, j, i) +
                        mu_tilde_j * tau * MYjm1 +
                        gamma_tilde_j * tau * MY0(b, v, k, j, i);
        // Then shuffle vars for next step
        Yjm2(b, v, k, j, i) = Yjm1(b, v, k, j, i);
        Yjm1(b, v, k, j, i) = Yj;
      });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
// TEMPORARY DIAGNOSTICS (2026-07-29, CT-vs-GLM ambipolar under-diffusion investigation).
// Not gated behind a flag; each is an extra blocking global reduction + rank-0 print, added
// as ordinary Task nodes (same pattern as MagneticTowerReducePowerContribs above) so they
// run at actual task-graph EXECUTION time (not graph-build time, when "base" would still
// hold the previous stage's stale data). Intended to be reverted once the investigation
// concludes; harmless to leave in (pure diagnostics, no state mutation).

// eta_A cache stats (max, mean) over the interior. Same code path/cache for CT and GLM, so
// this establishes whether the ambipolar diffusivity ITSELF differs between the two schemes
// at a given (identical, restart-derived) state.
TaskStatus DiagEtaA(MeshData<Real> *md, const std::string label) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto hydro_pkg_diag = pmb->packages.Get("Hydro");
  if (!hydro_pkg_diag->Param<bool>("nonideal_eta_cache")) {
    if (parthenon::Globals::my_rank == 0) {
      std::cout << "DIAG_ETA_A " << label << " (cache off, skipped)" << std::endl;
    }
    return TaskStatus::complete;
  }
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
  auto eta_pack = md->PackVariables(std::vector<std::string>{"nonideal_eta"});
  const int nb = eta_pack.GetDim(5);
  Real lmax = 0.0, lsum = 0.0;
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "DiagEtaA_max", parthenon::DevExecSpace(), 0,
      nb - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &m) {
        m = std::max(m, eta_pack(b, NonidealEtaIdx::A, k, j, i));
      },
      Kokkos::Max<Real>(lmax));
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "DiagEtaA_sum", parthenon::DevExecSpace(), 0,
      nb - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &s) {
        s += eta_pack(b, NonidealEtaIdx::A, k, j, i);
      },
      Kokkos::Sum<Real>(lsum));
  Real lcount = static_cast<Real>(nb) * (kb.e - kb.s + 1) * (jb.e - jb.s + 1) *
               (ib.e - ib.s + 1);
#ifdef MPI_PARALLEL
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, &lmax, 1, MPI_PARTHENON_REAL, MPI_MAX, MPI_COMM_WORLD));
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, &lsum, 1, MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD));
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, &lcount, 1, MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD));
#endif
  if (parthenon::Globals::my_rank == 0) {
    std::cout << "DIAG_ETA_A " << label << " max=" << lmax
              << " mean=" << (lsum / std::max(lcount, 1.0)) << " ncell=" << lcount
              << std::endl;
  }
  return TaskStatus::complete;
}

// Stage-init M(Y0) magnitude restricted to the B-components, GLM path: reads MY0's "cons"
// IB1..IB3 (the flux divergence of AmbipolarDiffFluxIsoFixed + OhmicDiffFluxIsoFixed, etc.,
// evaluated on the Y0/restart state).
TaskStatus DiagMOpB_GLM(MeshData<Real> *md_MY0, const std::string label) {
  auto pmb = md_MY0->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
  auto MY0 = md_MY0->PackVariables(std::vector<std::string>{"cons"});
  const int nb = MY0.GetDim(5);
  Real lsum = 0.0, lmax = 0.0;
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "DiagMOpB_GLM_sum", parthenon::DevExecSpace(), 0,
      nb - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &s) {
        s += std::abs(MY0(b, IB1, k, j, i)) + std::abs(MY0(b, IB1 + 1, k, j, i)) +
             std::abs(MY0(b, IB1 + 2, k, j, i));
      },
      Kokkos::Sum<Real>(lsum));
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "DiagMOpB_GLM_max", parthenon::DevExecSpace(), 0,
      nb - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &m) {
        const Real v = std::abs(MY0(b, IB1, k, j, i)) + std::abs(MY0(b, IB1 + 1, k, j, i)) +
                       std::abs(MY0(b, IB1 + 2, k, j, i));
        m = std::max(m, v);
      },
      Kokkos::Max<Real>(lmax));
#ifdef MPI_PARALLEL
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, &lsum, 1, MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD));
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, &lmax, 1, MPI_PARTHENON_REAL, MPI_MAX, MPI_COMM_WORLD));
#endif
  if (parthenon::Globals::my_rank == 0) {
    std::cout << "DIAG_MOP_B " << label << " L1sum=" << lsum << " Linf=" << lmax
              << std::endl;
  }
  return TaskStatus::complete;
}

// Stage-init M(Y0) magnitude, CT path: reads MY0's face field Bf (already holds
// -curl(E_diff) from CT_CurlEMFToBf at the point this is called).
TaskStatus DiagMOpB_CT(MeshData<Real> *md_MY0, const std::string label) {
  const int ndim = md_MY0->GetMeshPointer()->ndim;
  const bool three_d = ndim > 2;
  static auto desc = parthenon::MakePackDescriptor<Hydro::CT::Bf>(md_MY0);
  auto MY0 = desc.GetPack(md_MY0);
  const int nb = MY0.GetNBlocks();
  using TE = parthenon::TopologicalElement;
  using parthenon::CellLevel;
  Real lsum = 0.0, lmax = 0.0;
  auto accum = [&](TE te) {
    IndexRange ib = md_MY0->GetBoundsI(CellLevel::same, IndexDomain::interior, te);
    IndexRange jb = md_MY0->GetBoundsJ(CellLevel::same, IndexDomain::interior, te);
    IndexRange kb = md_MY0->GetBoundsK(CellLevel::same, IndexDomain::interior, te);
    Real s = 0.0, m = 0.0;
    parthenon::par_reduce(
        parthenon::loop_pattern_mdrange_tag, "DiagMOpB_CT_sum", parthenon::DevExecSpace(), 0,
        nb - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &acc) {
          acc += std::abs(MY0(b, te, Hydro::CT::Bf(), k, j, i));
        },
        Kokkos::Sum<Real>(s));
    parthenon::par_reduce(
        parthenon::loop_pattern_mdrange_tag, "DiagMOpB_CT_max", parthenon::DevExecSpace(), 0,
        nb - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &acc) {
          acc = std::max(acc, std::abs(MY0(b, te, Hydro::CT::Bf(), k, j, i)));
        },
        Kokkos::Max<Real>(m));
    lsum += s;
    lmax = std::max(lmax, m);
  };
  accum(TE::F1);
  accum(TE::F2);
  if (three_d) accum(TE::F3);
#ifdef MPI_PARALLEL
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, &lsum, 1, MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD));
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, &lmax, 1, MPI_PARTHENON_REAL, MPI_MAX, MPI_COMM_WORLD));
#endif
  if (parthenon::Globals::my_rank == 0) {
    std::cout << "DIAG_MOP_B " << label << " L1sum=" << lsum << " Linf=" << lmax
              << std::endl;
  }
  return TaskStatus::complete;
}

// Total magnetic energy density integral proxy: sum(B^2) over interior cells from "cons"
// IB1..IB3 -- valid/comparable for both CT (post CT_ProjectBfToCC) and GLM at any task
// boundary. Called before/after/during the STS group to track the actual suppression
// achieved, apples-to-apples.
TaskStatus DiagTotalB2(MeshData<Real> *md, const std::string label) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
  auto cons = md->PackVariables(std::vector<std::string>{"cons"});
  const int nb = cons.GetDim(5);
  Real lsum = 0.0;
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "DiagTotalB2", parthenon::DevExecSpace(), 0,
      nb - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &s) {
        s += SQR(cons(b, IB1, k, j, i)) + SQR(cons(b, IB1 + 1, k, j, i)) +
             SQR(cons(b, IB1 + 2, k, j, i));
      },
      Kokkos::Sum<Real>(lsum));
#ifdef MPI_PARALLEL
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, &lsum, 1, MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD));
#endif
  if (parthenon::Globals::my_rank == 0) {
    std::cout << "DIAG_TOTAL_B2 " << label << " sumB2=" << lsum << std::endl;
  }
  return TaskStatus::complete;
}

// PER-LEVEL magnetic-energy dissipation rate integral sum_V V * B . M(B), where M(B)=dB/dt of
// the (Ohmic+AD) parabolic operator evaluated on Y0. A discretely DISSIPATIVE operator must
// give sum V B.M <= 0 (modulo boundary terms); a positive/near-zero value flags a
// NON-dissipative (spurious) spatial operator. Split by AMR level to distinguish a coarse-fine
// (edge-restriction) defect [non-dissipation localized to blocks at level boundaries] from a
// local edge/curl/projection defect [non-dissipation throughout, incl. fine-block interiors].
// CT variant: M = MY0.Bf (the -curl(E_diff) face field after CT_CurlEMFToBf) projected to cell
// centres; B = Y0 (base) cons IB1..3. Reports rate = (sum V B.M)/(sum V B^2) per level.
TaskStatus DiagBdotM_CT(MeshData<Real> *md_base, MeshData<Real> *md_MY0,
                        const std::string label) {
  const bool three_d = md_base->GetMeshPointer()->ndim > 2;
  auto cons = md_base->PackVariables(std::vector<std::string>{"cons"});
  static auto desc = parthenon::MakePackDescriptor<Hydro::CT::Bf>(md_MY0);
  auto MY0 = desc.GetPack(md_MY0);
  static auto desc_b = parthenon::MakePackDescriptor<Hydro::CT::Bf>(md_base);
  auto Bf0 = desc_b.GetPack(md_base); // Y0 face field (base.Bf = u1.Bf snapshot)
  using TE = parthenon::TopologicalElement;
  constexpr int NLEV = 24;
  double vb2[NLEV] = {0.0}, vbm[NLEV] = {0.0};
  // FACE-NATIVE (no cell projection): dot co-located face Bf with face M(Bf). The projected
  // measure above 0.5-averages Bf[i]+Bf[i+1] -> NULLS the Nyquist/grid-scale mode; this one
  // keeps it, so it reveals AD acting on the sharp face field that the projection can't see.
  double fb2[NLEV] = {0.0}, fbm[NLEV] = {0.0};
  const int nblk = md_base->NumBlocks();
  for (int bb = 0; bb < nblk; ++bb) {
    auto pmb = md_base->GetBlockData(bb)->GetBlockPointer();
    const int lev = pmb->loc.level();
    IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
    IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
    IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
    Real lb2 = 0.0, lbm = 0.0, lf2 = 0.0, lfm = 0.0;
    parthenon::par_reduce(
        parthenon::loop_pattern_mdrange_tag, "DiagBdotM_CT", parthenon::DevExecSpace(), kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int k, const int j, const int i, Real &s2, Real &sm, Real &f2,
                      Real &fm) {
          const auto &co = cons.GetCoords(bb);
          const Real V = co.CellVolume(k, j, i);
          const Real bx = cons(bb, IB1, k, j, i), by = cons(bb, IB2, k, j, i),
                     bz = cons(bb, IB3, k, j, i);
          const Real Mx = 0.5 * (MY0(bb, TE::F1, Hydro::CT::Bf(), k, j, i) +
                                 MY0(bb, TE::F1, Hydro::CT::Bf(), k, j, i + 1));
          const Real My = 0.5 * (MY0(bb, TE::F2, Hydro::CT::Bf(), k, j, i) +
                                 MY0(bb, TE::F2, Hydro::CT::Bf(), k, j + 1, i));
          const Real Mz = three_d ? 0.5 * (MY0(bb, TE::F3, Hydro::CT::Bf(), k, j, i) +
                                           MY0(bb, TE::F3, Hydro::CT::Bf(), k + 1, j, i))
                                  : MY0(bb, TE::F3, Hydro::CT::Bf(), k, j, i);
          s2 += V * (bx * bx + by * by + bz * bz);
          sm += V * (bx * Mx + by * My + bz * Mz);
          // face-native (lower faces of the cell; co-located Bf.M, no averaging)
          const Real f1b = Bf0(bb, TE::F1, Hydro::CT::Bf(), k, j, i);
          const Real f2b = Bf0(bb, TE::F2, Hydro::CT::Bf(), k, j, i);
          const Real f3b = three_d ? Bf0(bb, TE::F3, Hydro::CT::Bf(), k, j, i) : 0.0;
          const Real f1m = MY0(bb, TE::F1, Hydro::CT::Bf(), k, j, i);
          const Real f2m = MY0(bb, TE::F2, Hydro::CT::Bf(), k, j, i);
          const Real f3m = three_d ? MY0(bb, TE::F3, Hydro::CT::Bf(), k, j, i) : 0.0;
          f2 += V * (f1b * f1b + f2b * f2b + f3b * f3b);
          fm += V * (f1b * f1m + f2b * f2m + f3b * f3m);
        },
        Kokkos::Sum<Real>(lb2), Kokkos::Sum<Real>(lbm), Kokkos::Sum<Real>(lf2),
        Kokkos::Sum<Real>(lfm));
    if (lev >= 0 && lev < NLEV) {
      vb2[lev] += lb2;
      vbm[lev] += lbm;
      fb2[lev] += lf2;
      fbm[lev] += lfm;
    }
  }
#ifdef MPI_PARALLEL
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, vb2, NLEV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, vbm, NLEV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, fb2, NLEV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, fbm, NLEV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
#endif
  if (parthenon::Globals::my_rank == 0) {
    for (int L = 0; L < NLEV; ++L)
      if (fb2[L] != 0.0)
        std::cout << "DIAG_BDOTM_FACE " << label << " L" << L << " FB2=" << fb2[L]
                  << " FBdotM=" << fbm[L] << " rate=" << (fbm[L] / fb2[L]) << std::endl;
  }
  if (parthenon::Globals::my_rank == 0) {
    double t2 = 0.0, tm = 0.0;
    for (int L = 0; L < NLEV; ++L) {
      t2 += vb2[L];
      tm += vbm[L];
      if (vb2[L] != 0.0)
        std::cout << "DIAG_BDOTM " << label << " L" << L << " VB2=" << vb2[L]
                  << " VBdotM=" << vbm[L] << " rate=" << (vbm[L] / vb2[L]) << std::endl;
    }
    std::cout << "DIAG_BDOTM " << label << " TOTAL VB2=" << t2 << " VBdotM=" << tm
              << " rate=" << (tm / t2) << " (rate<0 dissipative)" << std::endl;
  }
  return TaskStatus::complete;
}

// GLM variant of DiagBdotM: M = MY0 cons IB1..3 (the FluxDivergence of the parabolic fluxes);
// B = Y0 (base) cons IB1..3. Same per-level sum V B.M reference (should be strongly negative).
TaskStatus DiagBdotM_GLM(MeshData<Real> *md_base, MeshData<Real> *md_MY0,
                         const std::string label) {
  auto cons = md_base->PackVariables(std::vector<std::string>{"cons"});
  auto M = md_MY0->PackVariables(std::vector<std::string>{"cons"});
  constexpr int NLEV = 24;
  double vb2[NLEV] = {0.0}, vbm[NLEV] = {0.0};
  const int nblk = md_base->NumBlocks();
  for (int bb = 0; bb < nblk; ++bb) {
    auto pmb = md_base->GetBlockData(bb)->GetBlockPointer();
    const int lev = pmb->loc.level();
    IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
    IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
    IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
    Real lb2 = 0.0, lbm = 0.0;
    parthenon::par_reduce(
        parthenon::loop_pattern_mdrange_tag, "DiagBdotM_GLM", parthenon::DevExecSpace(), kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int k, const int j, const int i, Real &s2, Real &sm) {
          const auto &co = cons.GetCoords(bb);
          const Real V = co.CellVolume(k, j, i);
          const Real bx = cons(bb, IB1, k, j, i), by = cons(bb, IB2, k, j, i),
                     bz = cons(bb, IB3, k, j, i);
          s2 += V * (bx * bx + by * by + bz * bz);
          sm += V * (bx * M(bb, IB1, k, j, i) + by * M(bb, IB2, k, j, i) +
                     bz * M(bb, IB3, k, j, i));
        },
        Kokkos::Sum<Real>(lb2), Kokkos::Sum<Real>(lbm));
    if (lev >= 0 && lev < NLEV) {
      vb2[lev] += lb2;
      vbm[lev] += lbm;
    }
  }
#ifdef MPI_PARALLEL
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, vb2, NLEV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  PARTHENON_MPI_CHECK(
      MPI_Allreduce(MPI_IN_PLACE, vbm, NLEV, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
#endif
  if (parthenon::Globals::my_rank == 0) {
    double t2 = 0.0, tm = 0.0;
    for (int L = 0; L < NLEV; ++L) {
      t2 += vb2[L];
      tm += vbm[L];
      if (vb2[L] != 0.0)
        std::cout << "DIAG_BDOTM " << label << " L" << L << " VB2=" << vb2[L]
                  << " VBdotM=" << vbm[L] << " rate=" << (vbm[L] / vb2[L]) << std::endl;
    }
    std::cout << "DIAG_BDOTM " << label << " TOTAL VB2=" << t2 << " VBdotM=" << tm
              << " rate=" << (tm / t2) << " (rate<0 dissipative)" << std::endl;
  }
  return TaskStatus::complete;
}

// Assumes that prim and cons are in sync initially.
// Guarantees that prim and cons are in sync at the end.
void AddSTSTasks(TaskCollection *ptask_coll, Mesh *pmesh, BlockList_t &blocks,
                 const Real tau) {

  auto hydro_pkg = blocks[0]->packages.Get("Hydro");
  auto mindt_diff = hydro_pkg->Param<Real>("dt_diff");
  // With diffusion/rkl2_freeze_eta the stage-init below refreshes the non-ideal eta
  // cache from the half-start state and the stages jj=2..s reuse it (refill=false).
  const bool refill_eta_stages = !hydro_pkg->Param<bool>("rkl2_freeze_eta");
  // Under Constrained Transport the cons.flux(IBn) diffusive deposits are gated off, so the
  // cons recurrence super-time-steps only the gas energy (Ohmic/AD Poynting) and leaves B
  // inert; the divergence-free induction is super-time-stepped in PARALLEL on the face field
  // Bf with M(Bf) = -curl(E_diff). Each Bf sub-stage is projected onto cons IB1..IB3 so the
  // AD eta / edge EMF of the next sub-stage sees the updated field. See ct.hpp / ct.cpp.
  const bool use_ct = hydro_pkg->Param<bool>("use_ct");
  // TEMPORARY DIAGNOSTICS tag (see DiagEtaA/DiagMOpB_*/DiagTotalB2 above).
  const std::string diag_tag = use_ct ? "CT" : "GLM";
  // TEMPORARY: compile-time gate for the CT-vs-GLM AD diagnostics (per-substep blocking
  // reductions + rank-0 prints). false = production (no overhead, no log spam); flip to true
  // to re-enable the DiagEtaA/DiagMOpB/DiagTotalB2/DiagBdotM instrumentation. Remove with the
  // Diag* definitions once the CT-AD investigation is fully closed.
  constexpr bool kCTDiag = false;

  // get number of RKL steps
  // eq (21) using half hyperbolic timestep due to Strang split
  int s_rkl =
      static_cast<int>(0.5 * (std::sqrt(9.0 + 16.0 * tau / mindt_diff) - 1.0)) + 1;
  // ensure odd number of stages
  if (s_rkl % 2 == 0) s_rkl += 1;

  if (parthenon::Globals::my_rank == 0) {
    const auto ratio = 2.0 * tau / mindt_diff;
    std::cout << "STS ratio: " << ratio << " Taking " << s_rkl << " steps." << std::endl;
    if (ratio > 400.1) {
      std::cout << "WARNING: ratio is > 400. Proceed at own risk." << std::endl;
    }
  }

  TaskID none(0);

  // Store initial u0 in u1 as "base" will continuously be updated but initial state Y0 is
  // required for each stage.
  TaskRegion &region_copy_out = ptask_coll->AddRegion(blocks.size());
  for (int i = 0; i < blocks.size(); i++) {
    auto &tl = region_copy_out[i];
    auto &Y0 = blocks[i]->meshblock_data.Get("u1");
    auto &base = blocks[i]->meshblock_data.Get();
    tl.AddTask(
        none,
        [](MeshBlockData<Real> *dst, MeshBlockData<Real> *src, bool use_ct) {
          dst->Get("cons").data.DeepCopy(src->Get("cons").data);
          dst->Get("prim").data.DeepCopy(src->Get("prim").data);
          // CT: the face field Bf is the primary magnetic variable, super-time-stepped on
          // its own register here, so its Y0 snapshot must be copied too (the cons/prim-only
          // copy leaves u1.Bf stale -> the Bf recurrence would read a wrong Y0).
          if (use_ct) {
            dst->Get("Bf").data.DeepCopy(src->Get("Bf").data);
          }
          return TaskStatus::complete;
        },
        Y0.get(), base.get(), use_ct);
  }

  TaskRegion &region_init = ptask_coll->AddRegion(blocks.size());
  for (int i = 0; i < blocks.size(); i++) {
    auto &pmb = blocks[i];
    auto &tl = region_init[i];
    auto &base = pmb->meshblock_data.Get();

    // Add extra registers. No-op for existing variables so it's safe to call every
    // time.
    // TODO(pgrete) this allocates all Variables, i.e., prim and cons vector, but only a
    // subset is actually needed. Streamline to allocate only required vars.
    pmb->meshblock_data.Add("MY0", base);
    pmb->meshblock_data.Add("Yjm2", base);
  }

  const int num_partitions = pmesh->DefaultNumPartitions();
  TaskRegion &region_calc_fluxes_step_init = ptask_coll->AddRegion(num_partitions);
  for (int i = 0; i < num_partitions; i++) {
    auto &tl = region_calc_fluxes_step_init[i];
    auto &base = pmesh->mesh_data.GetOrAdd("base", i);
    const auto any = parthenon::BoundaryType::any;
    auto start_bnd = tl.AddTask(none, parthenon::StartReceiveBoundBufs<any>, base);
    auto start_flxcor_recv =
        tl.AddTask(none, parthenon::StartReceiveFluxCorrections, base);

    // Reset flux arrays (not guaranteed to be zero)
    auto reset_fluxes = tl.AddTask(none, ResetFluxes, base.get());

    // TEMPORARY DIAGNOSTIC: total B^2 of the Y0/restart-derived state, before this STS
    // group touches anything.
    if (kCTDiag) tl.AddTask(reset_fluxes, DiagTotalB2, base.get(),
              diag_tag + " part" + std::to_string(i) + " pre(Y0)");

    // Calculate the diffusive fluxes for Y0 (here still "base" as nothing has been
    // updated yet) so that we can store the result as MY0 and reuse later
    // (in every subsetp). Parabolic terms only: in the mixed rkl2+Hall mode the
    // dispersive Hall EMF is applied unsplit with the hyperbolic fluxes, not here.
    // Always refreshes the eta cache (refill=true): this is the half-start state the
    // frozen stages reuse.
    auto hydro_diff_fluxes = tl.AddTask(reset_fluxes, CalcDiffFluxes, hydro_pkg.get(),
                                        base.get(), DiffTermSet::parabolic, true);

    // TEMPORARY DIAGNOSTIC: eta_A cache stats right after the stage-init fill.
    if (kCTDiag) tl.AddTask(hydro_diff_fluxes, DiagEtaA, base.get(),
              diag_tag + " part" + std::to_string(i) + " stage-init");

    // CT: build the diffusive-only edge EMF (Ohmic + AD) into base.Bf.flux BEFORE the flux
    // correction, so the coarse-fine reflux restricts it too (as in the unsplit step).
    TaskID diff_ready = hydro_diff_fluxes;
    if (use_ct) {
      auto zero_emf = tl.AddTask(hydro_diff_fluxes, Hydro::CT::CT_ZeroEMF, base.get());
      auto ohm_emf = tl.AddTask(zero_emf, Hydro::CT::CT_AddOhmicEMF, base.get());
      auto ad_emf = tl.AddTask(ohm_emf, Hydro::CT::CT_AddAmbipolarEMF, base.get());
      // Rebuild the diffusive Poynting energy flux from this SAME edge EMF so the energy
      // flux and the induction share one stencil (the face-based deposits in
      // resistivity.cpp/ambipolar.cpp are gated off by ct_edge_poynting). Before the flux
      // correction so the coarse-fine reflux restricts the IEN deposit too.
      diff_ready = tl.AddTask(ad_emf, Hydro::CT::CT_AddDiffusivePoynting, base.get());
    }

    auto send_flx = tl.AddTask(diff_ready, parthenon::LoadAndSendFluxCorrections, base);
    auto recv_flx =
        tl.AddTask(start_flxcor_recv, parthenon::ReceiveFluxCorrections, base);
    auto set_flx =
        tl.AddTask(recv_flx | diff_ready, parthenon::SetFluxCorrections, base);

    auto &Y0 = pmesh->mesh_data.GetOrAdd("u1", i);
    auto &MY0 = pmesh->mesh_data.GetOrAdd("MY0", i);
    auto &Yjm2 = pmesh->mesh_data.GetOrAdd("Yjm2", i);

    auto init_MY0 = tl.AddTask(set_flx, parthenon::Update::FluxDivergence<MeshData<Real>>,
                               base.get(), MY0.get());

    // TEMPORARY DIAGNOSTIC: stage-init M(Y0) restricted to the B-components, cons/GLM path
    // (informative for CT too: should show near-zero since the IBn deposit is gated off).
    if (kCTDiag) tl.AddTask(init_MY0, DiagMOpB_GLM, MY0.get(),
              diag_tag + " part" + std::to_string(i) + " cons-IB");
    // TEMPORARY DIAGNOSTIC: per-level dissipation integral sum_V V B.M(Y0). GLM path (cons-IB
    // M); for CT this shows the inert cons-IB (~0, IBn deposit gated off) -- the live CT
    // measure is DiagBdotM_CT after CT_CurlEMFToBf below.
    if (kCTDiag) tl.AddTask(init_MY0, DiagBdotM_GLM, Y0.get(), MY0.get(),
              diag_tag + " part" + std::to_string(i) + " consIB");

    // Initialize Y0 and Y1 and the recursion relation starting with j = 2 needs data from
    // the two preceeding stages.
    auto rkl2_step_first = tl.AddTask(init_MY0, RKL2StepFirst, Y0.get(), base.get(),
                                      Yjm2.get(), MY0.get(), s_rkl, tau);

    // CT: MY0.Bf = -curl(E_diff of Y0), then the RKL2 first sub-step on Bf, then project
    // the updated face field onto cons IB1..IB3 (overwriting the inert cons-IB recurrence).
    TaskID after_step = rkl2_step_first;
    if (use_ct) {
      auto curl_my0 =
          tl.AddTask(set_flx, Hydro::CT::CT_CurlEMFToBf, base.get(), MY0.get());
      // TEMPORARY DIAGNOSTIC: stage-init M(Y0) on the face field (the actual CT operator).
      if (kCTDiag) tl.AddTask(curl_my0, DiagMOpB_CT, MY0.get(),
                diag_tag + " part" + std::to_string(i) + " Bf");
      // TEMPORARY DIAGNOSTIC: per-level dissipation integral sum_V V B.M(Y0), CT path (M =
      // projected MY0.Bf). rate<0 => dissipative; split by level localizes any non-dissipation
      // to coarse-fine boundaries (restriction defect) vs fine-block interiors (local defect).
      if (kCTDiag) tl.AddTask(curl_my0, DiagBdotM_CT, Y0.get(), MY0.get(),
                diag_tag + " part" + std::to_string(i) + " Bf");
      auto ct_first =
          tl.AddTask(curl_my0 | rkl2_step_first, Hydro::CT::CT_RKL2FirstBf, Y0.get(),
                     base.get(), Yjm2.get(), MY0.get(), s_rkl, tau);
      after_step = tl.AddTask(ct_first, Hydro::CT::CT_ProjectBfToCC, base.get());
    }

    // Update ghost cells of Y1 (as MY1 is calculated for each Y_j).
    // Y1 stored in "base", see rkl2_step_first task.
    // Update ghost cells (local and non local), prolongate and apply bound cond.
    // TODO(someone) experiment with split (local/nonlocal) comms with respect to
    // performance for various tests (static, amr, block sizes) and then decide on the
    // best impl. Go with default call (split local/nonlocal) for now.
    // TODO(pgrete) optimize (in parthenon) to only send subset of updated vars
    auto bounds_exchange = parthenon::AddBoundaryExchangeTasks(
        after_step | start_bnd, tl, base, pmesh->multilevel);

    auto fill_derived_1 = tl.AddTask(
        bounds_exchange, parthenon::Update::FillDerived<MeshData<Real>>, base.get());

    // TEMPORARY DIAGNOSTIC: total B^2 after the RKL2 first sub-stage (Y1).
    if (kCTDiag) tl.AddTask(fill_derived_1, DiagTotalB2, base.get(),
              diag_tag + " part" + std::to_string(i) + " afterY1(s=1/" +
                  std::to_string(s_rkl) + ")");
  }

  // Compute coefficients. Meyer+2012 eq. (16)
  Real b_j = 1. / 3.;
  Real b_jm1 = 1. / 3.;
  Real b_jm2 = 1. / 3.;
  Real w1 = 4. / (static_cast<Real>(s_rkl) * static_cast<Real>(s_rkl) +
                  static_cast<Real>(s_rkl) - 2.);
  Real mu_j, nu_j, j, mu_tilde_j, gamma_tilde_j;

  // RKL loop
  for (int jj = 2; jj <= s_rkl; jj++) {
    j = static_cast<Real>(jj);
    b_j = (j * j + j - 2.0) / (2 * j * (j + 1.0));
    mu_j = (2.0 * j - 1.0) / j * b_j / b_jm1;
    nu_j = -(j - 1.0) / j * b_j / b_jm2;
    mu_tilde_j = mu_j * w1;
    gamma_tilde_j = -(1.0 - b_jm1) * mu_tilde_j; // -a_jm1*mu_tilde_j

    TaskRegion &region_calc_fluxes_step_other = ptask_coll->AddRegion(num_partitions);
    for (int i = 0; i < num_partitions; i++) {
      auto &tl = region_calc_fluxes_step_other[i];
      auto &base = pmesh->mesh_data.GetOrAdd("base", i);

      // Only need boundaries for base as it's the only "active" container exchanging
      // data/fluxes with neighbors. All other containers are passive (i.e., data is only
      // used but not exchanged).
      const auto any = parthenon::BoundaryType::any;
      auto start_bnd = tl.AddTask(none, parthenon::StartReceiveBoundBufs<any>, base);
      auto start_flxcor_recv =
          tl.AddTask(none, parthenon::StartReceiveFluxCorrections, base);

      // Reset flux arrays (not guaranteed to be zero)
      auto reset_fluxes = tl.AddTask(none, ResetFluxes, base.get());

      // Calculate the diffusive fluxes for Yjm1 (here u1); parabolic terms only (the
      // dispersive Hall EMF is never super-time-stepped, see AddSTSTasks stage-init).
      // refill_eta_stages=false (diffusion/rkl2_freeze_eta) reuses the half-start eta.
      auto hydro_diff_fluxes = tl.AddTask(reset_fluxes, CalcDiffFluxes, hydro_pkg.get(),
                                          base.get(), DiffTermSet::parabolic,
                                          refill_eta_stages);

      // CT: rebuild the diffusive edge EMF into base.Bf.flux for the current sub-stage Yjm1.
      TaskID diff_ready = hydro_diff_fluxes;
      if (use_ct) {
        auto zero_emf = tl.AddTask(hydro_diff_fluxes, Hydro::CT::CT_ZeroEMF, base.get());
        auto ohm_emf = tl.AddTask(zero_emf, Hydro::CT::CT_AddOhmicEMF, base.get());
        auto ad_emf = tl.AddTask(ohm_emf, Hydro::CT::CT_AddAmbipolarEMF, base.get());
        // Same stencil-consistent diffusive Poynting flux as the stage-init block above.
        diff_ready = tl.AddTask(ad_emf, Hydro::CT::CT_AddDiffusivePoynting, base.get());
      }

      auto send_flx = tl.AddTask(diff_ready, parthenon::LoadAndSendFluxCorrections, base);
      auto recv_flx =
          tl.AddTask(start_flxcor_recv, parthenon::ReceiveFluxCorrections, base);
      auto set_flx =
          tl.AddTask(recv_flx | diff_ready, parthenon::SetFluxCorrections, base);

      auto &Y0 = pmesh->mesh_data.GetOrAdd("u1", i);
      auto &MY0 = pmesh->mesh_data.GetOrAdd("MY0", i);
      auto &Yjm2 = pmesh->mesh_data.GetOrAdd("Yjm2", i);

      auto rkl2_step_other =
          tl.AddTask(set_flx, RKL2StepOther, Y0.get(), base.get(), Yjm2.get(), MY0.get(),
                     mu_j, nu_j, mu_tilde_j, gamma_tilde_j, tau);

      // CT: RKL2 sub-step on Bf (MYjm1 = -curl(edge EMF in base.Bf.flux), formed inline),
      // then project Bf onto cons IB1..IB3.
      TaskID after_step = rkl2_step_other;
      if (use_ct) {
        auto ct_other = tl.AddTask(set_flx | rkl2_step_other, Hydro::CT::CT_RKL2OtherBf,
                                   Y0.get(), base.get(), Yjm2.get(), MY0.get(), mu_j, nu_j,
                                   mu_tilde_j, gamma_tilde_j, tau);
        after_step = tl.AddTask(ct_other, Hydro::CT::CT_ProjectBfToCC, base.get());
      }

      // update ghost cells of base (currently storing Yj)
      // Update ghost cells (local and non local), prolongate and apply bound cond.
      // TODO(someone) experiment with split (local/nonlocal) comms with respect to
      // performance for various tests (static, amr, block sizes) and then decide on the
      // best impl. Go with default call (split local/nonlocal) for now.
      // TODO(pgrete) optimize (in parthenon) to only send subset of updated vars
      auto bounds_exchange = parthenon::AddBoundaryExchangeTasks(
          after_step | start_bnd, tl, base, pmesh->multilevel);

      auto fill_derived_j = tl.AddTask(
          bounds_exchange, parthenon::Update::FillDerived<MeshData<Real>>, base.get());

      // TEMPORARY DIAGNOSTIC: total B^2 after RKL2 sub-stage jj (only the LAST stage,
      // s_rkl, to keep output compact -- the intermediate curve isn't needed to answer
      // "how much did the whole STS group suppress B").
      if (jj == s_rkl) {
        if (kCTDiag) tl.AddTask(fill_derived_j, DiagTotalB2, base.get(),
                  diag_tag + " part" + std::to_string(i) + " afterFinal(s=" +
                      std::to_string(s_rkl) + "/" + std::to_string(s_rkl) + ")");
      }
    }

    b_jm2 = b_jm1;
    b_jm1 = b_j;
  }
}

// See the advection.hpp declaration for a description of how this function gets called.
TaskCollection HydroDriver::MakeTaskCollection(BlockList_t &blocks, int stage) {
  TaskCollection tc;
  auto hydro_pkg = blocks[0]->packages.Get("Hydro");
  const bool use_ct = hydro_pkg->Param<bool>("use_ct");

  TaskID none(0);
  // Number of task lists that can be executed indepenently and thus *may*
  // be executed in parallel and asynchronous.
  // Being extra verbose here in this example to highlight that this is not
  // required to be 1 or blocks.size() but could also only apply to a subset of blocks.
  auto num_task_lists_executed_independently = blocks.size();

  const int num_partitions = pmesh->DefaultNumPartitions();

  // calculate agn triggering accretion rate
  if ((stage == 1) &&
      hydro_pkg->AllParams().hasKey("agn_triggering_reduce_accretion_rate") &&
      hydro_pkg->Param<bool>("agn_triggering_reduce_accretion_rate")) {

    // need to make sure that there's only one region in order to MPI_reduce to work
    TaskRegion &single_task_region = tc.AddRegion(1);
    auto &tl = single_task_region[0];
    // First globally reset triggering quantities
    auto prev_task =
        tl.AddTask(none, cluster::AGNTriggeringResetTriggering, hydro_pkg.get());

    // Adding one task for each partition. Given that they're all in one task list
    // they'll be executed sequentially. Given that a par_reduce to a host var is
    // blocking it's also save to store the variable in the Params for now.
    for (int i = 0; i < num_partitions; i++) {
      auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);
      auto new_agn_triggering =
          tl.AddTask(prev_task, cluster::AGNTriggeringReduceTriggering, mu0.get(), tm.dt);
      prev_task = new_agn_triggering;
    }
#ifdef MPI_PARALLEL
    auto reduce_agn_triggering =
        tl.AddTask(prev_task, cluster::AGNTriggeringMPIReduceTriggering, hydro_pkg.get());
    prev_task = reduce_agn_triggering;
#endif

    // Remove accreted gas
    for (int i = 0; i < num_partitions; i++) {
      auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);
      auto new_remove_accreted_gas =
          tl.AddTask(prev_task, cluster::AGNTriggeringFinalizeTriggering, mu0.get(), tm);
      prev_task = new_remove_accreted_gas;
    }
  }

  for (int i = 0; i < blocks.size(); i++) {
    auto &pmb = blocks[i];
    // Using "base" as u0, which already exists (and returned by using plain Get())
    auto &u0 = pmb->meshblock_data.Get();

    // Create meshblock data for register u1.
    // This is a noop if u1 already exists.
    // TODO(pgrete) update to derive from other quanity as u1 does not require fluxes
    if (stage == 1) {
      pmb->meshblock_data.Add("u1", u0);
    }
  }

  // calculate magnetic tower scaling
  if ((stage == 1) && hydro_pkg->AllParams().hasKey("magnetic_tower_power_scaling") &&
      hydro_pkg->Param<bool>("magnetic_tower_power_scaling")) {
    const auto &magnetic_tower =
        hydro_pkg->Param<cluster::MagneticTower>("magnetic_tower");

    // need to make sure that there's only one region in order to MPI_reduce to work
    TaskRegion &single_task_region = tc.AddRegion(1);
    auto &tl = single_task_region[0];
    // First globally reset magnetic_tower_linear_contrib and
    // magnetic_tower_quadratic_contrib
    auto prev_task =
        tl.AddTask(none, cluster::MagneticTowerResetPowerContribs, hydro_pkg.get());

    // Adding one task for each partition. Given that they're all in one task list
    // they'll be executed sequentially. Given that a par_reduce to a host var is
    // blocking it's also save to store the variable in the Params for now.
    for (int i = 0; i < num_partitions; i++) {
      auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);
      auto new_magnetic_tower_power_contrib =
          tl.AddTask(prev_task, cluster::MagneticTowerReducePowerContribs, mu0.get(), tm);
      prev_task = new_magnetic_tower_power_contrib;
    }
#ifdef MPI_PARALLEL
    auto reduce_magnetic_tower_power_contrib = tl.AddTask(
        prev_task,
        [](StateDescriptor *hydro_pkg) {
          Real magnetic_tower_contribs[] = {
              hydro_pkg->Param<Real>("magnetic_tower_linear_contrib"),
              hydro_pkg->Param<Real>("magnetic_tower_quadratic_contrib")};
          PARTHENON_MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, magnetic_tower_contribs, 2,
                                            MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD));
          hydro_pkg->UpdateParam("magnetic_tower_linear_contrib",
                                 magnetic_tower_contribs[0]);
          hydro_pkg->UpdateParam("magnetic_tower_quadratic_contrib",
                                 magnetic_tower_contribs[1]);
          return TaskStatus::complete;
        },
        hydro_pkg.get());
#endif
  }

  // First add split sources before the main time integration
  if (stage == 1) {
    // If any tasks modify the conserved variables before this place, then
    // the STS tasks should be updated to not assume prim and cons are in sync.
    const auto &diffint = hydro_pkg->Param<DiffInt>("diffint");
    if (diffint == DiffInt::rkl2) {
      AddSTSTasks(&tc, pmesh, blocks, 0.5 * tm.dt);
    }
    TaskRegion &strang_init_region = tc.AddRegion(num_partitions);
    for (int i = 0; i < num_partitions; i++) {
      auto &tl = strang_init_region[i];
      auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);

      // Add initial Strang split source terms, i.e., a dt/2 update
      // IMPORTANT 1: This task must also update `prim` and `cons` variables so that
      // the source term is applied to all active registers in the flux calculation.
      // IMPORTANT 2: The tasks should work using `cons` variables as input as in the
      // final step, `prim` are not updated yet from the flux calculation.
      tl.AddTask(none, AddSplitSourcesStrang, mu0.get(), tm);
    }
  }

  // Now start the main time integration by resetting the registers
  TaskRegion &async_region_init_int = tc.AddRegion(num_task_lists_executed_independently);
  for (int i = 0; i < blocks.size(); i++) {
    auto &pmb = blocks[i];
    auto &tl = async_region_init_int[i];
    auto &u0 = pmb->meshblock_data.Get();
    // init u1, see (11) in Athena++ method paper
    if (stage == 1) {
      auto &u1 = pmb->meshblock_data.Get("u1");
      const bool use_ct_snap = hydro_pkg->Param<bool>("use_ct");
      auto init_u1 = tl.AddTask(
          none,
          [](MeshBlockData<Real> *u0, MeshBlockData<Real> *u1, bool copy_prim,
             bool use_ct) {
            u1->Get("cons").data.DeepCopy(u0->Get("cons").data);
            if (copy_prim) {
              u1->Get("prim").data.DeepCopy(u0->Get("prim").data);
            }
            // CT: the face field Bf is a separate Independent variable (not part of
            // "cons"), so it needs its own stage-start snapshot for the VL2 combine.
            if (use_ct) {
              u1->Get("Bf").data.DeepCopy(u0->Get("Bf").data);
            }
            return TaskStatus::complete;
          },
          // First order flux correction needs the original prim variables in the
          // during the correction.
          u0.get(), u1.get(), hydro_pkg->Param<bool>("first_order_flux_correct"),
          use_ct_snap);
    }
  }

  // note that task within this region that contains one tasklist per pack
  // could still be executed in parallel
  TaskRegion &single_tasklist_per_pack_region = tc.AddRegion(num_partitions);
  for (int i = 0; i < num_partitions; i++) {
    auto &tl = single_tasklist_per_pack_region[i];
    auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);
    auto &mu1 = pmesh->mesh_data.GetOrAdd("u1", i);

    const auto any = parthenon::BoundaryType::any;
    auto start_bnd = tl.AddTask(none, parthenon::StartReceiveBoundBufs<any>, mu0);
    auto start_flxcor_recv =
        tl.AddTask(none, parthenon::StartReceiveFluxCorrections, mu0);

    const auto flux_str = (stage == 1) ? "flux_first_stage" : "flux_other_stage";
    FluxFun_t *calc_flux_fun = hydro_pkg->Param<FluxFun_t *>(flux_str);
    auto calc_flux = tl.AddTask(none, calc_flux_fun, mu0);

    // TODO(pgrete) figure out what to do about the sources from the first stage
    // that are potentially disregarded when the (m)hd fluxes are corrected in the second
    // stage.
    TaskID first_order_flux_correct = calc_flux;
    if (hydro_pkg->Param<bool>("first_order_flux_correct")) {
      auto *first_order_flux_correct_fun =
          hydro_pkg->Param<FirstOrderFluxCorrectFun_t *>("first_order_flux_correct_fun");
      first_order_flux_correct =
          tl.AddTask(calc_flux, first_order_flux_correct_fun, mu0.get(), mu1.get(),
                     integrator->gam0[stage - 1], integrator->gam1[stage - 1],
                     integrator->beta[stage - 1] * integrator->dt);
    }

    // CT increment 2 (AMR): assemble the edge EMF from the primitive state BEFORE the
    // flux-correction round. The flux-correction buffers of a Face variable carry its
    // edge fluxes (GetFluxCorrectionElements), so the same Load/Receive/Set trio that
    // corrects the cons cell-fluxes also restricts the fine-block edge EMFs onto the
    // coarse neighbour's *shared* edge -- both sides then curl the identical EMF, keeping
    // div B at round-off across coarse-fine boundaries (reflux-curl). Must precede
    // LoadAndSendFluxCorrections; in increment 1 the EMF was assembled after set_flx,
    // which is correct single-level but skips the C-F correction. GLM path: emf = none.
    TaskID emf = none;
    if (use_ct) {
      if (hydro_pkg->Param<bool>("use_ct_gs05")) {
        // GS05 upwind EMF needs the transverse-B + mass face fluxes -> depends on
        // calc_flux; dt sets the contact-upwind switch sharpness (GetWeightForCT).
        emf = tl.AddTask(calc_flux, Hydro::CT::CT_AssembleEMF_GS05, mu0.get(),
                         integrator->beta[stage - 1] * integrator->dt);
      } else {
        emf = tl.AddTask(calc_flux, Hydro::CT::CT_AssembleEMF, mu0.get());
      }
      // CT increment 4 (non-ideal): add the PARABOLIC (Ohmic + ambipolar) diffusive edge
      // EMFs onto the ideal edge EMF, BEFORE the flux-correction round so the coarse-fine
      // reflux restricts them too. No-op unless the respective term is active; the matching
      // cons induction deposits are gated off under CT (no double count); their energy
      // (Poynting) terms stay on the FV flux. ONLY in the unsplit integrator: under RKL2 the
      // parabolic induction is super-time-stepped on Bf in AddSTSTasks, so adding it here
      // too would double-apply the diffusion (over-damping / instability).
      if (hydro_pkg->Param<DiffInt>("diffint") != DiffInt::rkl2) {
        emf = tl.AddTask(emf, Hydro::CT::CT_AddOhmicEMF, mu0.get());
        emf = tl.AddTask(emf, Hydro::CT::CT_AddAmbipolarEMF, mu0.get());
      }
      // CT increment 4 (non-ideal): dispersive Hall edge EMF, same face->edge averaging.
      // No-op unless Hall is active. Unsplit only (RKL2+CT forbidden). Its cons induction
      // deposit in HallDiffFluxIsoFixed is gated off; the Poynting term stays FV.
      emf = tl.AddTask(emf, Hydro::CT::CT_AddHallEMF, mu0.get());
    }
    auto send_flx = tl.AddTask(first_order_flux_correct | emf,
                               parthenon::LoadAndSendFluxCorrections, mu0);
    auto recv_flx = tl.AddTask(start_flxcor_recv, parthenon::ReceiveFluxCorrections, mu0);
    auto set_flx = tl.AddTask(recv_flx | first_order_flux_correct,
                              parthenon::SetFluxCorrections, mu0);

    // compute the divergence of fluxes of conserved variables
    auto update = tl.AddTask(
        set_flx, parthenon::Update::UpdateWithFluxDivergence<MeshData<Real>>, mu0.get(),
        mu1.get(), integrator->gam0[stage - 1], integrator->gam1[stage - 1],
        integrator->beta[stage - 1] * integrator->dt);

    // Add non-operator split source terms (cooling, MHD, problem-defined).
    auto source_unsplit = tl.AddTask(update, AddUnsplitSources, mu0.get(), tm,
                                     integrator->beta[stage - 1] * integrator->dt);
    // BE collapse barotropic cooling (only active if problem_id = collapse_be)
    auto after_cooling = source_unsplit;
    {
      auto hydro_pkg = pmesh->packages.Get("Hydro");
      // WP-13: this gate used to fail SILENTLY on every restart. It is now backed by a
      // startup assert in the HydroDriver constructor (which has `pin`), so a missing key is
      // a hard error rather than a quietly-dropped task. See the note there.
      if (hydro_pkg->AllParams().hasKey("collapse_be_rhocrit")) {
        after_cooling = tl.AddTask(source_unsplit, collapse_be::ApplyBarotropicCooling,
                                   mu0.get(), tm,
                                   integrator->beta[stage - 1] * integrator->dt);
      }
    }
    // Then change the next task's dependency from source_unsplit to after_cooling

    // Constrained Transport (Phase 2). Assemble the edge EMF from the primitive state,
    // advance the face field Bf by its curl (VL2 combine), and project Bf onto the
    // cell-centered IB1..IB3 -- overwriting the flux-divergence update of B done above,
    // so the induction half of the step is divergence-free by construction. The Bf
    // ghost exchange rides along in the AddBoundaryExchangeTasks call below (Bf is
    // FillGhost). GLM path (use_ct=false) is untouched and bit-identical.
    auto after_ct = after_cooling;
    if (use_ct) {
      const Real gam0 = integrator->gam0[stage - 1];
      const Real gam1 = integrator->gam1[stage - 1];
      const Real bdt = integrator->beta[stage - 1] * integrator->dt;
      // Bf.flux now holds the coarse-fine-corrected edge EMF (via set_flx above). Curl
      // it onto the faces (VL2 combine), then project Bf -> cell-centered IB1..IB3.
      // Depend on `emf` explicitly: set_flx only transitively covers the edge EMF through
      // the send/recv flux-correction handshake, which on a single/degenerate block can
      // complete before the *later* EMF-chain tasks (Hall is added after Ohmic/ambipolar),
      // letting CT_UpdateBf curl a pre-Hall edge EMF. Gating on emf guarantees the full
      // ideal+Ohmic+ambipolar+Hall edge EMF is assembled before the curl.
      auto ct_update =
          tl.AddTask(set_flx | update | emf, Hydro::CT::CT_UpdateBf, mu0.get(), mu1.get(),
                     gam0, gam1, bdt);
      after_ct = tl.AddTask(ct_update | after_cooling, Hydro::CT::CT_ProjectBfToCC,
                            mu0.get());
    }

    auto source_split_first_order = after_ct;

    if (stage == integrator->nstages) {
      // Audit fix #7: the Strang split source must depend on after_cooling, NOT
      // source_unsplit. Both cooling and the split sources mutate cons; branching the
      // split chain from source_unsplit left them as unordered siblings (a GPU-async
      // race on the ambient momentum that barotropic cooling zeroes). Single chain:
      // update -> source_unsplit -> (cooling) -> strang -> first_order -> bnd_exchange.
      auto source_split_strang_final =
          tl.AddTask(after_ct, AddSplitSourcesStrang, mu0.get(), tm);
      source_split_first_order =
          tl.AddTask(source_split_strang_final, AddSplitSourcesFirstOrder, mu0.get(), tm);
    }

    parthenon::AddBoundaryExchangeTasks(source_split_first_order | start_bnd, tl, mu0,
                                        pmesh->multilevel);
  }
  // --- END of single_tasklist_per_pack_region loop ---

  // ------------------------------------------------------------------
  // Self-gravity: solve Poisson and apply gravity on the final stage only.
  // Must come BEFORE the boundary exchange above has cleared fluxes...
  // actually, fluxes are held in mu0 and not cleared by the boundary exchange,
  // so ordering is fine as long as we do this before the NEXT stage's calc_flux.
  // ------------------------------------------------------------------
  auto self_gravity_pkg = pmesh->packages.AllPackages().count("self_gravity") > 0
                              ? pmesh->packages.Get("self_gravity")
                              : nullptr;
  // Audit fix #2: self-gravity is applied as a STAGE-CONSISTENT source. Poisson is solved
  // from each stage's updated density and the gravitational kick uses that stage's
  // beta*dt weight (matching the hydro flux update and the other unsplit sources), instead
  // of a single lagged full-dt kick on the final stage only. This restores ~2nd-order time
  // coupling for the VL2 integrator (predictor AND corrector feel gravity), and the
  // flux-weighted gravitational work now uses each stage's contemporaneous mass flux
  // (ApplyGravitySource already takes a per-stage beta_dt). Cost: one extra Poisson solve
  // per step. NOTE: this changes collapse timing/energetics vs the old scheme -- validate
  // against the gravity test hierarchy and re-baseline before quoting as production physics.
  if (self_gravity_pkg != nullptr) {
    // Solve Poisson: this adds its own TaskRegion with num_partitions task lists.
    SelfGravity::SolvePoisson(tc, pmesh);

    // Apply gravitational source term per partition, weighted by this stage's beta*dt.
    TaskRegion &sg_source_region = tc.AddRegion(num_partitions);
    for (int i = 0; i < num_partitions; i++) {
      auto &tl = sg_source_region[i];
      auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);
      tl.AddTask(none, SelfGravity::ApplyGravitySource, mu0.get(), tm,
                 integrator->beta[stage - 1] * integrator->dt);
    }
  }

  // Sink-particle gravity on the gas (WS-1 increment 2): gather the global sink list, then
  // add each sink's softened point-mass force as an operator-split source on the final stage.
  if (stage == integrator->nstages &&
      pmesh->packages.Get("sinks")->Param<bool>("enabled")) {
    TaskRegion &sk_grav_region = tc.AddRegion(num_partitions);
    for (int i = 0; i < num_partitions; i++) {
      auto &tl = sk_grav_region[i];
      auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);
      auto gather = tl.AddTask(none, Sinks::GatherSinks, mu0.get());
      auto grav = tl.AddTask(gather, Sinks::ApplySinkGravity, mu0.get(), tm, integrator->dt);
      // Creation runs after the gather (needs the existing-sink list for the 2*r_acc test);
      // accretion runs after creation so a just-formed sink can accrete this step.
      auto create = tl.AddTask(grav, Sinks::CreateSinks, mu0.get(), tm);
      tl.AddTask(create, Sinks::AccreteSinks, mu0.get(), tm, integrator->dt);
    }
  }

  // Radiation (M1): operator-split, sub-cycled transport on the final stage only.
  auto radiation_pkg = pmesh->packages.AllPackages().count("radiation") > 0
                           ? pmesh->packages.Get("radiation")
                           : nullptr;
  if (stage == integrator->nstages && radiation_pkg != nullptr) {
    // Adds its own TaskRegion(s); sub-cycles at the radiation CFL over the hydro dt.
    Radiation::AddRadiationTasks(tc, pmesh, integrator->dt);
  }

  // Chemistry: operator-split reaction source on the passive scalars, final stage
  // only. Runs after RT (which sets the gas temperature the rates depend on) and
  // before FillDerived (so the primitive scalars pick up the reacted abundances).
  auto chemistry_pkg = pmesh->packages.AllPackages().count("chemistry") > 0
                           ? pmesh->packages.Get("chemistry")
                           : nullptr;
  if (stage == integrator->nstages && chemistry_pkg != nullptr) {
    Chemistry::AddChemistryTasks(tc, pmesh, integrator->dt);
  }
  auto dust_pkg = pmesh->packages.AllPackages().count("dust") > 0
                      ? pmesh->packages.Get("dust")
                      : nullptr;
  if (stage == integrator->nstages && dust_pkg != nullptr) {
    Dust::AddDustTasks(tc, pmesh, integrator->dt);
  }

  TaskRegion &single_tasklist_per_pack_region_3 = tc.AddRegion(num_partitions);
  // Audit fix A3 (+ #2 stage-consistent gravity): operator-split sources modify INTERIOR
  // conserved variables AFTER the stage boundary exchange above -- self-gravity on EVERY
  // stage now (#2), and RT matter coupling + chemistry on the final stage. Refresh ghosts
  // before FillDerived so the primitive ghosts (notably the x_e scalar read one cell into
  // the ghost zone by the first RKL2 non-ideal flux, and the gravity-updated momentum used
  // by the next stage's reconstruction) reflect the sourced interior, not stale ghosts.
  const bool sourced_interior =
      (stage == integrator->nstages) || (self_gravity_pkg != nullptr);
  for (int i = 0; i < num_partitions; i++) {
    auto &tl = single_tasklist_per_pack_region_3[i];
    auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);
    TaskID pre = none;
    if (sourced_interior)
      pre = parthenon::AddBoundaryExchangeTasks(none, tl, mu0, pmesh->multilevel);
    auto fill_derived =
        tl.AddTask(pre, parthenon::Update::FillDerived<MeshData<Real>>, mu0.get());
  }
  const auto &diffint = hydro_pkg->Param<DiffInt>("diffint");
  // If any tasks modify the conserved variables before this place and after FillDerived,
  // then the STS tasks should be updated to not assume prim and cons are in sync.
  if (diffint == DiffInt::rkl2 && stage == integrator->nstages) {
    AddSTSTasks(&tc, pmesh, blocks, 0.5 * tm.dt);
  }

  // Single task in single (serial) region to reset global vars used in reductions in the
  // first stage.
  // TODO(pgrete) check if we logically need this reset or if we can reset within the
  // timestep task
  if (stage == integrator->nstages &&
      (hydro_pkg->Param<bool>("calc_c_h") ||
       hydro_pkg->Param<DiffInt>("diffint") != DiffInt::none)) {
    TaskRegion &reset_reduction_vars_region = tc.AddRegion(1);
    auto &tl = reset_reduction_vars_region[0];
    tl.AddTask(
        none,
        [](StateDescriptor *hydro_pkg) {
          hydro_pkg->UpdateParam("mindx", std::numeric_limits<Real>::max());
          hydro_pkg->UpdateParam("dt_hyp", std::numeric_limits<Real>::max());
          hydro_pkg->UpdateParam("dt_diff", std::numeric_limits<Real>::max());
          return TaskStatus::complete;
        },
        hydro_pkg.get());
  }

  if (stage == integrator->nstages) {
    TaskRegion &tr = tc.AddRegion(num_partitions);
    for (int i = 0; i < num_partitions; i++) {
      auto &tl = tr[i];
      auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);
      auto new_dt = tl.AddTask(none, parthenon::Update::EstimateTimestep<MeshData<Real>>,
                               mu0.get());
    }
  }

  auto tracers_pkg = pmesh->packages.Get("tracers");
  // First order operator split tracer advection
  if (stage == integrator->nstages && tracers_pkg->Param<bool>("enabled")) {
    const std::string swarm_name = "tracers";
    TaskRegion &sync_region_tr = tc.AddRegion(1);
    {
      for (auto &pmb : blocks) {
        auto &tl = sync_region_tr[0];
        auto &sd = pmb->meshblock_data.Get()->GetSwarmData();
        auto reset_comms =
            tl.AddTask(none, &SwarmContainer::ResetCommunication, sd.get());
      }
    }

    TaskRegion &async_region_tr = tc.AddRegion(blocks.size());
    for (int n = 0; n < blocks.size(); n++) {
      auto &tl = async_region_tr[n];
      auto &pmb = blocks[n];
      auto &sd = pmb->meshblock_data.Get()->GetSwarmData();
      auto &mbd0 = pmb->meshblock_data.Get("base");
      auto tracer_advect =
          tl.AddTask(none, Tracers::AdvectTracers, mbd0.get(), integrator->dt);

      auto send = tl.AddTask(tracer_advect, &SwarmContainer::Send, sd.get(),
                             BoundaryCommSubset::all);

      auto receive =
          tl.AddTask(send, &SwarmContainer::Receive, sd.get(), BoundaryCommSubset::all);
    }
    // TODO(pgrete) Fix/cleanup once we got swarm packs.
    // We need just a single region with a single task in order to be able to use plain
    // MPI reductions (rather than Parthenon provided reduction tasks that work with
    // arbitrary packs).
    PARTHENON_REQUIRE_THROWS(num_partitions == 1,
                             "Only packs_per_rank=1 currently supported for tracers.")
    TaskRegion &single_tasklist_per_pack_region_4 = tc.AddRegion(num_partitions);
    for (int i = 0; i < num_partitions; i++) {
      auto &tl = single_tasklist_per_pack_region_4[i];
      auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);
      auto fill = tl.AddTask(none, Tracers::FillTracers, mu0.get(), tm);
      if (Tracers::ProblemFillTracers != nullptr) {
        fill =
            tl.AddTask(fill, Tracers::ProblemFillTracers, mu0.get(), tm, integrator->dt);
      }
    }
  }

  // First-order operator-split sink-particle advance (WS-1 increment 1: inert ballistic
  // drift + swarm migration). Mirrors the tracer block above.
  auto sinks_pkg = pmesh->packages.Get("sinks");
  if (stage == integrator->nstages && sinks_pkg->Param<bool>("enabled")) {
    TaskRegion &sync_region_sk = tc.AddRegion(1);
    {
      for (auto &pmb : blocks) {
        auto &tl = sync_region_sk[0];
        auto &sd = pmb->meshblock_data.Get()->GetSwarmData();
        tl.AddTask(none, &SwarmContainer::ResetCommunication, sd.get());
      }
    }
    // N-body advance is a single mesh-level task (all sinks interact); it must run after the
    // per-block ResetCommunication above and before the per-block Send/Receive below.
    TaskRegion &advance_region_sk = tc.AddRegion(1);
    {
      auto &tl = advance_region_sk[0];
      auto &mu0 = pmesh->mesh_data.GetOrAdd("base", 0);
      tl.AddTask(none, Sinks::AdvanceSinksNBody, mu0.get(), integrator->dt);
    }
    TaskRegion &async_region_sk = tc.AddRegion(blocks.size());
    for (int n = 0; n < blocks.size(); n++) {
      auto &tl = async_region_sk[n];
      auto &pmb = blocks[n];
      auto &sd = pmb->meshblock_data.Get()->GetSwarmData();
      auto send = tl.AddTask(none, &SwarmContainer::Send, sd.get(), BoundaryCommSubset::all);
      tl.AddTask(send, &SwarmContainer::Receive, sd.get(), BoundaryCommSubset::all);
    }
  }

  if (stage == integrator->nstages && pmesh->adaptive) {
    TaskRegion &async_region_4 = tc.AddRegion(num_task_lists_executed_independently);
    for (int i = 0; i < blocks.size(); i++) {
      auto &tl = async_region_4[i];
      auto &u0 = blocks[i]->meshblock_data.Get("base");
      auto tag_refine =
          tl.AddTask(none, parthenon::Refinement::Tag<MeshBlockData<Real>>, u0.get());
    }
  }

  return tc;
}
} // namespace Hydro
