//========================================================================================
// AthenaPK - Self-gravity package
// Ported from Artemis (LANL, BSD-licensed).
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include <limits>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include <bvals/boundary_conditions_generic.hpp>
#include <coordinates/coordinates.hpp>
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <solvers/bicgstab_solver.hpp>
#include <solvers/mg_solver.hpp>
#include <solvers/solver_utils.hpp>

#include "../diagnostics/grav_diag.hpp"
#include "../main.hpp"          // for IDN
#include "self_gravity.hpp"
#include "poisson_equation.hpp"
#include "multipole.hpp"
#include <solvers/internal_prolongation.hpp>

namespace SelfGravity {

using namespace parthenon::BoundaryFunction;
using namespace parthenon::package::prelude;

// Selector for BC enrollment — matches any variable in the "grav" namespace.
struct any_grav : public parthenon::variable_names::base_t<true> {
  template <class... Ts>
  KOKKOS_INLINE_FUNCTION any_grav(Ts &&...args)
      : base_t<true>(std::forward<Ts>(args)...) {}
  static std::string name() { return "grav[.].*"; }
};

// Zero Dirichlet (phi = 0 on face): FixedFace BC enforces phi(face) = 0 via
// linear extrapolation into ghosts. Constant = 0 here.
template <parthenon::CoordinateDirection DIR, BCSide SIDE>
auto DirZ() {
  return [](std::shared_ptr<MeshBlockData<Real>> &rc, bool coarse) -> void {
    GenericBC<DIR, SIDE, BCType::FixedFace, any_grav>(rc, coarse, 0.0);
  };
}

// Neumann (dphi/dn = 0): just copy from interior, which is Outflow semantics.
template <parthenon::CoordinateDirection DIR, BCSide SIDE>
auto NeuZ() {
  return [](std::shared_ptr<MeshBlockData<Real>> &rc, bool coarse) -> void {
    GenericBC<DIR, SIDE, BCType::Outflow, any_grav>(rc, coarse, 0.0);
  };
}

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {
  auto pkg = std::make_shared<StateDescriptor>("self_gravity");
  const std::string block_name = "self_gravity";

  // --- Coordinate system check -----------------------------------------------
  // This port supports Cartesian only. AthenaPK currently only supports
  // UniformCartesian anyway, so we just assert and move on.
  // (Parthenon's Coordinates_t is a typedef; check it here at runtime.)
  PARTHENON_REQUIRE(
      typeid(parthenon::Coordinates_t) == typeid(parthenon::UniformCartesian),
      "Self-gravity currently only supports UniformCartesian coordinates.");

  // --- Multigrid required ----------------------------------------------------
  const bool mg_enabled =
      pin->GetOrAddBoolean("parthenon/mesh", "multigrid", false);
  PARTHENON_REQUIRE(
      mg_enabled,
      "Self-gravity requires parthenon/mesh/multigrid = true. Set it in your input file.");

  // --- 4 pi G ----------------------------------------------------------------
  // v1: require units_override. Integrating with AthenaPK's Units class is a
  // follow-on. For G in code units, user sets four_pi_G explicitly.
  const bool units_override =
      pin->GetOrAddBoolean(block_name, "units_override", true);
  PARTHENON_REQUIRE(units_override,
                    "self_gravity/units_override must be true in v1. "
                    "Compute 4*pi*G in your code units and set four_pi_G directly.");
  const Real four_pi_G = pin->GetOrAddReal(block_name, "four_pi_G", 1.0);
  pkg->AddParam("four_pi_G", four_pi_G);

  // --- Jeans swindle ---------------------------------------------------------
  auto is_periodic = [&](const std::string &face) {
    return pin->GetOrAddString("parthenon/mesh", face, "outflow") == "periodic";
  };
  const bool fully_periodic =
      is_periodic("ix1_bc") && is_periodic("ox1_bc") &&
      is_periodic("ix2_bc") && is_periodic("ox2_bc") &&
      is_periodic("ix3_bc") && is_periodic("ox3_bc");
  const bool use_swindle =
      pin->GetOrAddBoolean(block_name, "use_swindle", fully_periodic);
  PARTHENON_REQUIRE(!fully_periodic || use_swindle,
                    "Fully periodic mesh BCs require Jeans swindle (use_swindle = true).");
  if (use_swindle && !fully_periodic) {
    PARTHENON_WARN("Jeans swindle enabled on non-fully-periodic mesh. Proceed carefully.");
  }
  pkg->AddParam("use_swindle", use_swindle);

  // --- phi boundary condition enrollment -------------------------------------
  auto valid_bc = [](const std::string &s) {
    return s == "default" || s == "zero" || s == "neumann" || s == "multipole";
  };
  const std::string b_ix1 = pin->GetOrAddString(block_name, "ix1_bc", "default");
  const std::string b_ox1 = pin->GetOrAddString(block_name, "ox1_bc", "default");
  const std::string b_ix2 = pin->GetOrAddString(block_name, "ix2_bc", "default");
  const std::string b_ox2 = pin->GetOrAddString(block_name, "ox2_bc", "default");
  const std::string b_ix3 = pin->GetOrAddString(block_name, "ix3_bc", "default");
  const std::string b_ox3 = pin->GetOrAddString(block_name, "ox3_bc", "default");
  PARTHENON_REQUIRE(valid_bc(b_ix1), "Invalid self_gravity/ix1_bc: " + b_ix1);
  PARTHENON_REQUIRE(valid_bc(b_ox1), "Invalid self_gravity/ox1_bc: " + b_ox1);
  PARTHENON_REQUIRE(valid_bc(b_ix2), "Invalid self_gravity/ix2_bc: " + b_ix2);
  PARTHENON_REQUIRE(valid_bc(b_ox2), "Invalid self_gravity/ox2_bc: " + b_ox2);
  PARTHENON_REQUIRE(valid_bc(b_ix3), "Invalid self_gravity/ix3_bc: " + b_ix3);
  PARTHENON_REQUIRE(valid_bc(b_ox3), "Invalid self_gravity/ox3_bc: " + b_ox3);

  using BF = parthenon::BoundaryFace;
  constexpr auto LL = BCSide::Inner;
  constexpr auto RR = BCSide::Outer;
  if (b_ix1 == "zero") pkg->UserBoundaryFunctions[BF::inner_x1].push_back(DirZ<parthenon::X1DIR, LL>());
  if (b_ox1 == "zero") pkg->UserBoundaryFunctions[BF::outer_x1].push_back(DirZ<parthenon::X1DIR, RR>());
  if (b_ix2 == "zero") pkg->UserBoundaryFunctions[BF::inner_x2].push_back(DirZ<parthenon::X2DIR, LL>());
  if (b_ox2 == "zero") pkg->UserBoundaryFunctions[BF::outer_x2].push_back(DirZ<parthenon::X2DIR, RR>());
  if (b_ix3 == "zero") pkg->UserBoundaryFunctions[BF::inner_x3].push_back(DirZ<parthenon::X3DIR, LL>());
  if (b_ox3 == "zero") pkg->UserBoundaryFunctions[BF::outer_x3].push_back(DirZ<parthenon::X3DIR, RR>());
  if (b_ix1 == "neumann") pkg->UserBoundaryFunctions[BF::inner_x1].push_back(NeuZ<parthenon::X1DIR, LL>());
  if (b_ox1 == "neumann") pkg->UserBoundaryFunctions[BF::outer_x1].push_back(NeuZ<parthenon::X1DIR, RR>());
  if (b_ix2 == "neumann") pkg->UserBoundaryFunctions[BF::inner_x2].push_back(NeuZ<parthenon::X2DIR, LL>());
  if (b_ox2 == "neumann") pkg->UserBoundaryFunctions[BF::outer_x2].push_back(NeuZ<parthenon::X2DIR, RR>());
  if (b_ix3 == "neumann") pkg->UserBoundaryFunctions[BF::inner_x3].push_back(NeuZ<parthenon::X3DIR, LL>());
  if (b_ox3 == "neumann") pkg->UserBoundaryFunctions[BF::outer_x3].push_back(NeuZ<parthenon::X3DIR, RR>());
  // Multipole (WS-5a): exterior Cartesian expansion of interior mass through the quadrupole.
  // Enrolled as the per-block post-solve BC; the packed solver-internal path is in
  // PoissonEquation::SetBoundary (both keyed on face type 3).
  if (b_ix1 == "multipole") pkg->UserBoundaryFunctions[BF::inner_x1].push_back(MultipoleBC<parthenon::X1DIR, LL>);
  if (b_ox1 == "multipole") pkg->UserBoundaryFunctions[BF::outer_x1].push_back(MultipoleBC<parthenon::X1DIR, RR>);
  if (b_ix2 == "multipole") pkg->UserBoundaryFunctions[BF::inner_x2].push_back(MultipoleBC<parthenon::X2DIR, LL>);
  if (b_ox2 == "multipole") pkg->UserBoundaryFunctions[BF::outer_x2].push_back(MultipoleBC<parthenon::X2DIR, RR>);
  if (b_ix3 == "multipole") pkg->UserBoundaryFunctions[BF::inner_x3].push_back(MultipoleBC<parthenon::X3DIR, LL>);
  if (b_ox3 == "multipole") pkg->UserBoundaryFunctions[BF::outer_x3].push_back(MultipoleBC<parthenon::X3DIR, RR>);

  // Cache the per-face BC type (0=default/other, 1=zero/Dirichlet, 2=neumann/outflow,
  // 3=multipole) so the packed PoissonEquation::SetBoundary can apply them without the
  // per-block BC dispatch. Order matches parthenon::BoundaryFace:
  // inner_x1,outer_x1,inner_x2,outer_x2,inner_x3,outer_x3.
  auto bc_code = [](const std::string &s) {
    return s == "zero" ? 1 : (s == "neumann" ? 2 : (s == "multipole" ? 3 : 0));
  };
  std::array<int, 6> grav_bc_face_type = {bc_code(b_ix1), bc_code(b_ox1), bc_code(b_ix2),
                                          bc_code(b_ox2), bc_code(b_ix3), bc_code(b_ox3)};
  pkg->AddParam("grav_bc_face_type", grav_bc_face_type);
  // Runtime escape hatch / A-B switch: self_gravity/packed_bc=false forces the original
  // per-block ApplyBoundaryConditionsOnCoarseOrFineMD path inside SetBoundary.
  const bool packed_bc = pin->GetOrAddBoolean(block_name, "packed_bc", true);
  pkg->AddParam("grav_packed_bc", packed_bc);

  // Multipole moments (WS-5a): mutable POD recomputed each step in FillPoissonRHS and
  // captured by value into the (device) boundary kernels. Registered unconditionally so
  // both BC paths can Param<>() it; four_pi_G is fixed, the rest start at zero.
  const bool has_multipole = (b_ix1 == "multipole") || (b_ox1 == "multipole") ||
                             (b_ix2 == "multipole") || (b_ox2 == "multipole") ||
                             (b_ix3 == "multipole") || (b_ox3 == "multipole");
  // Multipole BCs rely on the packed SetBoundary applying HOMOGENEOUS Dirichlet inside the
  // solve (with the inhomogeneous data lifted into the RHS). packed_bc=false would route
  // the operator through the enrolled inhomogeneous per-block MultipoleBC, making the
  // operator affine and diverging the preconditioned solve -- so forbid that combination.
  PARTHENON_REQUIRE(!has_multipole || packed_bc,
                    "self_gravity: multipole gravity BCs require packed_bc=true.");
  // NEW-1 (audit 2026-07-19): the multipole moments are assembled from the FULL density
  // while the swindle subtracts rho_mean from the RHS -- combining them would solve an
  // equation whose boundary data and source are mutually inconsistent. The multipole BC
  // is for isolated (non-periodic) boxes where the swindle is off anyway.
  PARTHENON_REQUIRE(!has_multipole || !use_swindle,
                    "self_gravity: multipole gravity BCs are incompatible with "
                    "use_swindle=true (moments use full rho, swindle RHS uses rho - "
                    "rho_mean). Disable the swindle for isolated-box multipole BCs.");
  MultipoleMoments mm0{};
  mm0.four_pi_G = four_pi_G;
  pkg->AddParam("grav_multipole_moments", mm0, /*is_mutable=*/true);
  pkg->AddParam("grav_has_multipole", has_multipole);

  // --- Fields ----------------------------------------------------------------
  // phi: solution variable. Needs ghost fill, fluxes, GMG prolong/restrict for AMR MG.
  {
    std::vector<parthenon::MetadataFlag> flags{
        Metadata::Cell,       Metadata::Independent,   Metadata::FillGhost,
        Metadata::WithFluxes, Metadata::GMGRestrict,   Metadata::GMGProlongate};
    Metadata m(flags);
    // phi MG prolongation/restriction ops. Artemis uses
    // ArtemisUtils::EnrollArtemisRefinementOps; we use Parthenon's built-in
    // ProlongatePiecewiseConstant / RestrictAverage, the safe elliptic defaults.
    m.RegisterRefinementOps<parthenon::refinement_ops::ProlongatePiecewiseConstant,
                            parthenon::refinement_ops::RestrictAverage>();
    pkg->AddField<grav::phi>(m);
  }
  // rhs: source. Derived, OneCopy. Matches upstream Parthenon poisson_gmg example.
  // rhs carries FillGhost + refinement ops to match Artemis: the GMG V-cycle
  // restricts/exchanges rhs across ranks, so it needs ghost data and AMR ops.
  {
    Metadata m({Metadata::Cell, Metadata::Derived, Metadata::OneCopy,
                Metadata::FillGhost});
    m.RegisterRefinementOps<parthenon::refinement_ops::ProlongatePiecewiseConstant,
                            parthenon::refinement_ops::RestrictAverage>();
    pkg->AddField<grav::rhs>(m);
  }

  // --- FillDerived: compute rhs from density on every timestep ---------------
  pkg->FillDerivedMesh = FillPoissonRHS;

  // --- Solver construction ---------------------------------------------------
  using PoissEq = PoissonEquation<grav::phi>;
  using prolongator_t = parthenon::solvers::ProlongationBlockInteriorZeroDirichlet;
  using preconditioner_t = parthenon::solvers::MGSolver<PoissEq, prolongator_t>;
  const std::string solver_params_block = block_name + "/solver_params";
  // Runtime-selectable solver. Default "BiCGSTAB" reproduces production behaviour
  // bit-for-bit. "MG" / "Multigrid" uses the *pure* geometric multigrid solver,
  // which has no global inner-products (BiCGSTAB needs ~2 all-reduces/iteration),
  // attacking the latency-bound bottleneck of the GPU self-gravity solve.
  // NOTE: solver=MG (pure multigrid) needs an adequate smoother on AMR. SRJ1 (one
  // weighted-Jacobi sweep) is enough on a UNIFORM grid but NOT across AMR fine-coarse
  // boundaries, where the coarse-grid correction injects high-frequency error each
  // V-cycle (the standalone V-cycle then has spectral radius > 1 and diverges). Use
  // SRJ2/SRJ3 for MG on AMR; SRJ2 is the parthenon MGParams default and converges in
  // ~6 V-cycles. BiCGSTAB tolerates SRJ1 because its Krylov outer loop stabilises a
  // non-contractive preconditioner.
  const std::string solver_type =
      pin->GetOrAddString(solver_params_block, "solver", "BiCGSTAB");
  std::shared_ptr<parthenon::solvers::SolverBase> psolver;
  if (solver_type == "MG" || solver_type == "Multigrid") {
    psolver = std::make_shared<parthenon::solvers::MGSolver<PoissEq, prolongator_t>>(
        /*container_base=*/"base",
        /*container_u=*/"phi",
        /*container_rhs=*/"rhs",
        pin, solver_params_block, PoissEq(pin, block_name));
  } else {
    psolver =
        std::make_shared<parthenon::solvers::BiCGSTABSolver<PoissEq, preconditioner_t>>(
            /*container_base=*/"base",
            /*container_u=*/"phi",
            /*container_rhs=*/"rhs",
            pin, solver_params_block, PoissEq(pin, block_name));
  }
  pkg->AddParam("solver_pointer", psolver);

  // --- VALIDATION WP-5: Poisson-solver convergence reporting -----------------
  // The BiCGSTAB completion task returns TaskStatus::complete both when a tolerance is
  // met and when max_iterations is exhausted, and reports neither -- so a solve that
  // silently bails at the ceiling is indistinguishable from a converged one. Record the
  // iteration count / final residual and flag the ceiling case. See diagnostics/grav_diag.hpp.
  // Default OFF; when off nothing is registered and no task is added => bit-identical.
  const bool solver_diag = pin->GetOrAddBoolean(block_name, "solver_diag", false);
  pkg->AddParam<>("solver_diag", solver_diag);

  // B2 (2026-08-02): the NON-CONVERGENCE WARNING is unconditional, unlike the hst columns.
  // Rationale: a Poisson solve that quietly gave up corrupts the potential for that step, and
  // production runs do NOT set solver_diag (it adds hst columns, which would shift every column
  // index and break existing analysis scripts). Gating the warning on the diagnostic meant the
  // one failure mode that most needs to be loud was silent in exactly the runs that matter.
  // These two Params are always registered; the warn task reads solver state only and touches
  // no field data, so the OFF-state remains bit-identical.
  const int max_iters_always =
      pin->GetOrAddInteger(solver_params_block, "max_iterations", 1000);
  pkg->AddParam<int>("grav_max_iters", max_iters_always);
  pkg->AddParam<int>("grav_nonconv_count", 0, true);

  // B3 (2026-08-02): the convergence criterion is an ABSOLUTE residual unless the deck asks
  // otherwise. From bicgstab_solver.hpp:54-66, `relative_residual` defaults to FALSE, and in
  // that branch  absolute_residual_tolerance = residual_tolerance  while
  // relative_residual_tolerance = 0. The test is therefore on
  //     rms_res = sqrt( sum(r^2) / total_cells )
  // of  grad^2 phi = 4 pi G rho,  with a FIXED ceiling. Through a collapse the RHS grows by
  // many decades while that ceiling does not move, so the *relative* accuracy being demanded
  // tightens monotonically -- the run gets progressively harder to converge exactly when the
  // physics gets most interesting. Measured at smoke scale (WP-5): grav-res climbed
  // 3.9e-8 -> 9.5e-7 over ~1 t0, i.e. to 95% of a 1e-6 tolerance, before any deep collapse.
  // This is NOT auto-corrected here: switching to a relative criterion changes which solves
  // are accepted and is therefore result-changing. It is made LOUD instead.
  const bool rel_res = pin->GetOrAddBoolean(solver_params_block, "relative_residual", false);
  if (!rel_res && parthenon::Globals::my_rank == 0) {
    const Real res_tol = pin->GetOrAddReal(solver_params_block, "residual_tolerance", 1.0e-12);
    std::cout << "## NOTE [self-gravity] convergence uses an ABSOLUTE residual tolerance ("
              << res_tol << "); <" << solver_params_block
              << "> relative_residual is false. As the collapse deepens the RHS grows while "
                 "this ceiling does not, so the effective relative tolerance tightens "
                 "monotonically. Watch grav-res (self_gravity/solver_diag=true) and the "
                 "non-convergence warnings. See B3 in VALIDATION_PLAN.md."
              << std::endl;
  }

  if (solver_diag) {
    // Mirror the ceiling the solver itself parsed, so the "did it bail" test uses the
    // same number rather than a hard-coded default. The key was already consumed above
    // by BiCGSTABParams/MGParams, so GetOrAdd here just reads it back.
    pkg->AddParam<Real>("grav_last_iters", -1.0, true);
    pkg->AddParam<Real>("grav_last_res", -1.0, true);
    pkg->AddParam<Real>("grav_last_nonconv", 0.0, true);

    using Diagnostics::GravDiag;
    using Diagnostics::GravDiagReport;
    parthenon::HstVar_list hst_grav = {};
    auto add_grav = [&hst_grav](GravDiag w, const std::string &label) {
      // max, not sum: the underlying scalars come from a global all-reduce inside the
      // solver and are identical on every partition, so sum would multiply by the
      // partition count.
      hst_grav.emplace_back(parthenon::HistoryOutputVar(
          parthenon::UserHistoryOperation::max,
          [w](MeshData<Real> *md) { return GravDiagReport(md, w); }, label));
    };
    add_grav(GravDiag::iters, "grav-iters");
    add_grav(GravDiag::res, "grav-res");
    add_grav(GravDiag::nonconv, "grav-nonconv");
    pkg->AddParam<>(parthenon::hist_param_key, hst_grav);

    if (parthenon::Globals::my_rank == 0) {
      std::cout << "## Self-gravity solver diagnostics ON (self_gravity/solver_diag): "
                   "hst grav-iters, grav-res, grav-nonconv; max_iterations="
                << max_iters_always << "." << std::endl;
    }
  }

  return pkg;
}

// FillPoissonRHS: assemble rhs = 4 pi G (rho - rho_mean) from gas density.
// Runs over ENTIRE (including ghosts) so solver sees consistent ghost values.
void FillPoissonRHS(MeshData<Real> *md) {
  auto pm = md->GetParentPointer();
  auto &grav_pkg = pm->packages.Get("self_gravity");
  const bool use_swindle = grav_pkg->Param<bool>("use_swindle");
  const Real four_pi_G = grav_pkg->Param<Real>("four_pi_G");
  const bool has_multipole = grav_pkg->Param<bool>("grav_has_multipole");

  // WP-13 (2026-08-02): read the density from "cons", NOT "prim".
  //
  // This function is self_gravity's FillDerivedMesh. `prim` is produced by *Hydro's*
  // FillDerivedMesh (ConsToPrim), and Parthenon does not guarantee that one package's
  // FillDerived runs before another's. On a FRESH start that was harmless: initialization
  // performs several FillDerived passes and the last one sees a populated `prim`. On a
  // RESTART there is a single init pass, and self_gravity's ran first -- so `prim` was still
  // ALL ZEROS (measured: "prim rho entire[0,0] interior[0,0]"), the Poisson RHS came out
  // identically zero, and the first solve of the first step fed BiCGSTAB a zero right-hand
  // side. The Krylov recurrence then divided by a zero residual norm and returned phi = NaN
  // in all 512000 cells, which the gravity kick applied to momentum and energy, flooring the
  // entire domain in one step (mass 5.17e4 -> 0.1406, KE = nan). Silent, and fatal.
  //
  // cons(IDN) and prim(IDN) are the SAME number -- density is not transformed by
  // ConsToPrim -- so this is numerically a no-op on every ordinary step, while removing the
  // inter-package FillDerived ordering dependency entirely. `cons` is Independent and is
  // restored directly from the restart file, so it is valid before any derived quantity is.
  // The density floor ConsToPrim applies (adiabatic_glmmhd.hpp:152 /
  // adiabatic_hydro.hpp:81 do `u_d = max(u_d, density_floor_)`). Applying it here is what
  // makes reading `cons` EXACTLY equivalent to reading `prim`: without it the RHS uses
  // unfloored (possibly sub-floor or negative) densities and the trajectory changes -- which
  // it measurably did, diverging from the pre-fix fresh run by cycle 2. A no-op fix must
  // reproduce prim bit-for-bit, floor included.
  const Real rho_floor = pm->packages.Get("Hydro")->Param<Real>("grav_rho_floor");

  const auto &cons_pack = md->PackVariables(std::vector<std::string>{"cons"});
  auto &resolved = pm->resolved_packages;
  auto desc_rhs = parthenon::MakePackDescriptor<grav::rhs>(resolved.get());
  auto rhs_pack = desc_rhs.GetPack(md);

  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  IndexRange ibe = md->GetBoundsI(IndexDomain::entire);
  IndexRange jbe = md->GetBoundsJ(IndexDomain::entire);
  IndexRange kbe = md->GetBoundsK(IndexDomain::entire);
  const int nblocks = md->NumBlocks();

  // --- Mean density (for Jeans swindle) via par_reduce + MPI Allreduce -------
  Real grav_mean_rho = 0.0;
  if (use_swindle) {
    Real total_mass = 0.0, total_volume = 0.0;
    parthenon::par_reduce(
        parthenon::loop_pattern_mdrange_tag, "SG::TotalMass",
        parthenon::DevExecSpace(), 0, nblocks - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i,
                      Real &lmass, Real &lvol) {
          const auto &coords = cons_pack.GetCoords(b);
          const Real vv = coords.CellVolume(k, j, i);
          lvol += vv;
          lmass += (cons_pack(b, IDN, k, j, i) > rho_floor ? cons_pack(b, IDN, k, j, i) : rho_floor) * vv;
        },
        Kokkos::Sum<Real>(total_mass), Kokkos::Sum<Real>(total_volume));
    Kokkos::fence();
#ifdef MPI_PARALLEL
    Real buf[2] = {total_mass, total_volume};
    PARTHENON_MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, buf, 2, MPI_PARTHENON_REAL,
                                      MPI_SUM, MPI_COMM_WORLD));
    total_mass = buf[0];
    total_volume = buf[1];
#endif
    grav_mean_rho = total_mass / total_volume;
  }

  // --- Multipole moments (for multipole exterior BC) -------------------------
  // Ten raw moments about the coordinate origin (mass, mass-weighted x_i, and
  // mass-weighted x_i x_j) summed over interior cells, then a global MPI_Allreduce.
  // From these derive the center of mass and the traceless quadrupole about it.
  // Uses the FULL density (not rho - rho_mean): the multipole BC is only used on
  // isolated (non-periodic) boxes where the swindle is off, so RHS = 4piG*rho and
  // the exterior expansion is of that same rho -- consistent by construction.
  if (has_multipole) {
    Kokkos::View<Real *> mom("SG::moments", 10);
    Kokkos::deep_copy(mom, 0.0);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "SG::Moments", parthenon::DevExecSpace(), 0, nblocks - 1,
        kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &coords = cons_pack.GetCoords(b);
          const Real vv = coords.CellVolume(k, j, i);
          const Real dm = (cons_pack(b, IDN, k, j, i) > rho_floor ? cons_pack(b, IDN, k, j, i) : rho_floor) * vv;
          const Real x = coords.Xc<parthenon::X1DIR>(k, j, i);
          const Real y = coords.Xc<parthenon::X2DIR>(k, j, i);
          const Real z = coords.Xc<parthenon::X3DIR>(k, j, i);
          Kokkos::atomic_add(&mom(0), dm);
          Kokkos::atomic_add(&mom(1), dm * x);
          Kokkos::atomic_add(&mom(2), dm * y);
          Kokkos::atomic_add(&mom(3), dm * z);
          Kokkos::atomic_add(&mom(4), dm * x * x);
          Kokkos::atomic_add(&mom(5), dm * y * y);
          Kokkos::atomic_add(&mom(6), dm * z * z);
          Kokkos::atomic_add(&mom(7), dm * x * y);
          Kokkos::atomic_add(&mom(8), dm * x * z);
          Kokkos::atomic_add(&mom(9), dm * y * z);
        });
    Kokkos::fence();
    auto mom_h = Kokkos::create_mirror_view(mom);
    Kokkos::deep_copy(mom_h, mom);
    Real h[10];
    for (int t = 0; t < 10; ++t) h[t] = mom_h(t);
#ifdef MPI_PARALLEL
    PARTHENON_MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, h, 10, MPI_PARTHENON_REAL, MPI_SUM,
                                      MPI_COMM_WORLD));
#endif
    MultipoleMoments mm;
    mm.four_pi_G = four_pi_G;
    mm.M = h[0];
    const Real invM = (h[0] != 0.0) ? 1.0 / h[0] : 0.0;
    mm.cx = h[1] * invM;
    mm.cy = h[2] * invM;
    mm.cz = h[3] * invM;
    // Second moments of mass about the COM: mu_ij = <x_i x_j> - M c_i c_j.
    const Real muxx = h[4] - mm.M * mm.cx * mm.cx;
    const Real muyy = h[5] - mm.M * mm.cy * mm.cy;
    const Real muzz = h[6] - mm.M * mm.cz * mm.cz;
    const Real muxy = h[7] - mm.M * mm.cx * mm.cy;
    const Real muxz = h[8] - mm.M * mm.cx * mm.cz;
    const Real muyz = h[9] - mm.M * mm.cy * mm.cz;
    // Traceless Cartesian quadrupole Q_ij = 3 mu_ij - delta_ij tr(mu).
    const Real tr = muxx + muyy + muzz;
    mm.Qxx = 3.0 * muxx - tr;
    mm.Qyy = 3.0 * muyy - tr;
    mm.Qzz = 3.0 * muzz - tr;
    mm.Qxy = 3.0 * muxy;
    mm.Qxz = 3.0 * muxz;
    mm.Qyz = 3.0 * muyz;
    grav_pkg->UpdateParam("grav_multipole_moments", mm);
  }

  // --- Fill rhs over entire domain (incl. ghosts) ----------------------------
  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "SG::SetRHS", parthenon::DevExecSpace(), 0, nblocks - 1,
      kbe.s, kbe.e, jbe.s, jbe.e, ibe.s, ibe.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const Real rho = (cons_pack(b, IDN, k, j, i) > rho_floor ? cons_pack(b, IDN, k, j, i) : rho_floor);
        rhs_pack(b, te, grav::rhs(), k, j, i) = four_pi_G * (rho - grav_mean_rho);
      });


  // --- Multipole boundary lift ----------------------------------------------
  // The Poisson solver runs a LINEAR operator with HOMOGENEOUS Dirichlet ghosts on the
  // multipole faces (see PoissonEquation::SetBoundary). The inhomogeneous exterior-
  // expansion face value Phi_mp is folded into the RHS here: for an interior cell adjacent
  // to a boundary face, the true ghost 2*Phi_mp - phi_int vs the homogeneous ghost -phi_int
  // differ by 2*Phi_mp, which the discrete Laplacian adds as +2*Phi_mp/dx_n^2. Moving it to
  // the RHS:  b -= 2*Phi_mp(face) / dx_n^2  at each boundary-adjacent interior cell (edges/
  // corners get one term per touching face). Post-solve, the enrolled per-block MultipoleBC
  // restores the physical ghost values for ApplyGravitySource.
  if (has_multipole) {
    const MultipoleMoments mm =
        grav_pkg->Param<MultipoleMoments>("grav_multipole_moments");
    const auto &bctype = grav_pkg->Param<std::array<int, 6>>("grav_bc_face_type");
    for (int b = 0; b < nblocks; ++b) {
      auto *pmb = md->GetBlockData(b)->GetBlockPointer();
      for (int f = 0; f < 6; ++f) {
        if (bctype[f] != 3) continue;
        if (pmb->boundary_flag[f] != parthenon::BoundaryFlag::user) continue;
        const int dd = f / 2;         // 0->x1, 1->x2, 2->x3
        const bool inner = (f % 2 == 0);
        const int foff = inner ? 0 : 1; // Xf face-index offset for this side
        // Interior cell slab adjacent to this face (one cell thick along dd).
        const int i0 = (dd == 0) ? (inner ? ib.s : ib.e) : ib.s;
        const int i1 = (dd == 0) ? (inner ? ib.s : ib.e) : ib.e;
        const int j0 = (dd == 1) ? (inner ? jb.s : jb.e) : jb.s;
        const int j1 = (dd == 1) ? (inner ? jb.s : jb.e) : jb.e;
        const int k0 = (dd == 2) ? (inner ? kb.s : kb.e) : kb.s;
        const int k1 = (dd == 2) ? (inner ? kb.s : kb.e) : kb.e;
        parthenon::par_for(
            DEFAULT_LOOP_PATTERN, "SG::MultipoleRHSLift", parthenon::DevExecSpace(),
            b, b, k0, k1, j0, j1, i0, i1,
            KOKKOS_LAMBDA(const int bb, const int k, const int j, const int i) {
              const auto &coords = cons_pack.GetCoords(bb);
              Real px = coords.Xc<parthenon::X1DIR>(k, j, i);
              Real py = coords.Xc<parthenon::X2DIR>(k, j, i);
              Real pz = coords.Xc<parthenon::X3DIR>(k, j, i);
              Real dxn;
              if (dd == 0) {
                px = coords.Xf<parthenon::X1DIR>(k, j, i + foff);
                dxn = coords.Dxc<parthenon::X1DIR>(k, j, i);
              } else if (dd == 1) {
                py = coords.Xf<parthenon::X2DIR>(k, j + foff, i);
                dxn = coords.Dxc<parthenon::X2DIR>(k, j, i);
              } else {
                pz = coords.Xf<parthenon::X3DIR>(k + foff, j, i);
                dxn = coords.Dxc<parthenon::X3DIR>(k, j, i);
              }
              const Real val = MultipolePhi(mm, px, py, pz);
              rhs_pack(bb, te, grav::rhs(), k, j, i) -= 2.0 * val / (dxn * dxn);
            });
      }
    }
  }
}

// ApplyGravitySource: momentum += rho * g * dt, energy += flux-weighted work.
// Flux-weighted energy (Artemis style) uses the hydro mass flux across each face
// to compute gravitational work. Requires hydro fluxes still be in memory.
TaskStatus ApplyGravitySource(MeshData<Real> *md, const parthenon::SimTime &tm,
                              const Real beta_dt) {
  auto pm = md->GetParentPointer();

  // Pack hydro cons + prim + phi, and cons-with-fluxes for mass flux access.
  const auto &cons_pack = md->PackVariables(std::vector<std::string>{"cons"});
  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});
  auto &resolved = pm->resolved_packages;
  auto desc_phi = parthenon::MakePackDescriptor<grav::phi>(resolved.get());
  auto phi_pack = desc_phi.GetPack(md);

  // Mass flux: pack the "cons" variable AND its fluxes BY NAME (density flux =
  // component IDN). Packing by name rather than the {Independent} metadata flag
  // is required here: grav::phi is also Independent + WithFluxes, so selecting by
  // flag would pull phi into the pack and make the flat IDN index point at phi's
  // flux instead of the gas mass flux whenever phi sorts ahead of cons.
  auto cons_flx_pack = md->PackVariablesAndFluxes(std::vector<std::string>{"cons"},
                                                  std::vector<std::string>{"cons"});

  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);

  const int nblocks = md->NumBlocks();
  const int ndim = pm->ndim;
  const bool multi_d = (ndim > 1);
  const bool three_d = (ndim > 2);

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "SG::ApplyGravity", parthenon::DevExecSpace(),
      0, nblocks - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        auto &cons = cons_pack(b);
        auto &cons_flx = cons_flx_pack(b);
        const auto &prim = prim_pack(b);
        const auto &coords = cons_pack.GetCoords(b);

        const Real dx1 = coords.Dxc<1>(k, j, i);
        const Real dx2 = multi_d ? coords.Dxc<2>(k, j, i) : 1.0;
        const Real dx3 = three_d ? coords.Dxc<3>(k, j, i) : 1.0;

        const Real hdtodx1 = 0.5 * beta_dt / dx1;
        const Real hdtodx2 = multi_d ? 0.5 * beta_dt / dx2 : 0.0;
        const Real hdtodx3 = three_d ? 0.5 * beta_dt / dx3 : 0.0;

        // phi differences (Artemis convention: dpl = -(phi_c - phi_{c-1}))
        const Real phic = phi_pack(b, te, grav::phi(), k, j, i);
        const Real dpl1 = -(phic - phi_pack(b, te, grav::phi(), k, j, i - 1));
        const Real dpr1 = -(phi_pack(b, te, grav::phi(), k, j, i + 1) - phic);
        const Real dpl2 = multi_d
            ? -(phic - phi_pack(b, te, grav::phi(), k, j - 1, i))
            : 0.0;
        const Real dpr2 = multi_d
            ? -(phi_pack(b, te, grav::phi(), k, j + 1, i) - phic)
            : 0.0;
        const Real dpl3 = three_d
            ? -(phic - phi_pack(b, te, grav::phi(), k - 1, j, i))
            : 0.0;
        const Real dpr3 = three_d
            ? -(phi_pack(b, te, grav::phi(), k + 1, j, i) - phic)
            : 0.0;

        // Momentum update: rho * g * dt, centered difference of phi.
        const Real rho = prim(IDN, k, j, i);
        cons(IM1, k, j, i) += rho * hdtodx1 * (dpl1 + dpr1);
        if (multi_d) cons(IM2, k, j, i) += rho * hdtodx2 * (dpl2 + dpr2);
        if (three_d) cons(IM3, k, j, i) += rho * hdtodx3 * (dpl3 + dpr3);

        // Energy update: flux-weighted work.
        // mass_flux(face) = cons.flux(IVn, IDN, k, j, i).
        // Work done by gravity = sum over faces of (mass_flux * dphi/face) * (0.5 dt / dx).
        Real de = hdtodx1 * (cons_flx.flux(X1DIR, IDN, k, j, i)     * dpl1 +
                             cons_flx.flux(X1DIR, IDN, k, j, i + 1) * dpr1);
        if (multi_d) {
          de += hdtodx2 * (cons_flx.flux(X2DIR, IDN, k, j,     i) * dpl2 +
                           cons_flx.flux(X2DIR, IDN, k, j + 1, i) * dpr2);
        }
        if (three_d) {
          de += hdtodx3 * (cons_flx.flux(X3DIR, IDN, k,     j, i) * dpl3 +
                           cons_flx.flux(X3DIR, IDN, k + 1, j, i) * dpr3);
        }
        cons(IEN, k, j, i) += de;
      });

  return TaskStatus::complete;
}

} // namespace SelfGravity
