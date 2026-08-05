//========================================================================================
// AthenaPK - physics-based AMR for magnetized non-ideal collapse (flagship Phase 7).
//
// Combined refinement criterion `jeans_nonideal`: Truelove Jeans PLUS current-sheet resolution.
// The flux-retention result is set in the thin current layers / non-ideal diffusion regions --
// if a current sheet reverses the field over ~1 cell it is under-resolved and the reconnection /
// flux loss there is a numerical artifact. This adds a criterion that keeps such layers resolved.
//
// OPT-IN (type = jeans_nonideal); the default `type = jeans` path is untouched (production
// bit-identical). Refine if EITHER Jeans OR current-sheet demands it; derefine only if BOTH allow,
// so it never de-refines a Jeans-required region.
//========================================================================================

#include <cmath>
#include <limits>

#include "../eos/adiabatic_glmmhd.hpp"
#include "../hydro/hydro.hpp"
#include "../main.hpp"
#include "refinement.hpp"

namespace refinement {
namespace nonideal {

using parthenon::IndexDomain;
using parthenon::IndexRange;

// L_B = |B| / |curl B| is the field-reversal scale; refine when it is resolved by fewer than
// `curr_nsheet` cells (a thin current sheet). Derefine when comfortably resolved. Combined with
// the Jeans criterion (min Jeans-length-in-cells njmin) via the refine/derefine logic below.
parthenon::AmrTag JeansNonideal(MeshBlockData<Real> *rc) {
  auto pmb = rc->GetBlockPointer();
  auto w = rc->Get("prim").data;

  auto hydro_pkg = pmb->packages.Get("Hydro");
  const Real njeans = hydro_pkg->Param<Real>("refinement/njeans");
  const Real curr_nsheet = hydro_pkg->Param<Real>("refinement/curr_nsheet");
  const Real curr_rho_thresh = hydro_pkg->Param<Real>("refinement/curr_rho_thresh");
  const int curr_max_level = hydro_pkg->Param<int>("refinement/curr_max_level");
  const int block_level = pmb->loc.level();
  const bool mhd = (hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd);

  const Real dx = pmb->coords.Dxc<1>(0);
  const Real dy = pmb->coords.Dxc<2>(0);
  const Real dz = pmb->coords.Dxc<3>(0);

  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
  // Jeans: lambda_J/dx = fac * v/sqrt(rho), with lambda_J = 2 pi v / sqrt(rho) at 4 pi G = 1.
  // AUDIT 2026-08-05 (B2): jeans.cpp documents that it assumes cubic cells and uses dx alone;
  // this file silently inherited that while using dx, dy, dz properly a few lines below in the
  // current-sheet stencil. Use the SMALLEST edge, which is the conservative choice: on a cubic
  // mesh it is identical to the old dx (so production is bit-identical), and on a stretched
  // one it refines rather than under-resolving.
  const Real fac = 2.0 * M_PI / std::min(dx, std::min(dy, dz));

  // Jeans length in cells (min over block).
  //
  // AUDIT 2026-08-05 (A1) -- MIRRORS refinement/jeans.cpp. The sound speed comes from the
  // EOS the integrator uses, not from a hard-coded gamma-law: with <hydro> eos = hydrogen
  // the tabulated Gamma_1 is ~5/3 in the cold envelope (H2 rotation frozen below 85 K) and
  // ~1.14 through dissociation, so a gamma=1.4 formula over-refines in the first regime and
  // UNDER-refines in the second. jeans_nonideal is MHD-only (the current-sheet half needs B),
  // and eos=hydrogen is glmmhd-only, so fetching the GLMMHD EOS here is always valid.
  // Measurement and full rationale: the header note in refinement/jeans.cpp.
  const auto eos = hydro_pkg->Param<AdiabaticGLMMHDEOS>("eos");
  Real njmin = std::numeric_limits<Real>::max();
  pmb->par_reduce(
      "jeans_ni: jeans", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int k, const int j, const int i, Real &lnjmin) {
        const Real rho = w(IDN, k, j, i);
        Real prim_c[NHYDRO];
        prim_c[IDN] = rho;
        prim_c[IPR] = w(IPR, k, j, i);
        Real v = eos.SoundSpeed(prim_c);
        if (mhd) {
          const Real bsq = w(IB1, k, j, i) * w(IB1, k, j, i) +
                           w(IB2, k, j, i) * w(IB2, k, j, i) +
                           w(IB3, k, j, i) * w(IB3, k, j, i);
          v += std::sqrt(bsq / rho);
        }
        lnjmin = std::min(lnjmin, v / std::sqrt(rho));
      },
      Kokkos::Min<Real>(njmin));
  njmin *= fac;

  // Current-sheet resolution: min over block of (L_B / dx), L_B = |B|/|curl B| (central diffs on
  // the FillGhost prim B). MHD only; skipped (=huge) for hydro.
  //
  // DENSITY-GATED: only cells with rho > curr_rho_thresh contribute. Current sheets set the
  // flux-loss/reconnection ONLY where there is dense gas to lose flux from; the diffuse turbulent
  // ENVELOPE is full of numerical current sheets (field reversals everywhere) that are irrelevant
  // to the fossil-flux result. Without this gate curr_nsheet refines the whole turbulent volume
  // (block runaway 526->1786 pre-collapse -> OOM, validated 2026-07-27). The gate restricts
  // current-sheet refinement to the collapsing region (near/above the initial core density).
  // LEVEL CAP: current-sheet refinement is an enhancement, not the essential collapse criterion.
  // Jeans (Truelove) drives the deep core to a converged collapse; letting the current-sheet
  // criterion ALSO drive the deepest levels explodes the block count in the dense collapsing core
  // (validated: jeans_nonideal makes 1.5x [128^3] to 3x [256^3] more blocks than plain jeans ->
  // deep-core regrid OOM). Capping current-sheet refinement at curr_max_level keeps it in the
  // moderate-depth flux-loss region while Jeans owns the deepest core, bounding memory.
  // AUDIT 2026-08-05 (B1) -- THE CAP NOW GATES THE ACTION, NOT THE MEASUREMENT.
  // Previously the reduce itself was skipped at block_level >= curr_max_level, leaving
  // lbmin = numeric_limits::max(). That value reads as "the field reversal is infinitely
  // well resolved", so curr_deref came out unconditionally TRUE at the cap. Consequence:
  // a parent one level BELOW the cap sees a thin sheet and refines; its children sit AT the
  // cap, never measure the sheet, and therefore vote to derefine; the parent re-refines.
  // A refine/derefine limit cycle, damped by derefine_count=50 but not removed.
  //
  // The fix is to keep measuring at every level and apply the cap only to the refine vote
  // (below). Children at the cap then see the sheet they were refined for, curr_deref is
  // false, and the cycle cannot close -- while the level cap still bounds block growth
  // exactly as before. Cost: one extra par_reduce on blocks at/above the cap.
  Real lbmin = std::numeric_limits<Real>::max();
  if (mhd) {
    pmb->par_reduce(
        "jeans_ni: current sheet", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int k, const int j, const int i, Real &llbmin) {
          if (w(IDN, k, j, i) < curr_rho_thresh) return; // envelope: skip (leave llbmin huge)
          // central differences (ghosts valid: prim is FillGhost)
          const Real dBz_dy = (w(IB3, k, j + 1, i) - w(IB3, k, j - 1, i)) / (2.0 * dy);
          const Real dBy_dz = (w(IB2, k + 1, j, i) - w(IB2, k - 1, j, i)) / (2.0 * dz);
          const Real dBx_dz = (w(IB1, k + 1, j, i) - w(IB1, k - 1, j, i)) / (2.0 * dz);
          const Real dBz_dx = (w(IB3, k, j, i + 1) - w(IB3, k, j, i - 1)) / (2.0 * dx);
          const Real dBy_dx = (w(IB2, k, j, i + 1) - w(IB2, k, j, i - 1)) / (2.0 * dx);
          const Real dBx_dy = (w(IB1, k, j + 1, i) - w(IB1, k, j - 1, i)) / (2.0 * dy);
          const Real Jx = dBz_dy - dBy_dz;
          const Real Jy = dBx_dz - dBz_dx;
          const Real Jz = dBy_dx - dBx_dy;
          const Real Jmag = std::sqrt(Jx * Jx + Jy * Jy + Jz * Jz);
          const Real Bmag = std::sqrt(w(IB1, k, j, i) * w(IB1, k, j, i) +
                                      w(IB2, k, j, i) * w(IB2, k, j, i) +
                                      w(IB3, k, j, i) * w(IB3, k, j, i));
          // L_B/dx = (|B|/|J|)/dx ; guard |J|~0 (uniform field) => huge => never triggers.
          const Real LB_cells = (Bmag / (Jmag + 1.0e-30)) / dx;
          llbmin = std::min(llbmin, LB_cells);
        },
        Kokkos::Min<Real>(lbmin));
  }

  const bool jeans_refine = (njmin < njeans);
  const bool jeans_deref = (njmin > 2.5 * njeans);
  // Level cap applied HERE (see the B1 note above): the sheet criterion may not DRIVE
  // refinement at or beyond curr_max_level, but its derefinement veto is always live.
  const bool curr_refine = (lbmin < curr_nsheet) && (block_level < curr_max_level);
  const bool curr_deref = (lbmin > 2.5 * curr_nsheet);

  if (jeans_refine || curr_refine) return parthenon::AmrTag::refine;
  if (jeans_deref && curr_deref) return parthenon::AmrTag::derefine;
  return parthenon::AmrTag::same;
}

} // namespace nonideal
} // namespace refinement
