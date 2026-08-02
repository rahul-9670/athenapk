//========================================================================================
// AthenaPK - angular-momentum history diagnostics (VALIDATION_PLAN.md WP-4).
//
// WHY THIS EXISTS. The production history file carries 31 columns and NOT ONE of them is an
// angular momentum. For a magnetic-braking paper L(t) and the torque budget ARE the result,
// so every statement about braking in the current output is an inference from rho_max and
// ME/E rather than a measurement. This file adds the measurement.
//
// Gated on hydro/angmom_diag (default OFF). Pure read-only volume/surface reductions that
// write nothing back into the state => bit-identical when off.
//
// UNITS: Heaviside-Lorentz code units throughout (CLAUDE.md). The Lorentz force density is
// J x B with no 4pi/c, and J = curl B by the SAME 2*dx cell-centered central difference the
// mag_diag dissipation columns and the jeans_nonideal refinement criterion use -- so the
// refinement trigger, the dissipation diagnostic and the torque budget all see one current.
//
// ORIGIN. Moments are taken about the CENTRE OF THE BOX, computed from mesh_size as
// 0.5*(xmin+xmax) per direction, NOT about the origin of coordinates and NOT about the
// instantaneous centre of mass. For the FHC decks the box is symmetric so the box centre is
// (0,0,0) and this agrees with the z-axis convention mag_diag's MEtor/MEpol already use.
// The origin actually used is echoed in the rank-0 banner so it can never be ambiguous.
//
// Columns (all UserHistoryOperation::sum):
//   am-Lx/Ly/Lz        = int rho (r x v) dV          angular momentum
//   am-FTx/FTy/FTz     = surf (r x (T.nhat)) dA      TOTAL L flux out -- THE BUDGET TERM
//   am-Tgravx/y/z      = int rho r x (-grad Phi) dV  gravitational torque
//   am-FLx/FLy/FLz     = surf rho (r x v)(v.nhat) dA advective part of am-FT (interpretive)
//   am-Tmagx/y/z       = int r x (J x B) dV          volume Lorentz torque (interpretive)
//
// with the momentum stress tensor (Heaviside-Lorentz)
//
//     T_ij = rho v_i v_j + (P + B^2/2) delta_ij - B_i B_j
//
// so the budget to close is
//
//     d/dt (am-Lx) = -am-FTx + am-Tgravx
//
// The sign convention on the flux columns is OUTWARD-POSITIVE: am-FT* > 0 means angular
// momentum is LEAVING the box, hence the minus sign above.
//
// DO NOT ADD am-Tmag TO am-FT -- THEY DOUBLE-COUNT. The Lorentz torque is the volume form
// of exactly the Maxwell part of the surface term, since r x (J x B) = -div(r x T_Maxwell)
// up to the boundary contribution. The budget is complete in the SURFACE form above; the
// volume form is the different (equivalent) grouping
//
//     d/dt L = int r x (J x B) dV + int r x (-grad P) dV + am-Tgrav - am-FL
//
// and mixing terms from the two groupings is the easiest way to get a budget that appears
// not to close. am-Tmag is registered because the magnetic torque IS the physical quantity
// a braking paper reports -- but it is a diagnostic, not a budget term.
//
// MEASURED 2026-07-31, and the reason the extra columns exist: on the eos_smoke deck the
// naive budget d/dt L = am-Tmag - am-FL misses by more than two orders of magnitude.
// dLz/dt = +12.5 code, while max|am-Tmagz| = 8.1e-2 and max|am-FLz| = 7.5e-4. The
// dominant term is the SURFACE PRESSURE torque -surf P (r x nhat) dA, which under outflow
// boundaries is not remotely negligible and which the advective-only flux column omits
// entirely. That is why am-FT (full stress) exists and why it, not am-FL, is the budget
// term.
//
// WHY THE FLUX COLUMNS ARE NOT OPTIONAL. Production runs use OUTFLOW fluid boundaries, so
// angular momentum is NOT a conserved quantity of the simulation by construction -- L can
// and does leave through the faces. A bare L(t) curve therefore cannot distinguish
// "magnetic braking transported L out of the core" from "L advected off the grid", which
// are completely different physical claims. Only L together with the flux term can.
//
// ACCEPTANCE TEST -- and a CORRECTION to the one VALIDATION_PLAN.md WP-4 specifies.
//
// The plan asks for "L conserved to machine precision" on a rotating, B = 0, reflecting
// run. A Cartesian finite-volume scheme CANNOT deliver that, and expecting it would make
// the instrument look broken when it is not. Linear momentum telescopes exactly, because
// the flux leaving cell i through a face is exactly the flux entering cell i+1. Angular
// momentum does not: summation by parts leaves the residual
//
//     sum_faces [ A_x F^x_{rho v_y} dx  -  A_y F^y_{rho v_x} dy ]
//
// which cancels only if the numerical stress is symmetric, T_xy == T_yx. Those two numbers
// come from DIFFERENT Riemann problems at DIFFERENT locations (an x-face and a y-face), so
// they agree only to O(dx^2). Angular momentum is therefore conserved to TRUNCATION error,
// not to roundoff, on any Cartesian grid.
//
// The acceptance test is accordingly split into a part that is exact and a part that is not:
//
//   (a) INSTRUMENT correctness -- exact. Recompute int rho (r x v) dV from a saved snapshot
//       in numpy, sharing no code path with the kernel below, and require agreement to the
//       output precision (float32 snapshots, ~6 significant figures in the .hst).
//       Script: docs/validation/scripts/wp4_angmom_check.py. This is the real gate: it
//       validates the diagnostic independently of whether the SCHEME conserves anything.
//   (b) SCHEME conservation -- measured, not asserted. On a rotating, B = 0, reflecting,
//       turbulence-free run, report the drift and require it to FALL under refinement at
//       the scheme's order. A drift that does not converge is a defect; a drift that does
//       is arithmetic.
//
// BOUNDARY DETECTION. A cell is on a physical domain face when its own face coordinate
// coincides with the mesh extent, tested as |Xf - xmesh| < 1e-9 * dx. This deliberately
// avoids MeshBlock boundary flags, which are not reachable from inside the reduction
// kernel. A corner cell sits on two or three faces and contributes to each independently,
// which is correct -- the surface integral runs over all of them.
//
// GRAVITATIONAL TORQUE. am-Tgrav uses grad(Phi) by the same 2*dx central difference as the
// current, on the solver's own potential (grav.phi is FillGhost, so ghosts are valid). In
// the continuum the TOTAL self-gravitational torque on an isolated system vanishes exactly
// by Newton's third law -- the pair forces are central. Discretely it does not, and its
// measured size is a direct probe of how well the Poisson solve plus the gradient conserve
// angular momentum. A large am-Tgrav is a finding about the gravity solver, not noise.
//
// Registered only when self-gravity is on; without it there is no potential to difference.
//========================================================================================
#ifndef DIAGNOSTICS_ANGMOM_DIAG_HPP_
#define DIAGNOSTICS_ANGMOM_DIAG_HPP_

#include <parthenon/package.hpp>

using namespace parthenon::package::prelude;

namespace Diagnostics {

//! Which angular-momentum reduction to perform. The L* and FL* variants need only `prim`;
//! the Tmag* variants additionally require an MHD field (fluid == glmmhd), so the caller
//! (hydro.cpp) registers them only in that case.
enum class AngMom {
  Lx, // int rho (r x v)_x dV
  Ly,
  Lz,
  TmagX, // int (r x (J x B))_x dV      -- interpretive, NOT additive with FT*
  TmagY,
  TmagZ,
  FLx, // surf rho (r x v)_x (v.nhat) dA, outward-positive -- advective part of FT*
  FLy,
  FLz,
  FTx, // surf (r x (T.nhat))_x dA, outward-positive -- THE budget flux term
  FTy,
  FTz,
  TgravX, // int rho (r x (-grad Phi))_x dV
  TgravY,
  TgravZ,
  // ---- density-split L, registered ONLY when hydro/angmom_diag_rho_split > 0 -------------
  // WP-4 finding: global L is ill-conditioned in the same way WP-8's dissipation is --
  // 1.6% of the mass carries 64% of Lz, so the 10% "drift" is dominated by a tiny, rapidly
  // evolving high-density core while the bulk integral is well behaved. These columns make
  // that split explicit instead of inferring it from phdf post-processing.
  //
  // Gating discipline (identical to mag_diag_rho_split, and it is load-bearing): appending
  // history columns SHIFTS every downstream index, which would silently invalidate the
  // parsers reading .hst files already on disk. With rho_split = 0 (the default) NONE of
  // these are registered, so the column set -- and therefore every existing analysis
  // script -- is bit-identical. Always read column indices from the .hst header, never
  // hardcode them (docs/validation/scripts/wp4_angmom_check.py does this correctly).
  Lxhi, // int rho (r x v)_x dV over cells with rho >  hydro/angmom_diag_rho_split
  Lxlo, // ... over cells with rho <= hydro/angmom_diag_rho_split
  Lyhi,
  Lylo,
  Lzhi,
  Lzlo,
  Mhi // int rho dV over cells with rho > rho_split -- the denominator for "x% of the mass
      // carries y% of L". Without it the split columns cannot be normalised.
};

//! Angular-momentum reduction over interior cells (volume terms) or over the cells adjacent
//! to a physical domain face (flux terms). Moments are about the centre of the box; see the
//! file header.
Real AngMomReduce(MeshData<Real> *md, AngMom which);

} // namespace Diagnostics

#endif // DIAGNOSTICS_ANGMOM_DIAG_HPP_
