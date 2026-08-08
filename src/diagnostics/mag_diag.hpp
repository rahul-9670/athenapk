//========================================================================================
// AthenaPK - magnetic-transport (fossil-field) history diagnostics.
//
// FLAGSHIP AUDIT ITEM 2. Peak ME/E and rho_max are safety diagnostics: they cannot say
// WHERE the magnetic flux went. For a fossil-field result the decisive question is whether
// the PHYSICAL transport coefficients dominate a MEASURED numerical transport residual, so
// the run has to report the terms of the induction/energy budget, not just field strength.
//
// Everything here is a volume-summed history reduction over interior cells, gated on
// hydro/mag_diag (default OFF, and writing nothing back into the state => bit-identical
// when off). J = curl B by central differences on the cell-centered primitive B -- the SAME
// stencil the jeans_nonideal current-sheet refinement criterion uses, so the refinement
// trigger and the dissipation diagnostic see one current.
//
// UNITS: Heaviside-Lorentz code units throughout (P_mag = B^2/2, J = curl B with no 4pi/c;
// see CLAUDE.md). Ohmic heating per volume is eta_O |J|^2 and ambipolar heating per volume
// is eta_A |J_perp|^2 with J_perp = J - (J.bhat) bhat, matching the flux kernels.
//
// Columns (all UserHistoryOperation::sum; form ratios in analysis):
//   mag-Jsq      = int |J|^2 dV                total current content
//   mag-Hc       = int B.J dV                  CURRENT helicity -- gauge-invariant, unlike
//                                              magnetic helicity int A.B dV, so it needs no
//                                              vector potential and no gauge choice. It is
//                                              the sign/handedness tracer of the field.
//   mag-MEtor    = int B_phi^2/2 dV            toroidal ME about the z axis (the rotation
//   mag-MEpol    = int (B^2 - B_phi^2)/2 dV    axis and the initial B0z direction)
//   mag-dissO    = int eta_O |J|^2 dV          resolved OHMIC dissipation      (eta cache)
//   mag-dissA    = int eta_A |J_perp|^2 dV     resolved AMBIPOLAR dissipation  (eta cache)
//   mag-dissOcap = int eta_O |J|^2 dV over cells where the eta_O ceiling is ACTIVE
//   mag-dissAcap = int eta_A |J_perp|^2 dV over cells where the eta_A ceiling is ACTIVE
//
// The last two are the direct answer to "what fraction of the magnetic dissipation happens
// in capped cells" -- i.e. how much of the flux loss is a numerical stabilizer rather than
// physics. They require diffusion/cap_diag=true (see CapDiagIdx in diffusion.hpp).
//
//----------------------------------------------------------------------------------------
// WP-8 FINDING (2026-07-31): mag-dissO and mag-dissA are ILL-CONDITIONED as global scalars.
//
// Measured on the njeans ladder at matched epoch rho=1e-12: mag-Jsq falls 52% then 39% per
// refinement rung and mag-dissO is not even monotonic (-76.5%, then +40.5%), while mag-ME
// converges to 0.2%. The cause is NOT this kernel. Two candidate defects were tested and
// BOTH falsified:
//   * the 2*dx central-difference curl below was compared against the flux kernel's own
//     one-cell face difference (resistivity.cpp:139) on identical saved snapshots: they
//     agree to 1.5% / 0.7% / 0.1% and give the same non-convergence. The stencil is fine.
//   * eta-cache staleness: dissO moves opposite to Jsq in only 3-7% of history rows.
//
// The actual cause is CONCENTRATION. 90% of int |J|^2 dV comes from a volume fraction of
// 7.6e-7 (nj4) / 1.6e-6 (nj8) / 4.7e-6 (nj16), and 90% of the rho-weighted integrand (a
// proxy for eta_O |J|^2, since eta_O climbs steeply with density) comes from 3.6e-8 /
// 7.7e-8 / 1.3e-7 -- a few thousand cells out of ~3.5e7. These "volume integrals" are
// effectively POINT SAMPLES of the innermost core, taken in exactly the region whose
// resolution changes between rungs. A single number summing an integrand that spans ~7
// decades across the domain cannot converge until that tiny region is itself converged.
//
// Consequences for analysis:
//   * NEVER quote dissA/dissO as a physical ratio. It is a quotient of two integrals
//     dominated by DIFFERENT regions (ambipolar in the diffuse envelope, Ohmic in the dense
//     core), and dissO is exactly 0 at t=0 (uniform B0z => J=0), so the ratio is undefined
//     on the first row. The qualitative statement "AD dominates Ohmic by 4-5 decades" is
//     robust; the number is not.
//   * Use the density-split columns below instead: each is an integral over a region
//     defined by PHYSICS rather than by the grid, so the core and envelope budgets converge
//     (or fail to) independently and visibly.
//   * mag-Vhi / mag-dissOV / mag-dissAV report the volume actually carrying the
//     dissipation, so this ill-conditioning is never again invisible in the output.
//
// Density-split columns (registered only when hydro/mag_diag_rho_split > 0; default 0 keeps
// the OFF state and the original column set bit-identical):
//   mag-dissO-hi / mag-dissO-lo   Ohmic dissipation above / below the density split
//   mag-dissA-hi / mag-dissA-lo   ambipolar dissipation above / below the split
//   mag-Vhi                       volume with rho > split
//   mag-dissOV / mag-dissAV       volume of all cells with nonzero Ohmic / ambipolar heating
//========================================================================================
#ifndef DIAGNOSTICS_MAG_DIAG_HPP_
#define DIAGNOSTICS_MAG_DIAG_HPP_

#include <parthenon/package.hpp>

using namespace parthenon::package::prelude;

namespace Diagnostics {

//! Which volume-summed magnetic-transport integral to reduce.
enum class MagDiag {
  Jsq,      // int |J|^2 dV
  Hc,       // int B.J dV            (current helicity)
  MEtor,    // int B_phi^2/2 dV      (toroidal, about the z axis)
  MEpol,    // int B_pol^2/2 dV
  dissO,    // int eta_O |J|^2 dV
  dissA,    // int eta_A |J_perp|^2 dV
  dissOcap, // dissO restricted to cells where the eta_O ceiling is active
  dissAcap, // dissA restricted to cells where the eta_A ceiling is active
  // --- density-split variants (2026-07-31; see the WP-8 note in the file header) ---
  dissOhi, // dissO over cells with rho >  hydro/mag_diag_rho_split
  dissOlo, // dissO over cells with rho <= hydro/mag_diag_rho_split
  dissAhi, // dissA over cells with rho >  hydro/mag_diag_rho_split
  dissAlo, // dissA over cells with rho <= hydro/mag_diag_rho_split
  Vhi,     // int dV over cells with rho > split -- the volume actually carrying dissOhi
  // Concentration probes. NOT volumes: these are int q^2 dV with q the local heating-rate
  // density, so analysis can form the inverse participation ratio
  //     f_eff = (int q dV)^2 / (V_box * int q^2 dV)
  // = the fraction of the box that would carry the whole integral if q were uniform over
  // it. f_eff ~ 1 means a genuinely volume-filling integral; f_eff ~ 1e-7 means the number
  // is a point sample of a handful of cells and cannot converge under refinement.
  // (A literal "volume where q>0" was implemented first and is USELESS -- eta and J are
  // nonzero almost everywhere, so it just returns the box volume. Measured: 1.40608e5 on
  // the L=52 smoke deck, i.e. exactly 52^3.)
  dissOsq, // int (eta_O |J|^2)^2 dV
  dissAsq, // int (eta_A |J_perp|^2)^2 dV
  // --- WP-8 REOPENED 2026-08-06, ADDRESSED 2026-08-08: the CURRENT-SHEET split ------------
  // The density split above was built from a smoke-deck measurement showing 89.6 % of dissO
  // sits above rho = 1 code. That is right for dissO, which is weighted by eta_O and so is
  // core-dominated. It is WRONG for Jsq, which is unweighted: re-measured on the ladder at
  // production resolution, the low-density bin carries 97-99 % of Jsq AND keeps f_eff ~ 1e-7,
  // so `Jsq-lo` is simply the original pathology under a new name (its convergence history,
  // -51.6 %/-39.9 %, is the global one to within a percent). WP08_dissipation_nonconvergence.md
  // states the consequence directly: "Jsq needs a split on a current-sheet indicator, not on
  // density". A density threshold cannot separate grid-scale current sheets by construction,
  // because those sheets live in the diffuse envelope.
  //
  // The indicator is the DIMENSIONLESS grid-scale current
  //     s = |J| * dx_min / |B|
  // i.e. the fraction of the local field that reverses across one cell. s -> 1 means B flips
  // over a single zone: the current is at the grid scale and is a numerical-resolution
  // artefact as much as a physical structure. s << 1 means a current spread over many cells,
  // which is a resolved, physically meaningful sheet. This is the right variable precisely
  // because it is measured IN UNITS OF THE GRID -- the thing that changes between ladder
  // rungs -- whereas density is blind to it.
  //
  // Split at hydro/mag_diag_sheet_thresh (default 0 => these columns are NEVER REGISTERED and
  // the OFF state is bit-identical, same convention as the density split). Cells with |B| = 0
  // are counted as "sheet": s is then formally infinite, and a nonzero J with zero B is a
  // pure grid artefact by definition.
  Jsqsheet, // int |J|^2 dV over cells with s >  thresh   (grid-scale current)
  Jsqsmth,  // int |J|^2 dV over cells with s <= thresh   (resolved current)
  Vsheet,   // int dV      over cells with s >  thresh   -- the volume carrying Jsqsheet
  // Concentration probe for Jsq itself. Until now Jsq had NO `sq` companion, so unlike dissO
  // and dissA its f_eff could not be formed from the history file at all -- it had to be
  // recomputed offline from 986 GB of ladder snapshots, which is why the reopened WP-8 test
  // could only be run once and only for Jsq. With this column,
  //     f_eff(Jsq) = (int |J|^2 dV)^2 / (V_box * int |J|^4 dV)
  // is available every history row, for free, in every run.
  Jsqsq // int |J|^4 dV
};

//! Volume-summed magnetic-transport reduction over interior cells. `need_eta` variants
//! read the cached "nonideal_eta" field; the `*cap` variants additionally read
//! "diff.capdiag". The caller (hydro.cpp) only registers the variants whose inputs exist.
Real MagDiagReduce(MeshData<Real> *md, MagDiag which);

} // namespace Diagnostics

#endif // DIAGNOSTICS_MAG_DIAG_HPP_
