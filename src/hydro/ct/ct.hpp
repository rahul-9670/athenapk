//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2020-2025, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
//! \file ct.hpp
//! \brief Constrained Transport (CT) magnetic-field update for AthenaPK.
//!
//! Phase 2 increment 1 (ideal MHD, single level). See docs/CT_DESIGN.md.
//!
//! Data model
//! ----------
//! On the CT path (``<hydro> divergence_control = ct``) a face-centered magnetic
//! field ``Bf`` (Metadata::Face, Independent, WithFluxes, FillGhost) is the
//! *primary* magnetic variable. The cell-centered components IB1..IB3 inside the
//! ``cons``/``prim`` containers become a *projection* of ``Bf`` (arithmetic
//! average of the two bounding faces) that is refreshed every substage, so every
//! existing consumer of cell-centered B (EOS, HLLD wave speeds, non-ideal eta,
//! gravity, radiation, analysis) is untouched.
//!
//! Induction update
//! ----------------
//!  1. CT_AssembleEMF : build the edge-centered EMF ``E = -v x B`` from the
//!     ghost-filled primitive velocity/field (Balsara-Spicer arithmetic average
//!     of the four cell-centered EMFs surrounding each edge) and store it in the
//!     edge-flux slots of ``Bf`` (Bf.flux(E1/E2/E3)).
//!  2. CT_UpdateBf  : advance the face field by the discrete curl of the edge EMF
//!     (generalized Stokes), combined with the VL2 low-storage weights.
//!  3. CT_ProjectBfToCC : write IB1..IB3 = 1/2 (Bf_i + Bf_{i+1}) into cons.
//!
//! The discrete curl telescopes so that the cell-centered *face* divergence
//! div B = sum_d (Bf_{d,+} - Bf_{d,-})/dx_d is preserved at round-off for any
//! single-valued edge EMF (this is the CT guarantee; CT_MaxRelFaceDivB measures
//! it). The sign convention below was derived independently and cross-checked
//! against Parthenon's generalized-Stokes face-curl indicator
//! (example/fine_advection/stokes.hpp): dB_x/dt = -(d_y E_z - d_z E_y), cyclic.
//========================================================================================
#ifndef HYDRO_CT_CT_HPP_
#define HYDRO_CT_CT_HPP_

#include <string>
#include <utility>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

using namespace parthenon::package::prelude;
using namespace parthenon::driver::prelude;

namespace Hydro {
namespace CT {

// Face-centered magnetic field. A single Metadata::Face field carries the three
// face-normal components accessed through TE::F1/F2/F3. The HDF5/restart name is
// "Bf".
struct Bf : public parthenon::variable_names::base_t<false> {
  template <class... Ts>
  KOKKOS_INLINE_FUNCTION Bf(Ts &&...args)
      : parthenon::variable_names::base_t<false>(std::forward<Ts>(args)...) {}
  static std::string name() { return "Bf"; }
};

// Assemble the edge-centered EMF (E = -v x B) from the primitive state and store
// it in the edge-flux slots of Bf. Reads "prim"; writes Bf.flux(E1,E2,E3).
//
// Two averaging schemes select the edge EMF value (div B is round-off for either,
// since the curl of any single-valued edge EMF telescopes):
//   CT_AssembleEMF       -- Balsara & Spicer (1999) arithmetic average of the four
//                           cell-centered EMFs around the edge (increment 1). Simple,
//                           but non-upwinded / under-dissipative on grid-aligned modes.
//   CT_AssembleEMF_GS05  -- Gardiner & Stone (2005) upwind CT (increment 3). The edge
//                           EMF is built from the HLLD *face* fluxes of transverse B
//                           (which carry the Riemann dissipation) plus a contact-mode-
//                           upwinded gradient correction toward the cell-centered
//                           reference EMF. This is Athena++'s scheme (src/field/
//                           calculate_corner_e.cpp), chosen to tighten the cross-code
//                           comparison. Reads "prim" + the transverse-B and mass face
//                           fluxes of "cons"; dt sets the upwind-switch sharpness.
TaskStatus CT_AssembleEMF(MeshData<Real> *md);
TaskStatus CT_AssembleEMF_GS05(MeshData<Real> *md, const Real dt);

// Non-ideal (resistive) edge EMF for CT (Phase 2, increment 4). Adds the Ohmic EMF
// E_ohm = eta * J, with J = curl(Bf) evaluated *edge-centered* directly from the face
// field, onto the ideal edge EMF already stored in Bf.flux(E1,E2,E3). Because
// CT_UpdateBf advances dB/dt = -curl(E), accumulating E += eta*J reproduces the
// resistive induction dB/dt = -curl(eta curl B) = eta grad^2 B (const eta) while
// keeping div B at round-off (curl of a single-valued edge EMF telescopes).
//
// This is the CT counterpart of the *induction* half of OhmicDiffFluxIsoFixed
// (resistivity.cpp): under divergence_control=ct that kernel's cons.flux(IBn)
// induction deposits are gated OFF (so GS05 does not double-count them), while its
// cons.flux(IEN) resistive-Poynting/heating term stays on the finite-volume energy
// flux. Must run AFTER CT_AssembleEMF[_GS05] (ideal EMF present) and BEFORE the
// flux-correction round (so the C-F reflux restricts the non-ideal edge EMF too).
// eta at the edge: fixed coeff -> uniform; ionization/cache -> average of the four
// (E3) / (E1,E2) cells bounding the edge. Unsplit only (STS+CT is incompatible).
TaskStatus CT_AddOhmicEMF(MeshData<Real> *md);

// Non-ideal (ambipolar) edge EMF for CT (Phase 2, increment 4; edge-current design fixed
// 2026-07-28). Adds the perpendicular-current EMF E_A = eta_A (J - (J.bhat) bhat) onto the
// ideal edge EMF in Bf.flux(E*). Unlike the Ohmic term (whose EMF is eta*J, a single curl
// component that lives naturally on the edge), the AD EMF needs the *full* current and field
// vectors at the edge. COMPACT edge-current construction (mirrors CT_AddHallEMF): the edge's
// own-direction current component comes from the tight 1dx curl of the face field Bf (same
// stencil Ohmic uses); the two transverse components are interpolated from neighboring
// SAME-TYPE edges (each already tight), and B is interpolated from Bf with short local
// averages. eta at the edge is a 4-cell average (NonidealEdgeEta), matching Ohmic/Hall.
// SUPERSEDES an earlier "assemble the perp-current EMF at 4 cell-prim-based faces (mixing a
// tight 1dx and a wide 2dx stencil for the same current derivative), then arithmetic-average
// to the edge" design: a first-principles comparison against the validated GLM face-EMF
// showed that design under-diffuses a diffusivity feature sharp on the grid scale (~1-2 dx,
// e.g. a deeply-refined AMR current sheet) by up to ~13% per evaluation -- a per-step deficit
// that compounds over many diffusion times into an order-of-magnitude under-diffusion (the
// flagship CT core-edge ME/E runaway that blocks first-core formation under CT+AD). See
// runs/ct_tests/diffusion_ad_bump_*.in and the AMBIPOLAR_CT_UNDERDIFFUSION note. div B stays
// at round-off (the edge value is still single-valued, so its curl telescopes); the uniform-
// eta damped-Alfven gate (diffusion_ad_ct.in) is unaffected (0.055% error, unchanged). As
// with Ohmic, the matching cons.flux(IBn) induction deposits in AmbipolarDiffFluxIsoFixed are
// gated OFF under CT while its cons.flux(IEN) ambipolar-Poynting term stays on the
// finite-volume energy flux. Must run after CT_AssembleEMF[_GS05] and before the
// flux-correction round. Unsplit only.
TaskStatus CT_AddAmbipolarEMF(MeshData<Real> *md);

// Non-ideal (Hall) edge EMF for CT (Phase 2, increment 4). Adds the dispersive Hall EMF
// E_H = eta_H (J x B)/|B| (plus the optional Ohmic stabilizer floor eta_O J) onto the ideal
// edge EMF, using the same COMPACT edge-current construction as CT_AddAmbipolarEMF above (the
// original template this fix was ported from: face-EMF-averaged Hall had the analogous
// whistler-dispersion-corrupting dilution and was replaced with this design first). Hall is
// dispersive (whistler) and unsplit-only under CT (RKL2+CT is forbidden), so the whistler
// part and the floor are both applied here. cons.flux(IBn) induction deposit in
// HallDiffFluxIsoFixed is gated OFF under CT; cons.flux(IEN) Poynting term stays on the FV
// energy flux. Must run after CT_AssembleEMF[_GS05] and before the flux-correction round.
TaskStatus CT_AddHallEMF(MeshData<Real> *md);

// Diffusive Poynting energy flux (S = E x B) built from the SAME edge EMF in Bf.flux(E*) that
// drives the CT induction, deposited into cons.flux(IEN). Replaces the face-based IEN deposits
// of resistivity.cpp/ambipolar.cpp, whose EMF comes from a different stencil -- a mismatch that
// lands directly in the recovered internal energy e = E - KE - ME. RKL2/STS path only (Bf.flux
// must hold the diffusive-only EMF); no-op unless "ct_edge_poynting". Must run after
// CT_AddOhmicEMF/CT_AddAmbipolarEMF and before the flux-correction round. See ct.cpp.
TaskStatus CT_AddDiffusivePoynting(MeshData<Real> *md);

// VL2 low-storage face update: Bf_base = gam0*Bf_base + gam1*Bf_u1
//                                       + beta_dt * curl(edge EMF stored in Bf.flux).
TaskStatus CT_UpdateBf(MeshData<Real> *md_base, MeshData<Real> *md_u1, const Real gam0,
                       const Real gam1, const Real beta_dt);

// ---- RKL2 super-time-stepping of the diffusive induction under CT (increment 4) --------
// The GLM STS (AddSTSTasks) super-time-steps the parabolic diffusion on the cell-centered
// cons via flux-divergence. Under CT the cons.flux(IBn) induction deposits are gated off,
// so the cons recurrence advances only the gas energy (resistive/AD Poynting) and leaves B
// inert; the divergence-free induction is super-time-stepped SEPARATELY on the face field
// Bf here, with the RKL2 operator M(Bf) = -curl(E_diff) (E_diff = the Ohmic+AD edge EMF).
// Because every sub-stage is a linear combination of Bf states plus curls of single-valued
// edge EMFs, div B stays at round-off through the whole super-step.

// Zero the edge-EMF flux slots Bf.flux(E1/E2/E3) so CT_AddOhmicEMF/CT_AddAmbipolarEMF (which
// accumulate with +=) build the *diffusive-only* edge EMF for the STS operator.
TaskStatus CT_ZeroEMF(MeshData<Real> *md);

// M-operator: write md_out's face field to the discrete curl of the edge EMF currently in
// md_emf's Bf.flux, i.e. Bf_out = -curl(E_diff). Used to precompute MY0 = M(Y0) once.
TaskStatus CT_CurlEMFToBf(MeshData<Real> *md_emf, MeshData<Real> *md_out);

// RKL2 first sub-step on Bf: base.Bf = Y0.Bf + mu_tilde_1*tau*MY0.Bf ; Yjm2.Bf = Y0.Bf.
TaskStatus CT_RKL2FirstBf(MeshData<Real> *md_Y0, MeshData<Real> *md_base,
                          MeshData<Real> *md_Yjm2, MeshData<Real> *md_MY0, const int s_rkl,
                          const Real tau);

// RKL2 sub-step j>=2 on Bf. MYjm1 = -curl(edge EMF in base.Bf.flux) is formed inline; then
// base.Bf = mu*base.Bf + nu*Yjm2.Bf + (1-mu-nu)*Y0.Bf + mu_tilde*tau*MYjm1
//           + gamma_tilde*tau*MY0.Bf, with Yjm2.Bf shuffled to the old base.Bf first.
TaskStatus CT_RKL2OtherBf(MeshData<Real> *md_Y0, MeshData<Real> *md_base,
                          MeshData<Real> *md_Yjm2, MeshData<Real> *md_MY0, const Real mu_j,
                          const Real nu_j, const Real mu_tilde_j, const Real gamma_tilde_j,
                          const Real tau);

// Project the face field onto the cell-centered slots IB1..IB3 of "cons".
TaskStatus CT_ProjectBfToCC(MeshData<Real> *md);

// History reductions over the per-cell projection diagnostic "ct.dEint" (filled by
// CT_ProjectBfToCC when hydro/ct_proj_diag=true; 0 otherwise). dEint_rel is the RELATIVE
// internal-energy change the face->cell projection imposes, measured BEFORE the eint guard
// acts -- i.e. the raw magnitude of the E-vs-ME bookkeeping mismatch. Min gives the damaging
// (cooling) direction; MaxAbs gives the largest excursion either way.
Real CT_ProjEintMin(MeshData<Real> *md);
Real CT_ProjEintMaxAbs(MeshData<Real> *md);

// History diagnostic: max over interior cells of |div B|_face * dx / |B|.
// For CT this stays at round-off; for GLM it grows at truncation level.
// NOTE: the /|B| normalization blows up in near-zero-field regions (e.g. a field-loop
// ambient), so a small round-off *absolute* divergence can read as a large *relative*
// one there -- use CT_MaxAbsFaceDivB to disambiguate a metric artifact from real growth.
Real CT_MaxRelFaceDivB(MeshData<Real> *md);

// History diagnostic: max over interior cells of |div B|_face * dx (absolute, un-normalized).
// The true CT invariant: div(curl E)=0 identically, so for any single-valued edge EMF this
// stays at round-off regardless of the field magnitude or the EMF averaging scheme.
Real CT_MaxAbsFaceDivB(MeshData<Real> *md);

} // namespace CT
} // namespace Hydro

#endif // HYDRO_CT_CT_HPP_
