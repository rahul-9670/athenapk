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

// Non-ideal (ambipolar) edge EMF for CT (Phase 2, increment 4). Adds the perpendicular-
// current EMF E_A = eta_A (J - (J.bhat) bhat) onto the ideal edge EMF in Bf.flux(E*).
// Unlike the Ohmic term (whose EMF is eta*J, a single curl component that lives naturally
// on the edge), the AD EMF needs the *full* current and field vectors, so it is built the
// way Athena++'s FieldDiffusion::AddEMF does: the perp-current EMF is evaluated at cell
// FACES with the exact same stencils as the GLM path (AmbipolarDiffFluxIsoFixed) and the
// relevant component is arithmetic-averaged from the four faces bounding each edge (the
// same four-face pattern GS05 uses for the ideal base EMF). div B stays at round-off (the
// edge value is single-valued, so its curl telescopes) and the operator is identical to
// the validated GLM AD operator by construction. As with Ohmic, the matching cons.flux(IBn)
// induction deposits in AmbipolarDiffFluxIsoFixed are gated OFF under CT while its
// cons.flux(IEN) ambipolar-Poynting term stays on the finite-volume energy flux. Must run
// after CT_AssembleEMF[_GS05] and before the flux-correction round. Unsplit only.
TaskStatus CT_AddAmbipolarEMF(MeshData<Real> *md);

// Non-ideal (Hall) edge EMF for CT (Phase 2, increment 4). Adds the dispersive Hall EMF
// E_H = eta_H (J x B)/|B| (plus the optional Ohmic stabilizer floor eta_O J) onto the ideal
// edge EMF, built by the same four-face-average construction as CT_AddAmbipolarEMF (perp
// EMF evaluated at faces with the GLM stencils, averaged to edges). Hall is dispersive
// (whistler) and unsplit-only under CT (RKL2+CT is forbidden), so the whistler part and the
// floor are both applied here. cons.flux(IBn) induction deposit in HallDiffFluxIsoFixed is
// gated OFF under CT; cons.flux(IEN) Poynting term stays on the FV energy flux. Must run
// after CT_AssembleEMF[_GS05] and before the flux-correction round.
TaskStatus CT_AddHallEMF(MeshData<Real> *md);

// VL2 low-storage face update: Bf_base = gam0*Bf_base + gam1*Bf_u1
//                                       + beta_dt * curl(edge EMF stored in Bf.flux).
TaskStatus CT_UpdateBf(MeshData<Real> *md_base, MeshData<Real> *md_u1, const Real gam0,
                       const Real gam1, const Real beta_dt);

// Project the face field onto the cell-centered slots IB1..IB3 of "cons".
TaskStatus CT_ProjectBfToCC(MeshData<Real> *md);

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
