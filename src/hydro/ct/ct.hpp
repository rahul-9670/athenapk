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
TaskStatus CT_AssembleEMF(MeshData<Real> *md);

// VL2 low-storage face update: Bf_base = gam0*Bf_base + gam1*Bf_u1
//                                       + beta_dt * curl(edge EMF stored in Bf.flux).
TaskStatus CT_UpdateBf(MeshData<Real> *md_base, MeshData<Real> *md_u1, const Real gam0,
                       const Real gam1, const Real beta_dt);

// Project the face field onto the cell-centered slots IB1..IB3 of "cons".
TaskStatus CT_ProjectBfToCC(MeshData<Real> *md);

// History diagnostic: max over interior cells of |div B|_face * dx / |B|.
// For CT this stays at round-off; for GLM it grows at truncation level.
Real CT_MaxRelFaceDivB(MeshData<Real> *md);

} // namespace CT
} // namespace Hydro

#endif // HYDRO_CT_CT_HPP_
