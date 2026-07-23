//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2020-2025, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
//! \file ct.cpp
//! \brief Constrained Transport magnetic-field update (Phase 2, increment 1).
//! See ct.hpp / docs/CT_DESIGN.md.
//========================================================================================

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

#include "../../main.hpp"
#include "../diffusion/diffusion.hpp"
#include "ct.hpp"

using namespace parthenon::package::prelude;

namespace Hydro {
namespace CT {

using parthenon::CellLevel;
using TE = parthenon::TopologicalElement;

//----------------------------------------------------------------------------------------
//! Cell-centered EMF components E = -v x B, evaluated from primitive (v, B).
//!   E_x = v_z B_y - v_y B_z ,  E_y = v_x B_z - v_z B_x ,  E_z = v_y B_x - v_x B_y
template <int comp, class Pack>
KOKKOS_INLINE_FUNCTION Real EmfCC(const Pack &prim, const int b, const int k, const int j,
                                  const int i) {
  const Real vx = prim(b, IV1, k, j, i);
  const Real vy = prim(b, IV2, k, j, i);
  const Real vz = prim(b, IV3, k, j, i);
  const Real bx = prim(b, IB1, k, j, i);
  const Real by = prim(b, IB2, k, j, i);
  const Real bz = prim(b, IB3, k, j, i);
  if (comp == 0) return vz * by - vy * bz; // E_x
  if (comp == 1) return vx * bz - vz * bx; // E_y
  return vy * bx - vx * by;                // E_z
}

//----------------------------------------------------------------------------------------
//! Assemble edge EMF (Balsara-Spicer arithmetic average of the 4 surrounding CC EMFs)
//! and store it in the edge-flux slots of Bf.
TaskStatus CT_AssembleEMF(MeshData<Real> *md) {
  const int ndim = md->GetMeshPointer()->ndim;
  auto prim = md->PackVariables(std::vector<std::string>{"prim"});
  static auto desc =
      parthenon::MakePackDescriptor<Bf>(md, {}, {parthenon::PDOpt::WithFluxes});
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();

  // E_z on E3 edges (corner in the x-y plane). Always present (needed in 2D & 3D).
  {
    IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E3);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_EMF_E3", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          pack.flux(b, TE::E3, Bf(), k, j, i) =
              0.25 * (EmfCC<2>(prim, b, k, j, i) + EmfCC<2>(prim, b, k, j, i - 1) +
                      EmfCC<2>(prim, b, k, j - 1, i) + EmfCC<2>(prim, b, k, j - 1, i - 1));
        });
  }

  if (ndim > 2) {
    // E_x on E1 edges (corner in the y-z plane).
    IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E1);
    IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E1);
    IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E1);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_EMF_E1", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          pack.flux(b, TE::E1, Bf(), k, j, i) =
              0.25 * (EmfCC<0>(prim, b, k, j, i) + EmfCC<0>(prim, b, k - 1, j, i) +
                      EmfCC<0>(prim, b, k, j - 1, i) + EmfCC<0>(prim, b, k - 1, j - 1, i));
        });
    // E_y on E2 edges (corner in the z-x plane).
    ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E2);
    jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E2);
    kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E2);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_EMF_E2", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          pack.flux(b, TE::E2, Bf(), k, j, i) =
              0.25 * (EmfCC<1>(prim, b, k, j, i) + EmfCC<1>(prim, b, k - 1, j, i) +
                      EmfCC<1>(prim, b, k, j, i - 1) + EmfCC<1>(prim, b, k - 1, j, i - 1));
        });
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Constrained-transport upwind weight (Gardiner & Stone 2005; identical to Athena++
//! Hydro::GetWeightForCT). Returns w in [0,1]: w->1 when the mass flux is strongly
//! positive (flow from the "left"/lower-index cell, take its gradient), w->0 when
//! strongly negative, w=1/2 at a stagnation face. The 1024 dt/dx scaling saturates the
//! switch except within a narrow transonic band around zero contact speed.
KOKKOS_INLINE_FUNCTION Real WeightForCT(const Real dflx, const Real rhol, const Real rhor,
                                        const Real dx, const Real dt) {
  const Real v_over_c = 1024.0 * dt * dflx / (dx * (rhol + rhor) + 1.0e-300);
  return 0.5 + std::max(-0.5, std::min(0.5, v_over_c));
}

//----------------------------------------------------------------------------------------
//! Gardiner & Stone (2005) upwind CT edge EMF. Mirrors Athena++'s ComputeCornerE
//! (src/field/calculate_corner_e.cpp): the corner EMF is 1/4 of the sum of the four
//! bounding face EMFs plus four contact-upwinded gradient corrections (face-EMF minus
//! cell-EMF), which reintroduces the directional dissipation the plain arithmetic
//! average discards. Face EMFs come from the HLLD transverse-B fluxes:
//!   E_z from x-face = -F_x(B_y) = -cons.flux(X1,IB2);  E_z from y-face = +F_y(B_x) = +cons.flux(X2,IB1)
//!   E_y from x-face = +cons.flux(X1,IB3);              E_y from z-face = -cons.flux(X3,IB1)
//!   E_x from y-face = -cons.flux(X2,IB3);              E_x from z-face = +cons.flux(X3,IB2)
//! (all cyclic; signs from F_i(B_j)=-E_k for (i,j,k) cyclic, see docs/CT_DESIGN.md 1.2).
//! The transverse fluxes are available one layer into the ghost zone because the flux
//! kernel computes X1 flux over j in [js-1,je+1], X2 flux over i in [is-1,ie+1], etc.,
//! so the wider GS05 stencil (and same-level block-boundary consistency) is satisfied
//! without any extra exchange.
TaskStatus CT_AssembleEMF_GS05(MeshData<Real> *md, const Real dt) {
  const int ndim = md->GetMeshPointer()->ndim;
  auto prim = md->PackVariables(std::vector<std::string>{"prim"});
  auto cons_pack = md->PackVariablesAndFluxes(std::vector<std::string>{"cons"},
                                              std::vector<std::string>{"cons"});
  static auto desc =
      parthenon::MakePackDescriptor<Bf>(md, {}, {parthenon::PDOpt::WithFluxes});
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();
  const bool three_d = ndim > 2;

  // ---- E3 on E3 edges (x-y corner). Needed in 2D and 3D. ----
  {
    IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E3);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_EMF_GS05_E3", parthenon::DevExecSpace(), 0, nb - 1,
        kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(b);
          const auto &cons = cons_pack(b);
          const Real dx1 = c.Dxc<X1DIR>(k, j, i);
          const Real dx2 = c.Dxc<X2DIR>(k, j, i);
          // face EMFs (E_z): x-face = -F_x(By); y-face = +F_y(Bx)
          const Real e3x1f_j = -cons.flux(X1DIR, IB2, k, j, i);
          const Real e3x1f_jm = -cons.flux(X1DIR, IB2, k, j - 1, i);
          const Real e3x2f_i = cons.flux(X2DIR, IB1, k, j, i);
          const Real e3x2f_im = cons.flux(X2DIR, IB1, k, j, i - 1);
          // cell-centered reference E_z = -(v x B)_z
          const Real cc_jj_ii = EmfCC<2>(prim, b, k, j, i);
          const Real cc_jj_im = EmfCC<2>(prim, b, k, j, i - 1);
          const Real cc_jm_ii = EmfCC<2>(prim, b, k, j - 1, i);
          const Real cc_jm_im = EmfCC<2>(prim, b, k, j - 1, i - 1);
          // contact-upwind weights from the mass flux across each bounding face
          const Real wx1_jm = WeightForCT(cons.flux(X1DIR, IDN, k, j - 1, i),
                                          prim(b, IDN, k, j - 1, i - 1),
                                          prim(b, IDN, k, j - 1, i), dx1, dt);
          const Real wx1_j =
              WeightForCT(cons.flux(X1DIR, IDN, k, j, i), prim(b, IDN, k, j, i - 1),
                          prim(b, IDN, k, j, i), dx1, dt);
          const Real wx2_im = WeightForCT(cons.flux(X2DIR, IDN, k, j, i - 1),
                                          prim(b, IDN, k, j - 1, i - 1),
                                          prim(b, IDN, k, j, i - 1), dx2, dt);
          const Real wx2_i =
              WeightForCT(cons.flux(X2DIR, IDN, k, j, i), prim(b, IDN, k, j - 1, i),
                          prim(b, IDN, k, j, i), dx2, dt);
          const Real de_l2 = (1.0 - wx1_jm) * (e3x2f_i - cc_jm_ii) +
                             wx1_jm * (e3x2f_im - cc_jm_im);
          const Real de_r2 =
              (1.0 - wx1_j) * (e3x2f_i - cc_jj_ii) + wx1_j * (e3x2f_im - cc_jj_im);
          const Real de_l1 = (1.0 - wx2_im) * (e3x1f_j - cc_jj_im) +
                             wx2_im * (e3x1f_jm - cc_jm_im);
          const Real de_r1 =
              (1.0 - wx2_i) * (e3x1f_j - cc_jj_ii) + wx2_i * (e3x1f_jm - cc_jm_ii);
          pack.flux(b, TE::E3, Bf(), k, j, i) =
              0.25 * (de_l1 + de_r1 + de_l2 + de_r2 + e3x2f_im + e3x2f_i + e3x1f_jm +
                      e3x1f_j);
        });
  }

  if (three_d) {
    // ---- E1 on E1 edges (y-z corner). E_x: y-face=-F_y(Bz), z-face=+F_z(Bx) ----
    IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E1);
    IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E1);
    IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E1);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_EMF_GS05_E1", parthenon::DevExecSpace(), 0, nb - 1,
        kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(b);
          const auto &cons = cons_pack(b);
          const Real dx2 = c.Dxc<X2DIR>(k, j, i);
          const Real dx3 = c.Dxc<X3DIR>(k, j, i);
          const Real e1x2f_k = -cons.flux(X2DIR, IB3, k, j, i);
          const Real e1x2f_km = -cons.flux(X2DIR, IB3, k - 1, j, i);
          const Real e1x3f_j = cons.flux(X3DIR, IB2, k, j, i);
          const Real e1x3f_jm = cons.flux(X3DIR, IB2, k, j - 1, i);
          const Real cc_kk_jj = EmfCC<0>(prim, b, k, j, i);
          const Real cc_kk_jm = EmfCC<0>(prim, b, k, j - 1, i);
          const Real cc_km_jj = EmfCC<0>(prim, b, k - 1, j, i);
          const Real cc_km_jm = EmfCC<0>(prim, b, k - 1, j - 1, i);
          const Real wx2_km = WeightForCT(cons.flux(X2DIR, IDN, k - 1, j, i),
                                          prim(b, IDN, k - 1, j - 1, i),
                                          prim(b, IDN, k - 1, j, i), dx2, dt);
          const Real wx2_k =
              WeightForCT(cons.flux(X2DIR, IDN, k, j, i), prim(b, IDN, k, j - 1, i),
                          prim(b, IDN, k, j, i), dx2, dt);
          const Real wx3_jm = WeightForCT(cons.flux(X3DIR, IDN, k, j - 1, i),
                                          prim(b, IDN, k - 1, j - 1, i),
                                          prim(b, IDN, k, j - 1, i), dx3, dt);
          const Real wx3_j =
              WeightForCT(cons.flux(X3DIR, IDN, k, j, i), prim(b, IDN, k - 1, j, i),
                          prim(b, IDN, k, j, i), dx3, dt);
          const Real de_l3 = (1.0 - wx2_km) * (e1x3f_j - cc_km_jj) +
                             wx2_km * (e1x3f_jm - cc_km_jm);
          const Real de_r3 =
              (1.0 - wx2_k) * (e1x3f_j - cc_kk_jj) + wx2_k * (e1x3f_jm - cc_kk_jm);
          const Real de_l2 = (1.0 - wx3_jm) * (e1x2f_k - cc_kk_jm) +
                             wx3_jm * (e1x2f_km - cc_km_jm);
          const Real de_r2 =
              (1.0 - wx3_j) * (e1x2f_k - cc_kk_jj) + wx3_j * (e1x2f_km - cc_km_jj);
          pack.flux(b, TE::E1, Bf(), k, j, i) =
              0.25 * (de_l3 + de_r3 + de_l2 + de_r2 + e1x2f_km + e1x2f_k + e1x3f_jm +
                      e1x3f_j);
        });
    // ---- E2 on E2 edges (z-x corner). E_y: x-face=+F_x(Bz), z-face=-F_z(Bx) ... ----
    // face EMFs (E_y): x-face = +F_x(Bz)=+cons.flux(X1,IB3); z-face = -F_z(Bx)=-cons.flux(X3,IB1)
    ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E2);
    jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E2);
    kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E2);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_EMF_GS05_E2", parthenon::DevExecSpace(), 0, nb - 1,
        kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(b);
          const auto &cons = cons_pack(b);
          const Real dx1 = c.Dxc<X1DIR>(k, j, i);
          const Real dx3 = c.Dxc<X3DIR>(k, j, i);
          const Real e2x1f_k = cons.flux(X1DIR, IB3, k, j, i);
          const Real e2x1f_km = cons.flux(X1DIR, IB3, k - 1, j, i);
          const Real e2x3f_i = -cons.flux(X3DIR, IB1, k, j, i);
          const Real e2x3f_im = -cons.flux(X3DIR, IB1, k, j, i - 1);
          const Real cc_kk_ii = EmfCC<1>(prim, b, k, j, i);
          const Real cc_kk_im = EmfCC<1>(prim, b, k, j, i - 1);
          const Real cc_km_ii = EmfCC<1>(prim, b, k - 1, j, i);
          const Real cc_km_im = EmfCC<1>(prim, b, k - 1, j, i - 1);
          const Real wx1_km = WeightForCT(cons.flux(X1DIR, IDN, k - 1, j, i),
                                          prim(b, IDN, k - 1, j, i - 1),
                                          prim(b, IDN, k - 1, j, i), dx1, dt);
          const Real wx1_k =
              WeightForCT(cons.flux(X1DIR, IDN, k, j, i), prim(b, IDN, k, j, i - 1),
                          prim(b, IDN, k, j, i), dx1, dt);
          const Real wx3_im = WeightForCT(cons.flux(X3DIR, IDN, k, j, i - 1),
                                          prim(b, IDN, k - 1, j, i - 1),
                                          prim(b, IDN, k, j, i - 1), dx3, dt);
          const Real wx3_i =
              WeightForCT(cons.flux(X3DIR, IDN, k, j, i), prim(b, IDN, k - 1, j, i),
                          prim(b, IDN, k, j, i), dx3, dt);
          const Real de_l3 = (1.0 - wx1_km) * (e2x3f_i - cc_km_ii) +
                             wx1_km * (e2x3f_im - cc_km_im);
          const Real de_r3 =
              (1.0 - wx1_k) * (e2x3f_i - cc_kk_ii) + wx1_k * (e2x3f_im - cc_kk_im);
          const Real de_l1 = (1.0 - wx3_im) * (e2x1f_k - cc_kk_im) +
                             wx3_im * (e2x1f_km - cc_km_im);
          const Real de_r1 =
              (1.0 - wx3_i) * (e2x1f_k - cc_kk_ii) + wx3_i * (e2x1f_km - cc_km_ii);
          pack.flux(b, TE::E2, Bf(), k, j, i) =
              0.25 * (de_l3 + de_r3 + de_l1 + de_r1 + e2x3f_im + e2x3f_i + e2x1f_km +
                      e2x1f_k);
        });
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! VL2 low-storage face update via the discrete curl of the edge EMF.
//!   dB_x/dt = -(d_y E_z - d_z E_y) , and cyclic. (E stored in Bf.flux(E1/E2/E3).)
TaskStatus CT_UpdateBf(MeshData<Real> *md_base, MeshData<Real> *md_u1, const Real gam0,
                       const Real gam1, const Real beta_dt) {
  const int ndim = md_base->GetMeshPointer()->ndim;
  static auto desc =
      parthenon::MakePackDescriptor<Bf>(md_base, {}, {parthenon::PDOpt::WithFluxes});
  auto pb = desc.GetPack(md_base);
  auto pu = desc.GetPack(md_u1);
  const int nb = pb.GetNBlocks();
  const bool three_d = ndim > 2;

  // ---- F1 (B_x): dB_x/dt = -(E_z(j+1)-E_z(j))/dy + (E_y(k+1)-E_y(k))/dz ----
  {
    IndexRange ib = md_base->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::F1);
    IndexRange jb = md_base->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::F1);
    IndexRange kb = md_base->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::F1);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_Bf_F1", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &coords = pb.GetCoordinates(b);
          Real curl = -(pb.flux(b, TE::E3, Bf(), k, j + 1, i) -
                        pb.flux(b, TE::E3, Bf(), k, j, i)) /
                      coords.Dxc<X2DIR>(k, j, i);
          if (three_d) {
            curl += (pb.flux(b, TE::E2, Bf(), k + 1, j, i) -
                     pb.flux(b, TE::E2, Bf(), k, j, i)) /
                    coords.Dxc<X3DIR>(k, j, i);
          }
          pb(b, TE::F1, Bf(), k, j, i) = gam0 * pb(b, TE::F1, Bf(), k, j, i) +
                                         gam1 * pu(b, TE::F1, Bf(), k, j, i) +
                                         beta_dt * curl;
        });
  }
  // ---- F2 (B_y): dB_y/dt = -(E_x(k+1)-E_x(k))/dz + (E_z(i+1)-E_z(i))/dx ----
  {
    IndexRange ib = md_base->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::F2);
    IndexRange jb = md_base->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::F2);
    IndexRange kb = md_base->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::F2);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_Bf_F2", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &coords = pb.GetCoordinates(b);
          Real curl = (pb.flux(b, TE::E3, Bf(), k, j, i + 1) -
                       pb.flux(b, TE::E3, Bf(), k, j, i)) /
                      coords.Dxc<X1DIR>(k, j, i);
          if (three_d) {
            curl -= (pb.flux(b, TE::E1, Bf(), k + 1, j, i) -
                     pb.flux(b, TE::E1, Bf(), k, j, i)) /
                    coords.Dxc<X3DIR>(k, j, i);
          }
          pb(b, TE::F2, Bf(), k, j, i) = gam0 * pb(b, TE::F2, Bf(), k, j, i) +
                                         gam1 * pu(b, TE::F2, Bf(), k, j, i) +
                                         beta_dt * curl;
        });
  }
  // ---- F3 (B_z): dB_z/dt = -(E_y(i+1)-E_y(i))/dx + (E_x(j+1)-E_x(j))/dy (3D only) ----
  if (three_d) {
    IndexRange ib = md_base->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::F3);
    IndexRange jb = md_base->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::F3);
    IndexRange kb = md_base->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::F3);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_Bf_F3", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &coords = pb.GetCoordinates(b);
          Real curl = -(pb.flux(b, TE::E2, Bf(), k, j, i + 1) -
                        pb.flux(b, TE::E2, Bf(), k, j, i)) /
                          coords.Dxc<X1DIR>(k, j, i) +
                      (pb.flux(b, TE::E1, Bf(), k, j + 1, i) -
                       pb.flux(b, TE::E1, Bf(), k, j, i)) /
                          coords.Dxc<X2DIR>(k, j, i);
          pb(b, TE::F3, Bf(), k, j, i) = gam0 * pb(b, TE::F3, Bf(), k, j, i) +
                                         gam1 * pu(b, TE::F3, Bf(), k, j, i) +
                                         beta_dt * curl;
        });
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Project the face field onto cell-centered IB1..IB3 in cons.
TaskStatus CT_ProjectBfToCC(MeshData<Real> *md) {
  const int ndim = md->GetMeshPointer()->ndim;
  auto cons = md->PackVariables(std::vector<std::string>{"cons"});
  static auto desc = parthenon::MakePackDescriptor<Bf>(md);
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();
  const int koff = (ndim > 2) ? 1 : 0;
  const int joff = (ndim > 1) ? 1 : 0;

  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "CT_ProjectBfToCC", parthenon::DevExecSpace(), 0, nb - 1,
      kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        cons(b, IB1, k, j, i) = 0.5 * (pack(b, TE::F1, Bf(), k, j, i) +
                                       pack(b, TE::F1, Bf(), k, j, i + 1));
        cons(b, IB2, k, j, i) = 0.5 * (pack(b, TE::F2, Bf(), k, j, i) +
                                       pack(b, TE::F2, Bf(), k, j + joff, i));
        cons(b, IB3, k, j, i) = 0.5 * (pack(b, TE::F3, Bf(), k, j, i) +
                                       pack(b, TE::F3, Bf(), k + koff, j, i));
      });
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! max over interior cells of |div B|_face * dx / |B_cc|.
Real CT_MaxRelFaceDivB(MeshData<Real> *md) {
  const int ndim = md->GetMeshPointer()->ndim;
  auto cons = md->PackVariables(std::vector<std::string>{"cons"});
  static auto desc = parthenon::MakePackDescriptor<Bf>(md);
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();
  const bool three_d = ndim > 2;

  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  Real maxdiv = 0.0;
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "CT_MaxRelFaceDivB", parthenon::DevExecSpace(),
      0, nb - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &lmax) {
        const auto &coords = pack.GetCoordinates(b);
        const Real dx = coords.Dxc<X1DIR>(k, j, i);
        Real divb = (pack(b, TE::F1, Bf(), k, j, i + 1) - pack(b, TE::F1, Bf(), k, j, i)) /
                    coords.Dxc<X1DIR>(k, j, i);
        divb += (pack(b, TE::F2, Bf(), k, j + 1, i) - pack(b, TE::F2, Bf(), k, j, i)) /
                coords.Dxc<X2DIR>(k, j, i);
        if (three_d) {
          divb += (pack(b, TE::F3, Bf(), k + 1, j, i) - pack(b, TE::F3, Bf(), k, j, i)) /
                  coords.Dxc<X3DIR>(k, j, i);
        }
        const Real bmag =
            std::sqrt(SQR(cons(b, IB1, k, j, i)) + SQR(cons(b, IB2, k, j, i)) +
                      SQR(cons(b, IB3, k, j, i))) +
            1.0e-30;
        const Real rel = std::abs(divb) * dx / bmag;
        lmax = std::max(lmax, rel);
      },
      Kokkos::Max<Real>(maxdiv));
  return maxdiv;
}

//----------------------------------------------------------------------------------------
//! max over interior cells of |div B|_face * dx  (absolute, no |B| normalization).
Real CT_MaxAbsFaceDivB(MeshData<Real> *md) {
  const int ndim = md->GetMeshPointer()->ndim;
  static auto desc = parthenon::MakePackDescriptor<Bf>(md);
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();
  const bool three_d = ndim > 2;

  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  Real maxdiv = 0.0;
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "CT_MaxAbsFaceDivB", parthenon::DevExecSpace(),
      0, nb - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &lmax) {
        const auto &coords = pack.GetCoordinates(b);
        Real divb = (pack(b, TE::F1, Bf(), k, j, i + 1) - pack(b, TE::F1, Bf(), k, j, i)) /
                    coords.Dxc<X1DIR>(k, j, i);
        divb += (pack(b, TE::F2, Bf(), k, j + 1, i) - pack(b, TE::F2, Bf(), k, j, i)) /
                coords.Dxc<X2DIR>(k, j, i);
        if (three_d) {
          divb += (pack(b, TE::F3, Bf(), k + 1, j, i) - pack(b, TE::F3, Bf(), k, j, i)) /
                  coords.Dxc<X3DIR>(k, j, i);
        }
        lmax = std::max(lmax, std::abs(divb) * coords.Dxc<X1DIR>(k, j, i));
      },
      Kokkos::Max<Real>(maxdiv));
  return maxdiv;
}

//----------------------------------------------------------------------------------------
//! Cell-centered Ohmic diffusivity eta_O, from the cache if on, else evaluated from the
//! local primitive state (for ResistivityCoeff::fixed this ignores the state and returns
//! the uniform coefficient). Used to average eta to the edge in CT_AddOhmicEMF.
template <class Pack>
KOKKOS_INLINE_FUNCTION Real CellEtaO(const bool use_cache, const OhmicDiffusivity &ohm,
                                     const Pack &eta_pack, const Pack &prim, const int b,
                                     const int k, const int j, const int i) {
  if (use_cache) {
    return eta_pack(b, NonidealEtaIdx::O, k, j, i);
  }
  const Real bm = std::sqrt(SQR(prim(b, IB1, k, j, i)) + SQR(prim(b, IB2, k, j, i)) +
                            SQR(prim(b, IB3, k, j, i)));
  const Real rho = prim(b, IDN, k, j, i);
  const Real prs = prim(b, IPR, k, j, i);
  const int i_xe = ohm.XeIndex();
  const Real xe = (i_xe >= 0) ? prim(b, i_xe, k, j, i) : -1.0;
  return ohm.Get(bm, rho, prs / rho, xe);
}

//----------------------------------------------------------------------------------------
//! Non-ideal (Ohmic) edge EMF for CT: E += eta * J with J = curl(Bf) evaluated directly
//! from the staggered face field, so J is naturally edge-centered on the same E1/E2/E3
//! edges that carry the ideal EMF. See ct.hpp. Sign: CT_UpdateBf does dB/dt=-curl(E),
//! so E_ohm = eta*J = eta*curl(B) gives dB/dt = -curl(eta curl B) = eta grad^2 B.
TaskStatus CT_AddOhmicEMF(MeshData<Real> *md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto hydro_pkg = pmb->packages.Get("Hydro");
  // No-op unless Ohmic resistivity is active (CT diffusion is unsplit-only).
  if (hydro_pkg->Param<Resistivity>("resistivity") == Resistivity::none) {
    return TaskStatus::complete;
  }
  const int ndim = md->GetMeshPointer()->ndim;
  const bool three_d = ndim > 2;
  const auto &ohm_diff = hydro_pkg->Param<OhmicDiffusivity>("ohm_diff");
  const bool use_cache = hydro_pkg->Param<bool>("nonideal_eta_cache");

  auto prim = md->PackVariables(std::vector<std::string>{"prim"});
  // When the cache is off, "prim" is packed as a dummy to keep the capture valid; it is
  // never indexed on that branch (mirrors OhmicDiffFluxIsoFixed).
  const auto eta_pack =
      md->PackVariables(std::vector<std::string>{use_cache ? "nonideal_eta" : "prim"});
  static auto desc =
      parthenon::MakePackDescriptor<Bf>(md, {}, {parthenon::PDOpt::WithFluxes});
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();

  // ---- E3 edge (x-y corner, i-1/2,j-1/2,k). J3 = dBy/dx - dBx/dy. Needed 2D & 3D. ----
  {
    IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E3);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_OhmicEMF_E3", parthenon::DevExecSpace(), 0, nb - 1,
        kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(b);
          const Real dBy_dx =
              (pack(b, TE::F2, Bf(), k, j, i) - pack(b, TE::F2, Bf(), k, j, i - 1)) /
              c.Dxc<X1DIR>(k, j, i);
          const Real dBx_dy =
              (pack(b, TE::F1, Bf(), k, j, i) - pack(b, TE::F1, Bf(), k, j - 1, i)) /
              c.Dxc<X2DIR>(k, j, i);
          const Real j3 = dBy_dx - dBx_dy;
          const Real eta =
              0.25 * (CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k, j, i) +
                      CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k, j, i - 1) +
                      CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k, j - 1, i) +
                      CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k, j - 1, i - 1));
          pack.flux(b, TE::E3, Bf(), k, j, i) += eta * j3;
        });
  }

  if (three_d) {
    // ---- E1 edge (y-z corner, i,j-1/2,k-1/2). J1 = dBz/dy - dBy/dz. ----
    {
      IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E1);
      IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E1);
      IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E1);
      parthenon::par_for(
          DEFAULT_LOOP_PATTERN, "CT_OhmicEMF_E1", parthenon::DevExecSpace(), 0, nb - 1,
          kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
          KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
            const auto &c = pack.GetCoordinates(b);
            const Real dBz_dy =
                (pack(b, TE::F3, Bf(), k, j, i) - pack(b, TE::F3, Bf(), k, j - 1, i)) /
                c.Dxc<X2DIR>(k, j, i);
            const Real dBy_dz =
                (pack(b, TE::F2, Bf(), k, j, i) - pack(b, TE::F2, Bf(), k - 1, j, i)) /
                c.Dxc<X3DIR>(k, j, i);
            const Real j1 = dBz_dy - dBy_dz;
            const Real eta =
                0.25 * (CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k, j, i) +
                        CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k, j - 1, i) +
                        CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k - 1, j, i) +
                        CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k - 1, j - 1, i));
            pack.flux(b, TE::E1, Bf(), k, j, i) += eta * j1;
          });
    }
    // ---- E2 edge (z-x corner, i-1/2,j,k-1/2). J2 = dBx/dz - dBz/dx. ----
    {
      IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E2);
      IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E2);
      IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E2);
      parthenon::par_for(
          DEFAULT_LOOP_PATTERN, "CT_OhmicEMF_E2", parthenon::DevExecSpace(), 0, nb - 1,
          kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
          KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
            const auto &c = pack.GetCoordinates(b);
            const Real dBx_dz =
                (pack(b, TE::F1, Bf(), k, j, i) - pack(b, TE::F1, Bf(), k - 1, j, i)) /
                c.Dxc<X3DIR>(k, j, i);
            const Real dBz_dx =
                (pack(b, TE::F3, Bf(), k, j, i) - pack(b, TE::F3, Bf(), k, j, i - 1)) /
                c.Dxc<X1DIR>(k, j, i);
            const Real j2 = dBx_dz - dBz_dx;
            const Real eta =
                0.25 * (CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k, j, i) +
                        CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k, j, i - 1) +
                        CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k - 1, j, i) +
                        CellEtaO(use_cache, ohm_diff, eta_pack, prim, b, k - 1, j, i - 1));
            pack.flux(b, TE::E2, Bf(), k, j, i) += eta * j2;
          });
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Perpendicular-current EMF for ambipolar diffusion: E = eta (J - (J.bhat) bhat).
//! Local copy of PerpCurrentEMF (ambipolar.cpp) -- identical math; kept here so the CT
//! edge routine does not reach into that translation unit. Keep the two in sync.
KOKKOS_INLINE_FUNCTION
void ADPerpEMF(const Real eta, const Real j1, const Real j2, const Real j3, const Real b1,
               const Real b2, const Real b3, Real &e1, Real &e2, Real &e3) {
  const Real bsq = SQR(b1) + SQR(b2) + SQR(b3) + TINY_NUMBER;
  const Real jdotb = j1 * b1 + j2 * b2 + j3 * b3;
  e1 = eta * (j1 - jdotb * b1 / bsq);
  e2 = eta * (j2 - jdotb * b2 / bsq);
  e3 = eta * (j3 - jdotb * b3 / bsq);
}

//----------------------------------------------------------------------------------------
//! Hall EMF E = eta_h (J x B)/|B| + eta_floor J. Local copy of HallEMF (hall.cpp); keep in
//! sync. eta_floor is the optional Ohmic stabilizer (parabolic, real dissipation).
KOKKOS_INLINE_FUNCTION
void HallEMFLocal(const Real eta_h, const Real eta_floor, const Real bmag, const Real j1,
                  const Real j2, const Real j3, const Real b1, const Real b2, const Real b3,
                  Real &e1, Real &e2, Real &e3) {
  const Real inv_b = 1.0 / (bmag + TINY_NUMBER);
  const Real jxb1 = j2 * b3 - j3 * b2;
  const Real jxb2 = j3 * b1 - j1 * b3;
  const Real jxb3 = j1 * b2 - j2 * b1;
  e1 = eta_h * jxb1 * inv_b + eta_floor * j1;
  e2 = eta_h * jxb2 * inv_b + eta_floor * j2;
  e3 = eta_h * jxb3 * inv_b + eta_floor * j3;
}

//----------------------------------------------------------------------------------------
//! Face-centered current density J and face-averaged field B at the X1/X2/X3 face. These
//! mirror the current stencils + B averages of the GLM non-ideal kernels (resistivity/
//! ambipolar/hall.cpp) EXACTLY, so the CT non-ideal operators equal the validated GLM
//! operators by construction. The AD/Hall edge routines pair these with the term-specific
//! eta and EMF combine. Keep in sync with the .cpp kernels.
template <class Prim, class Coords>
KOKKOS_INLINE_FUNCTION void
FaceCurrentAndB_X1(const Prim &prim, const Coords &coords, const int ndim, const int b,
                   const int k, const int j, const int i, Real &j1, Real &j2, Real &j3,
                   Real &b1, Real &b2, Real &b3) {
  const auto d3B1 =
      ndim > 2 ? (0.5 * (prim(b, IB1, k + 1, j, i - 1) + prim(b, IB1, k + 1, j, i)) -
                  0.5 * (prim(b, IB1, k - 1, j, i - 1) + prim(b, IB1, k - 1, j, i))) /
                     (coords.template Xf<3, 1>(k + 1, j, i) -
                      coords.template Xf<3, 1>(k - 1, j, i))
               : 0.0;
  const auto d1B3 =
      (prim(b, IB3, k, j, i) - prim(b, IB3, k, j, i - 1)) / coords.template Dxc<1>(k, j, i);
  j2 = d3B1 - d1B3;
  const auto d1B2 =
      (prim(b, IB2, k, j, i) - prim(b, IB2, k, j, i - 1)) / coords.template Dxc<1>(k, j, i);
  const auto d2B1 =
      ndim > 1 ? (0.5 * (prim(b, IB1, k, j + 1, i - 1) + prim(b, IB1, k, j + 1, i)) -
                  0.5 * (prim(b, IB1, k, j - 1, i - 1) + prim(b, IB1, k, j - 1, i))) /
                     (coords.template Xf<2, 1>(k, j + 1, i) -
                      coords.template Xf<2, 1>(k, j - 1, i))
               : 0.0;
  j3 = d1B2 - d2B1;
  const auto d2B3 =
      ndim > 1 ? (0.5 * (prim(b, IB3, k, j + 1, i - 1) + prim(b, IB3, k, j + 1, i)) -
                  0.5 * (prim(b, IB3, k, j - 1, i - 1) + prim(b, IB3, k, j - 1, i))) /
                     (coords.template Xf<2, 1>(k, j + 1, i) -
                      coords.template Xf<2, 1>(k, j - 1, i))
               : 0.0;
  const auto d3B2 =
      ndim > 2 ? (0.5 * (prim(b, IB2, k + 1, j, i - 1) + prim(b, IB2, k + 1, j, i)) -
                  0.5 * (prim(b, IB2, k - 1, j, i - 1) + prim(b, IB2, k - 1, j, i))) /
                     (coords.template Xf<3, 1>(k + 1, j, i) -
                      coords.template Xf<3, 1>(k - 1, j, i))
               : 0.0;
  j1 = d2B3 - d3B2;
  b1 = 0.5 * (prim(b, IB1, k, j, i - 1) + prim(b, IB1, k, j, i));
  b2 = 0.5 * (prim(b, IB2, k, j, i - 1) + prim(b, IB2, k, j, i));
  b3 = 0.5 * (prim(b, IB3, k, j, i - 1) + prim(b, IB3, k, j, i));
}
template <class Prim, class Coords>
KOKKOS_INLINE_FUNCTION void
FaceCurrentAndB_X2(const Prim &prim, const Coords &coords, const int ndim, const int b,
                   const int k, const int j, const int i, Real &j1, Real &j2, Real &j3,
                   Real &b1, Real &b2, Real &b3) {
  const auto d1B2 = (0.5 * (prim(b, IB2, k, j - 1, i + 1) + prim(b, IB2, k, j, i + 1)) -
                     0.5 * (prim(b, IB2, k, j - 1, i - 1) + prim(b, IB2, k, j, i - 1))) /
                    (coords.template Xf<1, 2>(k, j, i + 1) -
                     coords.template Xf<1, 2>(k, j, i - 1));
  const auto d2B1 =
      (prim(b, IB1, k, j, i) - prim(b, IB1, k, j - 1, i)) / coords.template Dxc<2>(k, j, i);
  j3 = d1B2 - d2B1;
  const auto d2B3 =
      (prim(b, IB3, k, j, i) - prim(b, IB3, k, j - 1, i)) / coords.template Dxc<2>(k, j, i);
  const auto d3B2 =
      ndim > 2 ? (0.5 * (prim(b, IB2, k + 1, j - 1, i) + prim(b, IB2, k + 1, j, i)) -
                  0.5 * (prim(b, IB2, k - 1, j - 1, i) + prim(b, IB2, k - 1, j, i))) /
                     (coords.template Xf<3, 2>(k + 1, j, i) -
                      coords.template Xf<3, 2>(k - 1, j, i))
               : 0.0;
  j1 = d2B3 - d3B2;
  const auto d3B1 =
      ndim > 2 ? (0.5 * (prim(b, IB1, k + 1, j - 1, i) + prim(b, IB1, k + 1, j, i)) -
                  0.5 * (prim(b, IB1, k - 1, j - 1, i) + prim(b, IB1, k - 1, j, i))) /
                     (coords.template Xf<3, 2>(k + 1, j, i) -
                      coords.template Xf<3, 2>(k - 1, j, i))
               : 0.0;
  const auto d1B3 = (0.5 * (prim(b, IB3, k, j - 1, i + 1) + prim(b, IB3, k, j, i + 1)) -
                     0.5 * (prim(b, IB3, k, j - 1, i - 1) + prim(b, IB3, k, j, i - 1))) /
                    (coords.template Xf<1, 2>(k, j, i + 1) -
                     coords.template Xf<1, 2>(k, j, i - 1));
  j2 = d3B1 - d1B3;
  b1 = 0.5 * (prim(b, IB1, k, j - 1, i) + prim(b, IB1, k, j, i));
  b2 = 0.5 * (prim(b, IB2, k, j - 1, i) + prim(b, IB2, k, j, i));
  b3 = 0.5 * (prim(b, IB3, k, j - 1, i) + prim(b, IB3, k, j, i));
}
template <class Prim, class Coords>
KOKKOS_INLINE_FUNCTION void
FaceCurrentAndB_X3(const Prim &prim, const Coords &coords, const int ndim, const int b,
                   const int k, const int j, const int i, Real &j1, Real &j2, Real &j3,
                   Real &b1, Real &b2, Real &b3) {
  const auto d2B3 = (0.5 * (prim(b, IB3, k - 1, j + 1, i) + prim(b, IB3, k, j + 1, i)) -
                     0.5 * (prim(b, IB3, k - 1, j - 1, i) + prim(b, IB3, k, j - 1, i))) /
                    (coords.template Xf<2, 3>(k, j + 1, i) -
                     coords.template Xf<2, 3>(k, j - 1, i));
  const auto d3B2 =
      (prim(b, IB2, k, j, i) - prim(b, IB2, k - 1, j, i)) / coords.template Dxc<3>(k, j, i);
  j1 = d2B3 - d3B2;
  const auto d3B1 =
      (prim(b, IB1, k, j, i) - prim(b, IB1, k - 1, j, i)) / coords.template Dxc<3>(k, j, i);
  const auto d1B3 = (0.5 * (prim(b, IB3, k - 1, j, i + 1) + prim(b, IB3, k, j, i + 1)) -
                     0.5 * (prim(b, IB3, k - 1, j, i - 1) + prim(b, IB3, k, j, i - 1))) /
                    (coords.template Xf<1, 3>(k, j, i + 1) -
                     coords.template Xf<1, 3>(k, j, i - 1));
  j2 = d3B1 - d1B3;
  const auto d1B2 = (0.5 * (prim(b, IB2, k - 1, j, i + 1) + prim(b, IB2, k, j, i + 1)) -
                     0.5 * (prim(b, IB2, k - 1, j, i - 1) + prim(b, IB2, k, j, i - 1))) /
                    (coords.template Xf<1, 3>(k, j, i + 1) -
                     coords.template Xf<1, 3>(k, j, i - 1));
  const auto d2B1 = (0.5 * (prim(b, IB1, k - 1, j + 1, i) + prim(b, IB1, k, j + 1, i)) -
                     0.5 * (prim(b, IB1, k - 1, j - 1, i) + prim(b, IB1, k, j - 1, i))) /
                    (coords.template Xf<2, 3>(k, j + 1, i) -
                     coords.template Xf<2, 3>(k, j - 1, i));
  j3 = d1B2 - d2B1;
  b1 = 0.5 * (prim(b, IB1, k - 1, j, i) + prim(b, IB1, k, j, i));
  b2 = 0.5 * (prim(b, IB2, k - 1, j, i) + prim(b, IB2, k, j, i));
  b3 = 0.5 * (prim(b, IB3, k - 1, j, i) + prim(b, IB3, k, j, i));
}

//----------------------------------------------------------------------------------------
//! Ambipolar perp-current EMF (e1,e2,e3) at the X{1,2,3} face. Pairs FaceCurrentAndB with
//! the face-averaged eta_A (cached or from the face state) and the perp-current combine.
template <class Prim, class Coords, class Eta>
KOKKOS_INLINE_FUNCTION void
AmbiFaceEMF_X1(const Prim &prim, const Coords &coords, const AmbipolarDiffusivity &ad_diff,
               const Eta &eta_pack, const bool use_cache, const int ndim, const int b,
               const int k, const int j, const int i, Real &e1, Real &e2, Real &e3) {
  Real j1, j2, j3, b1, b2, b3;
  FaceCurrentAndB_X1(prim, coords, ndim, b, k, j, i, j1, j2, j3, b1, b2, b3);
  Real eta;
  if (use_cache) {
    eta = 0.5 * (eta_pack(b, NonidealEtaIdx::A, k, j, i - 1) +
                 eta_pack(b, NonidealEtaIdx::A, k, j, i));
  } else {
    const Real rho = 0.5 * (prim(b, IDN, k, j, i - 1) + prim(b, IDN, k, j, i));
    const Real prs = 0.5 * (prim(b, IPR, k, j, i - 1) + prim(b, IPR, k, j, i));
    const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
    const int i_xe = ad_diff.XeIndex();
    const Real xe =
        (i_xe >= 0) ? 0.5 * (prim(b, i_xe, k, j, i - 1) + prim(b, i_xe, k, j, i)) : -1.0;
    eta = ad_diff.Get(bmag, rho, prs / rho, xe);
  }
  ADPerpEMF(eta, j1, j2, j3, b1, b2, b3, e1, e2, e3);
}
template <class Prim, class Coords, class Eta>
KOKKOS_INLINE_FUNCTION void
AmbiFaceEMF_X2(const Prim &prim, const Coords &coords, const AmbipolarDiffusivity &ad_diff,
               const Eta &eta_pack, const bool use_cache, const int ndim, const int b,
               const int k, const int j, const int i, Real &e1, Real &e2, Real &e3) {
  Real j1, j2, j3, b1, b2, b3;
  FaceCurrentAndB_X2(prim, coords, ndim, b, k, j, i, j1, j2, j3, b1, b2, b3);
  Real eta;
  if (use_cache) {
    eta = 0.5 * (eta_pack(b, NonidealEtaIdx::A, k, j - 1, i) +
                 eta_pack(b, NonidealEtaIdx::A, k, j, i));
  } else {
    const Real rho = 0.5 * (prim(b, IDN, k, j - 1, i) + prim(b, IDN, k, j, i));
    const Real prs = 0.5 * (prim(b, IPR, k, j - 1, i) + prim(b, IPR, k, j, i));
    const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
    const int i_xe = ad_diff.XeIndex();
    const Real xe =
        (i_xe >= 0) ? 0.5 * (prim(b, i_xe, k, j - 1, i) + prim(b, i_xe, k, j, i)) : -1.0;
    eta = ad_diff.Get(bmag, rho, prs / rho, xe);
  }
  ADPerpEMF(eta, j1, j2, j3, b1, b2, b3, e1, e2, e3);
}
template <class Prim, class Coords, class Eta>
KOKKOS_INLINE_FUNCTION void
AmbiFaceEMF_X3(const Prim &prim, const Coords &coords, const AmbipolarDiffusivity &ad_diff,
               const Eta &eta_pack, const bool use_cache, const int ndim, const int b,
               const int k, const int j, const int i, Real &e1, Real &e2, Real &e3) {
  Real j1, j2, j3, b1, b2, b3;
  FaceCurrentAndB_X3(prim, coords, ndim, b, k, j, i, j1, j2, j3, b1, b2, b3);
  Real eta;
  if (use_cache) {
    eta = 0.5 * (eta_pack(b, NonidealEtaIdx::A, k - 1, j, i) +
                 eta_pack(b, NonidealEtaIdx::A, k, j, i));
  } else {
    const Real rho = 0.5 * (prim(b, IDN, k - 1, j, i) + prim(b, IDN, k, j, i));
    const Real prs = 0.5 * (prim(b, IPR, k - 1, j, i) + prim(b, IPR, k, j, i));
    const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
    const int i_xe = ad_diff.XeIndex();
    const Real xe =
        (i_xe >= 0) ? 0.5 * (prim(b, i_xe, k - 1, j, i) + prim(b, i_xe, k, j, i)) : -1.0;
    eta = ad_diff.Get(bmag, rho, prs / rho, xe);
  }
  ADPerpEMF(eta, j1, j2, j3, b1, b2, b3, e1, e2, e3);
}

//----------------------------------------------------------------------------------------
//! Hall EMF (e1,e2,e3) at the X{1,2,3} face. Pairs FaceCurrentAndB with the face-averaged
//! eta_H (signed) + the optional Ohmic floor, and the (J x B)/|B| combine. Mirrors
//! HallDiffFluxIsoFixed (hall.cpp) with eta_h_on=floor_on=true (the unsplit CT case).
template <class Prim, class Coords, class Eta>
KOKKOS_INLINE_FUNCTION void
HallFaceEMF_X1(const Prim &prim, const Coords &coords, const HallDiffusivity &hall_diff,
               const Real eta_floor, const Eta &eta_pack, const bool use_cache,
               const int ndim, const int b, const int k, const int j, const int i, Real &e1,
               Real &e2, Real &e3) {
  Real j1, j2, j3, b1, b2, b3;
  FaceCurrentAndB_X1(prim, coords, ndim, b, k, j, i, j1, j2, j3, b1, b2, b3);
  const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
  Real eta_h;
  if (use_cache) {
    eta_h = 0.5 * (eta_pack(b, NonidealEtaIdx::H, k, j, i - 1) +
                   eta_pack(b, NonidealEtaIdx::H, k, j, i));
  } else {
    const Real rho = 0.5 * (prim(b, IDN, k, j, i - 1) + prim(b, IDN, k, j, i));
    const Real prs = 0.5 * (prim(b, IPR, k, j, i - 1) + prim(b, IPR, k, j, i));
    const int i_xe = hall_diff.XeIndex();
    const Real xe =
        (i_xe >= 0) ? 0.5 * (prim(b, i_xe, k, j, i - 1) + prim(b, i_xe, k, j, i)) : -1.0;
    eta_h = hall_diff.Get(bmag, rho, prs / rho, xe);
  }
  HallEMFLocal(eta_h, eta_floor, bmag, j1, j2, j3, b1, b2, b3, e1, e2, e3);
}
template <class Prim, class Coords, class Eta>
KOKKOS_INLINE_FUNCTION void
HallFaceEMF_X2(const Prim &prim, const Coords &coords, const HallDiffusivity &hall_diff,
               const Real eta_floor, const Eta &eta_pack, const bool use_cache,
               const int ndim, const int b, const int k, const int j, const int i, Real &e1,
               Real &e2, Real &e3) {
  Real j1, j2, j3, b1, b2, b3;
  FaceCurrentAndB_X2(prim, coords, ndim, b, k, j, i, j1, j2, j3, b1, b2, b3);
  const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
  Real eta_h;
  if (use_cache) {
    eta_h = 0.5 * (eta_pack(b, NonidealEtaIdx::H, k, j - 1, i) +
                   eta_pack(b, NonidealEtaIdx::H, k, j, i));
  } else {
    const Real rho = 0.5 * (prim(b, IDN, k, j - 1, i) + prim(b, IDN, k, j, i));
    const Real prs = 0.5 * (prim(b, IPR, k, j - 1, i) + prim(b, IPR, k, j, i));
    const int i_xe = hall_diff.XeIndex();
    const Real xe =
        (i_xe >= 0) ? 0.5 * (prim(b, i_xe, k, j - 1, i) + prim(b, i_xe, k, j, i)) : -1.0;
    eta_h = hall_diff.Get(bmag, rho, prs / rho, xe);
  }
  HallEMFLocal(eta_h, eta_floor, bmag, j1, j2, j3, b1, b2, b3, e1, e2, e3);
}
template <class Prim, class Coords, class Eta>
KOKKOS_INLINE_FUNCTION void
HallFaceEMF_X3(const Prim &prim, const Coords &coords, const HallDiffusivity &hall_diff,
               const Real eta_floor, const Eta &eta_pack, const bool use_cache,
               const int ndim, const int b, const int k, const int j, const int i, Real &e1,
               Real &e2, Real &e3) {
  Real j1, j2, j3, b1, b2, b3;
  FaceCurrentAndB_X3(prim, coords, ndim, b, k, j, i, j1, j2, j3, b1, b2, b3);
  const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
  Real eta_h;
  if (use_cache) {
    eta_h = 0.5 * (eta_pack(b, NonidealEtaIdx::H, k - 1, j, i) +
                   eta_pack(b, NonidealEtaIdx::H, k, j, i));
  } else {
    const Real rho = 0.5 * (prim(b, IDN, k - 1, j, i) + prim(b, IDN, k, j, i));
    const Real prs = 0.5 * (prim(b, IPR, k - 1, j, i) + prim(b, IPR, k, j, i));
    const int i_xe = hall_diff.XeIndex();
    const Real xe =
        (i_xe >= 0) ? 0.5 * (prim(b, i_xe, k - 1, j, i) + prim(b, i_xe, k, j, i)) : -1.0;
    eta_h = hall_diff.Get(bmag, rho, prs / rho, xe);
  }
  HallEMFLocal(eta_h, eta_floor, bmag, j1, j2, j3, b1, b2, b3, e1, e2, e3);
}

//----------------------------------------------------------------------------------------
//! Non-ideal (ambipolar) edge EMF for CT. See ct.hpp. The perp-current EMF is evaluated at
//! cell faces (AmbiFaceEMF_X{1,2,3}, same stencils as the GLM path) and the relevant
//! component is arithmetic-averaged from the four faces bounding each edge -- the exact
//! four-face index pattern GS05 uses for the ideal base EMF, so div B stays at round-off.
TaskStatus CT_AddAmbipolarEMF(MeshData<Real> *md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto hydro_pkg = pmb->packages.Get("Hydro");
  // No-op unless ambipolar diffusion is active (CT diffusion is unsplit-only).
  if (hydro_pkg->Param<Ambipolar>("ambipolar") == Ambipolar::none) {
    return TaskStatus::complete;
  }
  const int ndim = md->GetMeshPointer()->ndim;
  const bool three_d = ndim > 2;
  const auto &ad_diff = hydro_pkg->Param<AmbipolarDiffusivity>("ad_diff");
  const bool use_cache = hydro_pkg->Param<bool>("nonideal_eta_cache");

  auto prim = md->PackVariables(std::vector<std::string>{"prim"});
  const auto eta_pack =
      md->PackVariables(std::vector<std::string>{use_cache ? "nonideal_eta" : "prim"});
  static auto desc =
      parthenon::MakePackDescriptor<Bf>(md, {}, {parthenon::PDOpt::WithFluxes});
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();

  // ---- E3 edge (x-y corner, i-1/2,j-1/2,k): e3 from x-faces (rows j,j-1) + y-faces
  //      (cols i,i-1). Needed in 2D & 3D. ----
  {
    IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E3);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_AmbiEMF_E3", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(b);
          Real e1, e2, e3, acc = 0.0;
          AmbiFaceEMF_X1(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k, j, i, e1, e2,
                         e3);
          acc += e3;
          AmbiFaceEMF_X1(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k, j - 1, i, e1,
                         e2, e3);
          acc += e3;
          AmbiFaceEMF_X2(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k, j, i, e1, e2,
                         e3);
          acc += e3;
          AmbiFaceEMF_X2(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k, j, i - 1, e1,
                         e2, e3);
          acc += e3;
          pack.flux(b, TE::E3, Bf(), k, j, i) += 0.25 * acc;
        });
  }

  if (three_d) {
    // ---- E1 edge (y-z corner, i,j-1/2,k-1/2): e1 from y-faces (rows k,k-1) + z-faces
    //      (cols j,j-1). ----
    {
      IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E1);
      IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E1);
      IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E1);
      parthenon::par_for(
          DEFAULT_LOOP_PATTERN, "CT_AmbiEMF_E1", parthenon::DevExecSpace(), 0, nb - 1,
          kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
          KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
            const auto &c = pack.GetCoordinates(b);
            Real e1, e2, e3, acc = 0.0;
            AmbiFaceEMF_X2(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k, j, i, e1, e2,
                           e3);
            acc += e1;
            AmbiFaceEMF_X2(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k - 1, j, i, e1,
                           e2, e3);
            acc += e1;
            AmbiFaceEMF_X3(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k, j, i, e1, e2,
                           e3);
            acc += e1;
            AmbiFaceEMF_X3(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k, j - 1, i, e1,
                           e2, e3);
            acc += e1;
            pack.flux(b, TE::E1, Bf(), k, j, i) += 0.25 * acc;
          });
    }
    // ---- E2 edge (z-x corner, i-1/2,j,k-1/2): e2 from z-faces (cols i,i-1) + x-faces
    //      (rows k,k-1). ----
    {
      IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E2);
      IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E2);
      IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E2);
      parthenon::par_for(
          DEFAULT_LOOP_PATTERN, "CT_AmbiEMF_E2", parthenon::DevExecSpace(), 0, nb - 1,
          kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
          KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
            const auto &c = pack.GetCoordinates(b);
            Real e1, e2, e3, acc = 0.0;
            AmbiFaceEMF_X3(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k, j, i, e1, e2,
                           e3);
            acc += e2;
            AmbiFaceEMF_X3(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k, j, i - 1, e1,
                           e2, e3);
            acc += e2;
            AmbiFaceEMF_X1(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k, j, i, e1, e2,
                           e3);
            acc += e2;
            AmbiFaceEMF_X1(prim, c, ad_diff, eta_pack, use_cache, ndim, b, k - 1, j, i, e1,
                           e2, e3);
            acc += e2;
            pack.flux(b, TE::E2, Bf(), k, j, i) += 0.25 * acc;
          });
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Non-ideal (Hall) edge EMF for CT. Same four-face-average construction as
//! CT_AddAmbipolarEMF, but with the Hall EMF E_H = eta_H (J x B)/|B| (+ optional Ohmic
//! floor eta_O J). Hall is dispersive and unsplit-only under CT (RKL2+CT forbidden), so the
//! whistler part and the floor are both applied here (eta_h_on = floor_on = true). See
//! ct.hpp. The matching cons.flux(IBn) deposit in HallDiffFluxIsoFixed is gated off under
//! CT; its cons.flux(IEN) Poynting term stays on the FV energy flux.
TaskStatus CT_AddHallEMF(MeshData<Real> *md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto hydro_pkg = pmb->packages.Get("Hydro");
  if (hydro_pkg->Param<Hall>("hall") == Hall::none) {
    return TaskStatus::complete;
  }
  const int ndim = md->GetMeshPointer()->ndim;
  const bool three_d = ndim > 2;
  const auto &hall_diff = hydro_pkg->Param<HallDiffusivity>("hall_diff");
  const Real eta_floor = hall_diff.GetOhmicFloor();
  const bool use_cache = hydro_pkg->Param<bool>("nonideal_eta_cache");

  auto prim = md->PackVariables(std::vector<std::string>{"prim"});
  const auto eta_pack =
      md->PackVariables(std::vector<std::string>{use_cache ? "nonideal_eta" : "prim"});
  static auto desc =
      parthenon::MakePackDescriptor<Bf>(md, {}, {parthenon::PDOpt::WithFluxes});
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();

  // ---- E3 edge: e3 from x-faces (rows j,j-1) + y-faces (cols i,i-1). 2D & 3D. ----
  {
    IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E3);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_HallEMF_E3", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(b);
          Real e1, e2, e3, acc = 0.0;
          HallFaceEMF_X1(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b, k, j,
                         i, e1, e2, e3);
          acc += e3;
          HallFaceEMF_X1(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b, k,
                         j - 1, i, e1, e2, e3);
          acc += e3;
          HallFaceEMF_X2(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b, k, j,
                         i, e1, e2, e3);
          acc += e3;
          HallFaceEMF_X2(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b, k, j,
                         i - 1, e1, e2, e3);
          acc += e3;
          pack.flux(b, TE::E3, Bf(), k, j, i) += 0.25 * acc;
        });
  }

  if (three_d) {
    // ---- E1 edge: e1 from y-faces (rows k,k-1) + z-faces (cols j,j-1). ----
    {
      IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E1);
      IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E1);
      IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E1);
      parthenon::par_for(
          DEFAULT_LOOP_PATTERN, "CT_HallEMF_E1", parthenon::DevExecSpace(), 0, nb - 1,
          kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
          KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
            const auto &c = pack.GetCoordinates(b);
            Real e1, e2, e3, acc = 0.0;
            HallFaceEMF_X2(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b, k,
                           j, i, e1, e2, e3);
            acc += e1;
            HallFaceEMF_X2(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b,
                           k - 1, j, i, e1, e2, e3);
            acc += e1;
            HallFaceEMF_X3(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b, k,
                           j, i, e1, e2, e3);
            acc += e1;
            HallFaceEMF_X3(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b, k,
                           j - 1, i, e1, e2, e3);
            acc += e1;
            pack.flux(b, TE::E1, Bf(), k, j, i) += 0.25 * acc;
          });
    }
    // ---- E2 edge: e2 from z-faces (cols i,i-1) + x-faces (rows k,k-1). ----
    {
      IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E2);
      IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E2);
      IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E2);
      parthenon::par_for(
          DEFAULT_LOOP_PATTERN, "CT_HallEMF_E2", parthenon::DevExecSpace(), 0, nb - 1,
          kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
          KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
            const auto &c = pack.GetCoordinates(b);
            Real e1, e2, e3, acc = 0.0;
            HallFaceEMF_X3(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b, k,
                           j, i, e1, e2, e3);
            acc += e2;
            HallFaceEMF_X3(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b, k,
                           j, i - 1, e1, e2, e3);
            acc += e2;
            HallFaceEMF_X1(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b, k,
                           j, i, e1, e2, e3);
            acc += e2;
            HallFaceEMF_X1(prim, c, hall_diff, eta_floor, eta_pack, use_cache, ndim, b,
                           k - 1, j, i, e1, e2, e3);
            acc += e2;
            pack.flux(b, TE::E2, Bf(), k, j, i) += 0.25 * acc;
          });
    }
  }
  return TaskStatus::complete;
}

} // namespace CT
} // namespace Hydro
