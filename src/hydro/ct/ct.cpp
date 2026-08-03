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
//! Discrete curl of the edge EMF (stored in the Bf.flux slots of `emf`) for one face,
//! i.e. the M-operator value (-curl E) used by the CT RKL2 super-time-stepping. Mirrors the
//! per-face curl of CT_UpdateBf exactly. `comp` selects F1/F2/F3.
template <int comp, class Pack, class Coords>
KOKKOS_INLINE_FUNCTION Real CurlEMF(const Pack &emf, const Coords &coords, const int b,
                                    const int k, const int j, const int i,
                                    const bool three_d) {
  if (comp == 1) { // F1 (B_x): -(dE3/dy) + (dE2/dz)
    Real curl = -(emf.flux(b, TE::E3, Bf(), k, j + 1, i) -
                  emf.flux(b, TE::E3, Bf(), k, j, i)) /
                coords.template Dxc<X2DIR>(k, j, i);
    if (three_d) {
      curl += (emf.flux(b, TE::E2, Bf(), k + 1, j, i) -
               emf.flux(b, TE::E2, Bf(), k, j, i)) /
              coords.template Dxc<X3DIR>(k, j, i);
    }
    return curl;
  } else if (comp == 2) { // F2 (B_y): (dE3/dx) - (dE1/dz)
    Real curl = (emf.flux(b, TE::E3, Bf(), k, j, i + 1) -
                 emf.flux(b, TE::E3, Bf(), k, j, i)) /
                coords.template Dxc<X1DIR>(k, j, i);
    if (three_d) {
      curl -= (emf.flux(b, TE::E1, Bf(), k + 1, j, i) -
               emf.flux(b, TE::E1, Bf(), k, j, i)) /
              coords.template Dxc<X3DIR>(k, j, i);
    }
    return curl;
  } else { // F3 (B_z): -(dE2/dx) + (dE1/dy)  [3D only]
    return -(emf.flux(b, TE::E2, Bf(), k, j, i + 1) -
             emf.flux(b, TE::E2, Bf(), k, j, i)) /
               coords.template Dxc<X1DIR>(k, j, i) +
           (emf.flux(b, TE::E1, Bf(), k, j + 1, i) -
            emf.flux(b, TE::E1, Bf(), k, j, i)) /
               coords.template Dxc<X2DIR>(k, j, i);
  }
}

//----------------------------------------------------------------------------------------
//! Zero the edge-EMF flux slots so the += accumulators build the diffusive-only EMF.
TaskStatus CT_ZeroEMF(MeshData<Real> *md) {
  const int ndim = md->GetMeshPointer()->ndim;
  static auto desc =
      parthenon::MakePackDescriptor<Bf>(md, {}, {parthenon::PDOpt::WithFluxes});
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();
  auto zero = [&](TE te, const char *name) {
    IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, te);
    IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, te);
    IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, te);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, name, parthenon::DevExecSpace(), 0, nb - 1, kb.s, kb.e, jb.s,
        jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          pack.flux(b, te, Bf(), k, j, i) = 0.0;
        });
  };
  zero(TE::E3, "CT_ZeroEMF_E3");
  if (ndim > 2) {
    zero(TE::E1, "CT_ZeroEMF_E1");
    zero(TE::E2, "CT_ZeroEMF_E2");
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! MY0 = M(Y0) = -curl(E_diff): write md_out's face field to the curl of the edge EMF in
//! md_emf's Bf.flux. (Both containers carry the Bf variable; md_emf must have valid fluxes.)
TaskStatus CT_CurlEMFToBf(MeshData<Real> *md_emf, MeshData<Real> *md_out) {
  const int ndim = md_emf->GetMeshPointer()->ndim;
  const bool three_d = ndim > 2;
  static auto desc_e =
      parthenon::MakePackDescriptor<Bf>(md_emf, {}, {parthenon::PDOpt::WithFluxes});
  static auto desc_o = parthenon::MakePackDescriptor<Bf>(md_out);
  auto emf = desc_e.GetPack(md_emf);
  auto out = desc_o.GetPack(md_out);
  const int nb = emf.GetNBlocks();
  // Explicit per-face blocks (not a generic lambda: nvcc forbids an extended
  // __host__ __device__ lambda inside a generic lambda).
  {
    IndexRange ib = md_emf->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::F1);
    IndexRange jb = md_emf->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::F1);
    IndexRange kb = md_emf->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::F1);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_CurlEMF_F1", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = emf.GetCoordinates(b);
          out(b, TE::F1, Bf(), k, j, i) = CurlEMF<1>(emf, c, b, k, j, i, three_d);
        });
  }
  {
    IndexRange ib = md_emf->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::F2);
    IndexRange jb = md_emf->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::F2);
    IndexRange kb = md_emf->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::F2);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_CurlEMF_F2", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = emf.GetCoordinates(b);
          out(b, TE::F2, Bf(), k, j, i) = CurlEMF<2>(emf, c, b, k, j, i, three_d);
        });
  }
  if (three_d) {
    IndexRange ib = md_emf->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::F3);
    IndexRange jb = md_emf->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::F3);
    IndexRange kb = md_emf->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::F3);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_CurlEMF_F3", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = emf.GetCoordinates(b);
          out(b, TE::F3, Bf(), k, j, i) = CurlEMF<3>(emf, c, b, k, j, i, three_d);
        });
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! RKL2 first sub-step on the face field: base.Bf = Y0.Bf + mu_tilde_1*tau*MY0.Bf and
//! Yjm2.Bf = Y0.Bf. Mirrors RKL2StepFirst but for Bf (F1/F2/F3).
TaskStatus CT_RKL2FirstBf(MeshData<Real> *md_Y0, MeshData<Real> *md_base,
                          MeshData<Real> *md_Yjm2, MeshData<Real> *md_MY0, const int s_rkl,
                          const Real tau) {
  const Real mu_tilde_1 = 4. / 3. /
                          (static_cast<Real>(s_rkl) * static_cast<Real>(s_rkl) +
                           static_cast<Real>(s_rkl) - 2.);
  const int ndim = md_base->GetMeshPointer()->ndim;
  static auto desc = parthenon::MakePackDescriptor<Bf>(md_base);
  auto Y0 = desc.GetPack(md_Y0);
  auto base = desc.GetPack(md_base);
  auto Yjm2 = desc.GetPack(md_Yjm2);
  auto MY0 = desc.GetPack(md_MY0);
  const int nb = base.GetNBlocks();
  auto run = [&](TE te, const char *name) {
    IndexRange ib = md_base->GetBoundsI(CellLevel::same, IndexDomain::interior, te);
    IndexRange jb = md_base->GetBoundsJ(CellLevel::same, IndexDomain::interior, te);
    IndexRange kb = md_base->GetBoundsK(CellLevel::same, IndexDomain::interior, te);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, name, parthenon::DevExecSpace(), 0, nb - 1, kb.s, kb.e, jb.s,
        jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const Real y0 = Y0(b, te, Bf(), k, j, i);
          base(b, te, Bf(), k, j, i) = y0 + mu_tilde_1 * tau * MY0(b, te, Bf(), k, j, i);
          Yjm2(b, te, Bf(), k, j, i) = y0;
        });
  };
  run(TE::F1, "CT_RKL2First_F1");
  run(TE::F2, "CT_RKL2First_F2");
  if (ndim > 2) run(TE::F3, "CT_RKL2First_F3");
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! RKL2 sub-step j>=2 on the face field. MYjm1 = -curl(edge EMF in base.Bf.flux) is formed
//! inline; then base.Bf is advanced by the RKL2 recurrence with the Yjm2 shuffle. Mirrors
//! RKL2StepOther but for Bf, reading the edge EMF from base's own flux slots.
TaskStatus CT_RKL2OtherBf(MeshData<Real> *md_Y0, MeshData<Real> *md_base,
                          MeshData<Real> *md_Yjm2, MeshData<Real> *md_MY0, const Real mu_j,
                          const Real nu_j, const Real mu_tilde_j, const Real gamma_tilde_j,
                          const Real tau) {
  const int ndim = md_base->GetMeshPointer()->ndim;
  const bool three_d = ndim > 2;
  static auto desc_f =
      parthenon::MakePackDescriptor<Bf>(md_base, {}, {parthenon::PDOpt::WithFluxes});
  static auto desc = parthenon::MakePackDescriptor<Bf>(md_base);
  auto base = desc_f.GetPack(md_base); // has Bf + edge-EMF fluxes
  auto Y0 = desc.GetPack(md_Y0);
  auto Yjm2 = desc.GetPack(md_Yjm2);
  auto MY0 = desc.GetPack(md_MY0);
  const int nb = base.GetNBlocks();
  // Explicit per-face blocks (nvcc forbids an extended lambda inside a generic lambda).
  {
    IndexRange ib = md_base->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::F1);
    IndexRange jb = md_base->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::F1);
    IndexRange kb = md_base->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::F1);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_RKL2Other_F1", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = base.GetCoordinates(b);
          const Real MYjm1 = CurlEMF<1>(base, c, b, k, j, i, three_d);
          const Real yjm1 = base(b, TE::F1, Bf(), k, j, i);
          const Real yj = mu_j * yjm1 + nu_j * Yjm2(b, TE::F1, Bf(), k, j, i) +
                          (1.0 - mu_j - nu_j) * Y0(b, TE::F1, Bf(), k, j, i) +
                          mu_tilde_j * tau * MYjm1 +
                          gamma_tilde_j * tau * MY0(b, TE::F1, Bf(), k, j, i);
          Yjm2(b, TE::F1, Bf(), k, j, i) = yjm1;
          base(b, TE::F1, Bf(), k, j, i) = yj;
        });
  }
  {
    IndexRange ib = md_base->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::F2);
    IndexRange jb = md_base->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::F2);
    IndexRange kb = md_base->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::F2);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_RKL2Other_F2", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = base.GetCoordinates(b);
          const Real MYjm1 = CurlEMF<2>(base, c, b, k, j, i, three_d);
          const Real yjm1 = base(b, TE::F2, Bf(), k, j, i);
          const Real yj = mu_j * yjm1 + nu_j * Yjm2(b, TE::F2, Bf(), k, j, i) +
                          (1.0 - mu_j - nu_j) * Y0(b, TE::F2, Bf(), k, j, i) +
                          mu_tilde_j * tau * MYjm1 +
                          gamma_tilde_j * tau * MY0(b, TE::F2, Bf(), k, j, i);
          Yjm2(b, TE::F2, Bf(), k, j, i) = yjm1;
          base(b, TE::F2, Bf(), k, j, i) = yj;
        });
  }
  if (three_d) {
    IndexRange ib = md_base->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::F3);
    IndexRange jb = md_base->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::F3);
    IndexRange kb = md_base->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::F3);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_RKL2Other_F3", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = base.GetCoordinates(b);
          const Real MYjm1 = CurlEMF<3>(base, c, b, k, j, i, three_d);
          const Real yjm1 = base(b, TE::F3, Bf(), k, j, i);
          const Real yj = mu_j * yjm1 + nu_j * Yjm2(b, TE::F3, Bf(), k, j, i) +
                          (1.0 - mu_j - nu_j) * Y0(b, TE::F3, Bf(), k, j, i) +
                          mu_tilde_j * tau * MYjm1 +
                          gamma_tilde_j * tau * MY0(b, TE::F3, Bf(), k, j, i);
          Yjm2(b, TE::F3, Bf(), k, j, i) = yjm1;
          base(b, TE::F3, Bf(), k, j, i) = yj;
        });
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Project the face field onto cell-centered IB1..IB3 in cons.
//----------------------------------------------------------------------------------------
//! Diffusive Poynting energy flux built from the SAME edge EMF that drives the CT induction.
//!
//! WHY THIS EXISTS. Under CT the diffusive induction comes from the edge EMF in Bf.flux(E*),
//! but the matching energy flux was still being deposited by the face-based kernels in
//! ambipolar.cpp / resistivity.cpp, which build their own EMF from a DIFFERENT stencil (face-
//! averaged cell-centered B and a face-centered current). The gas internal energy is recovered
//! as e = E - KE - ME, so any disagreement between the energy flux and the field update lands
//! directly in e. For Ohmic the two stencils are close (edge EMF = eta*J from the tight 1dx
//! curl of Bf) and the error is harmless; for ambipolar they are not, and since the 2026-07-28
//! direct-edge AD stencil (AmbiEdgeEMF_E*, tight in-plane + wide 2-cell along-edge) they are
//! badly mismatched. Measured consequence: with the guard disabled, CT+RKL2 Orszag-Tang runs
//! CLEAN with Ohmic only (cycle 300, matches GLM) but aborts on negative pressure at cycle
//! ~100 with ambipolar. That isolates the defect to the AD energy/induction stencil pair.
//!
//! WHAT IT DOES. S = E x B, deposited on faces from the edge EMFs bounding each face:
//!   S_1 = E2 B3 - E3 B2,  S_2 = E3 B1 - E1 B3,  S_3 = E1 B2 - E2 B1
//! with each transverse EMF averaged from the two edges that bound the face, and the face-
//! tangential B taken as the two-cell average of the cell-centered field -- i.e. the identical
//! b1/b2/b3 the face kernels use, so ONLY the EMF source changes. Sign convention matches the
//! existing deposits (ambipolar.cpp: flux(X1DIR,IEN) += e2*b3 - e3*b2).
//!
//! Bf.flux must hold the DIFFUSIVE-ONLY edge EMF, so this is wired into the RKL2/STS path only
//! (CT_ZeroEMF -> Ohmic -> ambipolar). In the unsplit path Bf.flux carries ideal+diffusive and
//! the ideal Poynting is already in the HLLD energy flux; that path is left on the face-based
//! deposits (and is not affected by the defect -- CT+AD-unsplit runs clean). Must run before
//! the flux-correction round so the coarse-fine reflux restricts this deposit too.
TaskStatus CT_AddDiffusivePoynting(MeshData<Real> *md) {
  const int ndim = md->GetMeshPointer()->ndim;
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  if (!hydro_pkg->Param<bool>("ct_edge_poynting")) {
    return TaskStatus::complete;
  }
  const std::vector<std::string> cons_names{"cons"};
  auto cons_pack = md->PackVariablesAndFluxes(cons_names, cons_names);
  const auto &prim = md->PackVariables(std::vector<std::string>{"prim"});
  static auto desc =
      parthenon::MakePackDescriptor<Bf>(md, {}, {parthenon::PDOpt::WithFluxes});
  auto pb = desc.GetPack(md);
  const int nb = pb.GetNBlocks();
  const bool three_d = ndim > 2;
  const bool two_d = ndim > 1;

  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);

  // ---- X1 faces: S_1 = E2 B3 - E3 B2 ----
  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "CT_Poynting_X1", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
      kb.e, jb.s, jb.e, ib.s, ib.e + 1,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const Real e2 = three_d ? 0.5 * (pb.flux(b, TE::E2, Bf(), k, j, i) +
                                         pb.flux(b, TE::E2, Bf(), k + 1, j, i))
                                : pb.flux(b, TE::E2, Bf(), k, j, i);
        const Real e3 = two_d ? 0.5 * (pb.flux(b, TE::E3, Bf(), k, j, i) +
                                       pb.flux(b, TE::E3, Bf(), k, j + 1, i))
                              : pb.flux(b, TE::E3, Bf(), k, j, i);
        const Real b2 = 0.5 * (prim(b, IB2, k, j, i - 1) + prim(b, IB2, k, j, i));
        const Real b3 = 0.5 * (prim(b, IB3, k, j, i - 1) + prim(b, IB3, k, j, i));
        cons_pack(b).flux(X1DIR, IEN, k, j, i) += e2 * b3 - e3 * b2;
      });

  // ---- X2 faces: S_2 = E3 B1 - E1 B3 ----
  if (two_d) {
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_Poynting_X2", parthenon::DevExecSpace(), 0, nb - 1,
        kb.s, kb.e, jb.s, jb.e + 1, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const Real e3 = 0.5 * (pb.flux(b, TE::E3, Bf(), k, j, i) +
                                 pb.flux(b, TE::E3, Bf(), k, j, i + 1));
          const Real e1 = three_d ? 0.5 * (pb.flux(b, TE::E1, Bf(), k, j, i) +
                                           pb.flux(b, TE::E1, Bf(), k + 1, j, i))
                                  : pb.flux(b, TE::E1, Bf(), k, j, i);
          const Real b1 = 0.5 * (prim(b, IB1, k, j - 1, i) + prim(b, IB1, k, j, i));
          const Real b3 = 0.5 * (prim(b, IB3, k, j - 1, i) + prim(b, IB3, k, j, i));
          cons_pack(b).flux(X2DIR, IEN, k, j, i) += e3 * b1 - e1 * b3;
        });
  }

  // ---- X3 faces: S_3 = E1 B2 - E2 B1 ----
  if (three_d) {
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_Poynting_X3", parthenon::DevExecSpace(), 0, nb - 1,
        kb.s, kb.e + 1, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const Real e1 = 0.5 * (pb.flux(b, TE::E1, Bf(), k, j, i) +
                                 pb.flux(b, TE::E1, Bf(), k, j + 1, i));
          const Real e2 = 0.5 * (pb.flux(b, TE::E2, Bf(), k, j, i) +
                                 pb.flux(b, TE::E2, Bf(), k, j, i + 1));
          const Real b1 = 0.5 * (prim(b, IB1, k - 1, j, i) + prim(b, IB1, k, j, i));
          const Real b2 = 0.5 * (prim(b, IB2, k - 1, j, i) + prim(b, IB2, k, j, i));
          cons_pack(b).flux(X3DIR, IEN, k, j, i) += e1 * b2 - e2 * b1;
        });
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
TaskStatus CT_ProjectBfToCC(MeshData<Real> *md) {
  const int ndim = md->GetMeshPointer()->ndim;
  auto cons = md->PackVariables(std::vector<std::string>{"cons"});
  static auto desc = parthenon::MakePackDescriptor<Bf>(md);
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();
  const int koff = (ndim > 2) ? 1 : 0;
  const int joff = (ndim > 1) ? 1 : 0;

  // Internal-energy guard on the magnetic-energy replacement (see the long note below).
  // 0 disables it -> bit-identical to the unguarded projection.
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  const Real eint_guard = hydro_pkg->Param<Real>("ct_eint_guard_frac");
  const bool guard = eint_guard > 0.0;

  // DIAGNOSTIC (hydro/ct_proj_diag, default false): record per cell the RELATIVE internal-
  // energy change the face->cell projection imposes,
  //     dEint_rel = (eint_post - eint_pre)/eint_pre = -(me_post - me_pre)/eint_pre,
  // i.e. how much energy the projection silently moves into (positive) or out of (negative)
  // the gas, because e is recovered as E - KE - ME while E and ME are advanced by different
  // integrators. This is the DIRECT measurement of the mechanism that seven ablations have
  // only bounded: if |dEint_rel| is small inside the evacuated region, the projection is
  // exonerated and the cause lies in the source-term coupling; if it is O(1) there, the
  // mechanism is identified with a number. Recorded BEFORE the guard acts, so it measures the
  // raw discrepancy, not the guarded one. Field "ct.dEint" is Derived (not in restarts).
  const bool proj_diag = hydro_pkg->Param<bool>("ct_proj_diag");
  // Energy-neutral projection (see the long note in the kernel). Default ON under CT.
  const bool energy_neutral = hydro_pkg->Param<bool>("ct_energy_neutral_projection");
  auto diag = proj_diag ? md->PackVariables(std::vector<std::string>{"ct.dEint"})
                        : md->PackVariables(std::vector<std::string>{"cons"});

  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "CT_ProjectBfToCC", parthenon::DevExecSpace(), 0, nb - 1,
      kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        // Magnetic energy the cons/energy update was built against, BEFORE we overwrite it.
        const Real me_pre =
            0.5 * (SQR(cons(b, IB1, k, j, i)) + SQR(cons(b, IB2, k, j, i)) +
                   SQR(cons(b, IB3, k, j, i)));

        const Real b1 = 0.5 * (pack(b, TE::F1, Bf(), k, j, i) +
                               pack(b, TE::F1, Bf(), k, j, i + 1));
        const Real b2 = 0.5 * (pack(b, TE::F2, Bf(), k, j, i) +
                               pack(b, TE::F2, Bf(), k, j + joff, i));
        const Real b3 = 0.5 * (pack(b, TE::F3, Bf(), k, j, i) +
                               pack(b, TE::F3, Bf(), k + koff, j, i));
        cons(b, IB1, k, j, i) = b1;
        cons(b, IB2, k, j, i) = b2;
        cons(b, IB3, k, j, i) = b3;

        // The gas internal energy is recovered downstream as e = E - KE - ME. But E and ME
        // are advanced by two DIFFERENT integrators: E by the conservative flux divergence
        // (whose Poynting/energy flux is consistent with the cell-centered, Riemann-solved
        // field), ME by the CT curl of the edge EMF on the face field -- and in the RKL2
        // path by a second, independent super-time-stepping recurrence (RKL2StepFirst on
        // cons vs CT_RKL2FirstBf on Bf). Nothing makes their magnetic-energy bookkeeping
        // agree, so the whole discrepancy (me_pre - me_post) is silently dumped into the
        // internal energy. At a strong shock, or across an STS super-step, that transfer is
        // O(1) rather than truncation-small and of arbitrary SIGN; when it is negative and
        // exceeds e, the cell goes to negative pressure. Downstream the pressure floor then
        // "repairs" it by injecting an unrelated amount of energy, so the thermodynamics is
        // lost. (Cheap reproducer: runs/ct_tests/orszag_tang_ad_ct.in -- CT+AD+RKL2 aborts
        // on negative pressure at cycle ~100, while ideal-CT and CT+AD-unsplit are clean.)
        //
        // A magnetic-energy re-synchronization is physically only allowed to be a small
        // NON-NEGATIVE numerical dissipation. We therefore guard one-sidedly: the projection
        // may heat the gas (dissipation, stabilizing) but may not cool it below a fraction
        // eint_guard of the internal energy the conservative update produced. When the guard
        // does not fire, cons(IEN) is untouched and the scheme is exactly conservative and
        // bit-identical. When it does, we correct E rather than let the pressure floor fire:
        // a strictly more local and more accurate repair, using this cell's own pre-
        // projection internal energy instead of a global floor constant.
        if (guard || proj_diag || energy_neutral) {
          const Real rho = cons(b, IDN, k, j, i);
          const Real e_kin =
              (rho > 0.0) ? 0.5 *
                                (SQR(cons(b, IM1, k, j, i)) + SQR(cons(b, IM2, k, j, i)) +
                                 SQR(cons(b, IM3, k, j, i))) /
                                rho
                          : 0.0;
          const Real etot = cons(b, IEN, k, j, i);
          const Real eint_pre = etot - e_kin - me_pre;
          const Real me_post = 0.5 * (SQR(b1) + SQR(b2) + SQR(b3));
          const Real eint_post = etot - e_kin - me_post;

          // RAW relative internal-energy change the projection would impose. Recorded even
          // when the neutralization below cancels it, because it measures a property of the
          // SCHEME (the size of the E-vs-ME bookkeeping inconsistency), which stays
          // diagnostically useful after the fix.
          if (proj_diag) {
            diag(b, 0, k, j, i) =
                (eint_pre > 0.0) ? (eint_post - eint_pre) / eint_pre : 0.0;
          }

          if (energy_neutral) {
            // ENERGY-NEUTRAL PROJECTION. Replacing the cell-centered B by the projection of
            // the face field is a change of magnetic REPRESENTATION, not a physical process:
            // it must not add or remove gas heat. Carrying the magnetic-energy change into
            // the total energy holds e = E - KE - ME exactly fixed across the projection.
            //
            // WHY THIS IS THE RIGHT PLACE TO PUT THE DISCREPANCY. The same absolute quantity
            // D = me_post - me_pre must land somewhere. Currently it lands in e, where the
            // relative damage is D/e -- and in a magnetized cell e is the SMALLEST of the
            // three energies, so the error is amplified by ME/IE (measured up to ~150 in the
            // flagship). Sending it to E instead gives D/E, de-amplified because E >= ME >> D.
            // Measured on runs/flagship_integration/fc128diag at t=1.101883:
            //   cells with ME/IE>1: median relative error 1.331e-03 (into e) vs 2.976e-05
            //                       (into E)  ->  45x better
            //   cells with ME/IE>3: 7.916e-03 vs 7.682e-05          -> 103x better
            // The cost is exact total-energy conservation, and it is tiny: the net injection
            // is -7.6e-08 % of E_tot per application, i.e. ~6.6e6 cycles to drift 1% of the
            // box energy (the flagship runs O(1e3-1e4) cycles).
            //
            // This also removes the mechanism behind the flagship's evacuated "hole": the
            // projection was systematically HEATING magnetized low-density cells (+1.5e-3
            // relative per application at ME/IE>1), driving expansion -> lower rho -> higher
            // ME/IE -> more heating. The gas there ended up 154x hotter than adiabatic
            // expansion permits, which is what identified an unphysical energy source.
            // See DEV_LOG 2026-07-29.
            cons(b, IEN, k, j, i) = etot + (me_post - me_pre);
          } else if (guard && eint_pre > 0.0 && eint_post < eint_guard * eint_pre) {
            // Legacy one-sided guard (superseded by energy_neutral, which makes it inert by
            // construction since eint is then preserved exactly). Only act where the
            // pre-projection state was itself sane (eint_pre > 0); if the conservative update
            // already produced a negative internal energy this is not a CT bookkeeping
            // problem and the existing floor logic must keep ownership of it.
            cons(b, IEN, k, j, i) = e_kin + me_post + eint_guard * eint_pre;
          }
        }
      });
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! History reductions over the projection diagnostic "ct.dEint" (see CT_ProjectBfToCC).
//! ct_projEintMin = the most NEGATIVE relative internal-energy change the projection imposed
//!                  anywhere this cycle (the damaging direction: projection cooling the gas).
//! ct_projEintMax = the largest |relative change| in either direction.
//! Both are 0 when hydro/ct_proj_diag is false.
Real CT_ProjEintMin(MeshData<Real> *md) {
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  if (!hydro_pkg->Param<bool>("ct_proj_diag")) return 0.0;
  auto diag = md->PackVariables(std::vector<std::string>{"ct.dEint"});
  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  Real lmin = 0.0;
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "CT_ProjEintMin", parthenon::DevExecSpace(), 0,
      diag.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &m) {
        m = std::min(m, diag(b, 0, k, j, i));
      },
      Kokkos::Min<Real>(lmin));
  return lmin;
}

Real CT_ProjEintMaxAbs(MeshData<Real> *md) {
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  if (!hydro_pkg->Param<bool>("ct_proj_diag")) return 0.0;
  auto diag = md->PackVariables(std::vector<std::string>{"ct.dEint"});
  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  Real lmax = 0.0;
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "CT_ProjEintMaxAbs", parthenon::DevExecSpace(),
      0, diag.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &m) {
        m = std::max(m, std::abs(diag(b, 0, k, j, i)));
      },
      Kokkos::Max<Real>(lmax));
  return lmax;
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
// NOTE (superseded): this file used to carry FaceCurrentAndB_X{1,2,3} + AmbiFaceEMF_X{1,2,3}
// here -- a "assemble the ambipolar perp-current EMF at 4 cell-prim-based faces, then
// arithmetic-average into the shared edge" construction, mirroring the ideal Balsara-Spicer
// CT_AssembleEMF pattern. That design MIXES, at every edge, one pair of face estimates built
// from a tight 1-cell (dx) derivative of By (the "own-direction" faces) with another pair
// built from a wider 2-cell (2dx) central difference (the "cross-direction" faces, needed
// because those faces don't sample the derivative directly) -- diluting the accurate estimate
// with the coarser one. A first-principles comparison against the validated GLM face-EMF
// (see runs/ct_tests/AMBIPOLAR_CT_UNDERDIFFUSION.md) showed this systematically WEAKENS the
// diffusive EMF at grid-scale (~1-2 dx) diffusivity features by up to ~13%, saturating at
// that floor for arbitrarily large eta contrast and vanishing (<1%) once the feature spans
// >~15 dx -- a per-step deficit that compounds over the many diffusion times a persistent
// current sheet experiences into an order-of-magnitude under-diffusion (the flagship CT
// core-edge ME/E runaway). This is exactly the class of defect CT_AddHallEMF's "COMPACT
// edge-current" comment above already documents fixing for Hall (replacing an earlier
// face-EMF-averaged Hall implementation with a tight, Bf-based edge current). Ambipolar was
// never upgraded to that pattern; CT_AddAmbipolarEMF below now uses it (edge-tight current
// via HallJxE1/HallJyE2/HallJzE3, which are just "current at an edge" and not Hall-specific,
// plus a shared NonidealEdgeEta helper), eliminating the wide/narrow-stencil mix entirely.
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
//! COMPACT edge-current Hall (fixes the whistler-dispersion corruption of the earlier
//! face-EMF-averaged Hall). The current is evaluated with the tight curl of the face field
//! on the SAME edges the Ohmic term uses -- one value per edge, no averaging of the
//! derivative:
//!   J_x on E1 edges (i,j-1/2,k-1/2) = dBz/dy - dBy/dz
//!   J_y on E2 edges (i-1/2,j,k-1/2) = dBx/dz - dBz/dx
//!   J_z on E3 edges (i-1/2,j-1/2,k) = dBy/dx - dBx/dy
//! For the E_H component on a given edge we need the full J and B vectors THERE, so the two
//! transverse J components are interpolated from their own natural edges to the target edge
//! -- but ONLY in the directions transverse to each component's derivative. E.g. for the
//! whistler in x, the E3-edge EMF needs J_y = -dBz/dx (natural on E2 edges, which share the
//! x=i-1/2 location with the E3 edge), so J_y is interpolated only in y,z: the compact
//! dBz/dx stencil survives un-filtered, unlike the box-filter of the face-averaged scheme.
template <class Pack, class Coords>
KOKKOS_INLINE_FUNCTION Real HallJxE1(const Pack &p, const Coords &c, const int b,
                                     const int k, const int j, const int i,
                                     const bool three_d) {
  const Real dBz_dy = (p(b, TE::F3, Bf(), k, j, i) - p(b, TE::F3, Bf(), k, j - 1, i)) /
                      c.template Dxc<X2DIR>(k, j, i);
  const Real dBy_dz =
      three_d ? (p(b, TE::F2, Bf(), k, j, i) - p(b, TE::F2, Bf(), k - 1, j, i)) /
                    c.template Dxc<X3DIR>(k, j, i)
              : 0.0;
  return dBz_dy - dBy_dz;
}
template <class Pack, class Coords>
KOKKOS_INLINE_FUNCTION Real HallJyE2(const Pack &p, const Coords &c, const int b,
                                     const int k, const int j, const int i,
                                     const bool three_d) {
  const Real dBx_dz =
      three_d ? (p(b, TE::F1, Bf(), k, j, i) - p(b, TE::F1, Bf(), k - 1, j, i)) /
                    c.template Dxc<X3DIR>(k, j, i)
              : 0.0;
  const Real dBz_dx = (p(b, TE::F3, Bf(), k, j, i) - p(b, TE::F3, Bf(), k, j, i - 1)) /
                      c.template Dxc<X1DIR>(k, j, i);
  return dBx_dz - dBz_dx;
}
template <class Pack, class Coords>
KOKKOS_INLINE_FUNCTION Real HallJzE3(const Pack &p, const Coords &c, const int b,
                                     const int k, const int j, const int i) {
  const Real dBy_dx = (p(b, TE::F2, Bf(), k, j, i) - p(b, TE::F2, Bf(), k, j, i - 1)) /
                      c.template Dxc<X1DIR>(k, j, i);
  const Real dBx_dy = (p(b, TE::F1, Bf(), k, j, i) - p(b, TE::F1, Bf(), k, j - 1, i)) /
                      c.template Dxc<X2DIR>(k, j, i);
  return dBy_dx - dBx_dy;
}

//! eta_H at an edge as the mean of its four bounding cells (cached NonidealEtaIdx::H, or
//! evaluated from the cell-mean rho/T/x_e and the edge |B|). The four (kc,jc,ic) triplets
//! are the cells sharing the edge.
template <class Prim, class Eta>
KOKKOS_INLINE_FUNCTION Real
HallEdgeEta(const HallDiffusivity &hd, const Prim &prim, const Eta &eta_pack,
            const bool use_cache, const Real bmag, const int b, const int k0, const int j0,
            const int i0, const int k1, const int j1, const int i1, const int k2,
            const int j2, const int i2, const int k3, const int j3, const int i3) {
  if (use_cache) {
    return 0.25 * (eta_pack(b, NonidealEtaIdx::H, k0, j0, i0) +
                   eta_pack(b, NonidealEtaIdx::H, k1, j1, i1) +
                   eta_pack(b, NonidealEtaIdx::H, k2, j2, i2) +
                   eta_pack(b, NonidealEtaIdx::H, k3, j3, i3));
  }
  const Real rho = 0.25 * (prim(b, IDN, k0, j0, i0) + prim(b, IDN, k1, j1, i1) +
                           prim(b, IDN, k2, j2, i2) + prim(b, IDN, k3, j3, i3));
  const Real prs = 0.25 * (prim(b, IPR, k0, j0, i0) + prim(b, IPR, k1, j1, i1) +
                           prim(b, IPR, k2, j2, i2) + prim(b, IPR, k3, j3, i3));
  const int i_xe = hd.XeIndex();
  const Real xe = (i_xe >= 0) ? 0.25 * (prim(b, i_xe, k0, j0, i0) + prim(b, i_xe, k1, j1, i1) +
                                        prim(b, i_xe, k2, j2, i2) + prim(b, i_xe, k3, j3, i3))
                              : -1.0;
  return hd.Get(bmag, rho, prs / rho, xe);
}

//! Generic edge eta (mean of the four bounding cells; cached NonidealEtaIdx entry, or
//! evaluated from the cell-mean rho/T/x_e and the edge |B|) for any non-ideal term whose
//! Diffusivity type exposes Get(bmag,rho,temp,xe)/XeIndex() (Ohmic/Ambipolar/HallDiffusivity
//! all do). Identical in structure to HallEdgeEta above, generalized over eta_idx and the
//! diffusivity type so CT_AddAmbipolarEMF can reuse it without depending on HallDiffusivity.
template <int eta_idx, class Diff, class Prim, class Eta>
KOKKOS_INLINE_FUNCTION Real
NonidealEdgeEta(const Diff &diff, const Prim &prim, const Eta &eta_pack, const bool use_cache,
                const Real bmag, const int b, const int k0, const int j0, const int i0,
                const int k1, const int j1, const int i1, const int k2, const int j2,
                const int i2, const int k3, const int j3, const int i3) {
  if (use_cache) {
    return 0.25 * (eta_pack(b, eta_idx, k0, j0, i0) + eta_pack(b, eta_idx, k1, j1, i1) +
                   eta_pack(b, eta_idx, k2, j2, i2) + eta_pack(b, eta_idx, k3, j3, i3));
  }
  const Real rho = 0.25 * (prim(b, IDN, k0, j0, i0) + prim(b, IDN, k1, j1, i1) +
                           prim(b, IDN, k2, j2, i2) + prim(b, IDN, k3, j3, i3));
  const Real prs = 0.25 * (prim(b, IPR, k0, j0, i0) + prim(b, IPR, k1, j1, i1) +
                           prim(b, IPR, k2, j2, i2) + prim(b, IPR, k3, j3, i3));
  const int i_xe = diff.XeIndex();
  const Real xe = (i_xe >= 0) ? 0.25 * (prim(b, i_xe, k0, j0, i0) + prim(b, i_xe, k1, j1, i1) +
                                        prim(b, i_xe, k2, j2, i2) + prim(b, i_xe, k3, j3, i3))
                              : -1.0;
  return diff.Get(bmag, rho, prs / rho, xe);
}

//----------------------------------------------------------------------------------------
//! Ambipolar perp-current EMF E = eta_A (J - (J.bhat)bhat) evaluated AT A CELL FACE, byte-for-
//! byte the validated GLM construction (AmbipolarDiffFluxIsoFixed / PerpCurrentEMF): the full
//! current vector and B are CO-LOCATED at the face, with the face-normal derivative TIGHT
//! (1 dx -- grid-scale aware) and the two transverse derivatives wide-central (2 dx, symmetric
//! -> damps the grid-scale null mode). This is dissipative + stable (GLM reaches first core
//! with it). CT_AddAmbipolarEMF then builds each edge EMF as the mean of the 4 surrounding
//! FACE EMFs (flux-CT / Balsara-Spicer, same averaging ideal CT uses for the Riemann fluxes),
//! so the edge EMF is single-valued (div(curl E)=0 to round-off) yet keeps GLM's tight,
//! grid-scale-aware dissipation. This REPLACES (i) the "compact edge current + transverse J
//! interpolated from OFFSET edges" scheme (grid-scale dispersive INSTABILITY on tangled
//! fields: OT+AD crash, flagship ME/E runaway) and (ii) a cell-centred central-difference
//! variant (STABLE but grid-scale BLIND -> ~0 dissipation on the sharp core-edge sheets). The
//! 2 dx wide differences use 2*Dxc (exact on the uniform-within-block grid). Uniform-eta and
//! 1D gates are indifferent (all faces identical there). See runs/ct_tests + DEV_LOG 07-29.
template <class Prim, class Eta, class Coords, class Diff>
KOKKOS_INLINE_FUNCTION void
AmbiFaceEMF_X1(const Prim &prim, const Eta &eta_pack, const Coords &c, const Diff &ad_diff,
               const bool use_cache, const int ndim, const int b, const int k, const int j,
               const int i, Real &e1, Real &e2, Real &e3) {
  const Real inv2dy = (ndim > 1) ? 0.5 / c.template Dxc<X2DIR>(k, j, i) : 0.0;
  const Real inv2dz = (ndim > 2) ? 0.5 / c.template Dxc<X3DIR>(k, j, i) : 0.0;
  const Real invdx = 1.0 / c.template Dxc<X1DIR>(k, j, i);
  const auto d1B3 = (prim(b, IB3, k, j, i) - prim(b, IB3, k, j, i - 1)) * invdx;
  const auto d1B2 = (prim(b, IB2, k, j, i) - prim(b, IB2, k, j, i - 1)) * invdx;
  const auto d3B1 =
      ndim > 2 ? (0.5 * (prim(b, IB1, k + 1, j, i - 1) + prim(b, IB1, k + 1, j, i)) -
                  0.5 * (prim(b, IB1, k - 1, j, i - 1) + prim(b, IB1, k - 1, j, i))) *
                     inv2dz
               : 0.0;
  const auto d3B2 =
      ndim > 2 ? (0.5 * (prim(b, IB2, k + 1, j, i - 1) + prim(b, IB2, k + 1, j, i)) -
                  0.5 * (prim(b, IB2, k - 1, j, i - 1) + prim(b, IB2, k - 1, j, i))) *
                     inv2dz
               : 0.0;
  const auto d2B1 =
      ndim > 1 ? (0.5 * (prim(b, IB1, k, j + 1, i - 1) + prim(b, IB1, k, j + 1, i)) -
                  0.5 * (prim(b, IB1, k, j - 1, i - 1) + prim(b, IB1, k, j - 1, i))) *
                     inv2dy
               : 0.0;
  const auto d2B3 =
      ndim > 1 ? (0.5 * (prim(b, IB3, k, j + 1, i - 1) + prim(b, IB3, k, j + 1, i)) -
                  0.5 * (prim(b, IB3, k, j - 1, i - 1) + prim(b, IB3, k, j - 1, i))) *
                     inv2dy
               : 0.0;
  const auto j1 = d2B3 - d3B2;
  const auto j2 = d3B1 - d1B3;
  const auto j3 = d1B2 - d2B1;
  const Real b1 = 0.5 * (prim(b, IB1, k, j, i - 1) + prim(b, IB1, k, j, i));
  const Real b2 = 0.5 * (prim(b, IB2, k, j, i - 1) + prim(b, IB2, k, j, i));
  const Real b3 = 0.5 * (prim(b, IB3, k, j, i - 1) + prim(b, IB3, k, j, i));
  Real eta;
  if (use_cache) {
    eta = 0.5 * (eta_pack(b, NonidealEtaIdx::A, k, j, i - 1) +
                 eta_pack(b, NonidealEtaIdx::A, k, j, i));
  } else {
    const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
    const int i_xe = ad_diff.XeIndex();
    const Real xe = (i_xe >= 0) ? 0.5 * (prim(b, i_xe, k, j, i - 1) + prim(b, i_xe, k, j, i))
                                : -1.0;
    eta = ad_diff.Get(bmag, 0.5 * (prim(b, IDN, k, j, i - 1) + prim(b, IDN, k, j, i)),
                      0.5 * (prim(b, IPR, k, j, i - 1) + prim(b, IPR, k, j, i)) /
                          (0.5 * (prim(b, IDN, k, j, i - 1) + prim(b, IDN, k, j, i))),
                      xe);
  }
  ADPerpEMF(eta, j1, j2, j3, b1, b2, b3, e1, e2, e3);
}

//! Ambipolar perp-EMF at the j-1/2 (X2) face; GLM AmbipolarDiffFluxIsoFixed X2 stencil.
template <class Prim, class Eta, class Coords, class Diff>
KOKKOS_INLINE_FUNCTION void
AmbiFaceEMF_X2(const Prim &prim, const Eta &eta_pack, const Coords &c, const Diff &ad_diff,
               const bool use_cache, const int ndim, const int b, const int k, const int j,
               const int i, Real &e1, Real &e2, Real &e3) {
  const Real inv2dx = 0.5 / c.template Dxc<X1DIR>(k, j, i);
  const Real inv2dz = (ndim > 2) ? 0.5 / c.template Dxc<X3DIR>(k, j, i) : 0.0;
  const Real invdy = 1.0 / c.template Dxc<X2DIR>(k, j, i);
  const auto d2B1 = (prim(b, IB1, k, j, i) - prim(b, IB1, k, j - 1, i)) * invdy;
  const auto d2B3 = (prim(b, IB3, k, j, i) - prim(b, IB3, k, j - 1, i)) * invdy;
  const auto d1B2 = (0.5 * (prim(b, IB2, k, j - 1, i + 1) + prim(b, IB2, k, j, i + 1)) -
                     0.5 * (prim(b, IB2, k, j - 1, i - 1) + prim(b, IB2, k, j, i - 1))) *
                    inv2dx;
  const auto d1B3 = (0.5 * (prim(b, IB3, k, j - 1, i + 1) + prim(b, IB3, k, j, i + 1)) -
                     0.5 * (prim(b, IB3, k, j - 1, i - 1) + prim(b, IB3, k, j, i - 1))) *
                    inv2dx;
  const auto d3B1 =
      ndim > 2 ? (0.5 * (prim(b, IB1, k + 1, j - 1, i) + prim(b, IB1, k + 1, j, i)) -
                  0.5 * (prim(b, IB1, k - 1, j - 1, i) + prim(b, IB1, k - 1, j, i))) *
                     inv2dz
               : 0.0;
  const auto d3B2 =
      ndim > 2 ? (0.5 * (prim(b, IB2, k + 1, j - 1, i) + prim(b, IB2, k + 1, j, i)) -
                  0.5 * (prim(b, IB2, k - 1, j - 1, i) + prim(b, IB2, k - 1, j, i))) *
                     inv2dz
               : 0.0;
  const auto j1 = d2B3 - d3B2;
  const auto j2 = d3B1 - d1B3;
  const auto j3 = d1B2 - d2B1;
  const Real b1 = 0.5 * (prim(b, IB1, k, j - 1, i) + prim(b, IB1, k, j, i));
  const Real b2 = 0.5 * (prim(b, IB2, k, j - 1, i) + prim(b, IB2, k, j, i));
  const Real b3 = 0.5 * (prim(b, IB3, k, j - 1, i) + prim(b, IB3, k, j, i));
  Real eta;
  if (use_cache) {
    eta = 0.5 * (eta_pack(b, NonidealEtaIdx::A, k, j - 1, i) +
                 eta_pack(b, NonidealEtaIdx::A, k, j, i));
  } else {
    const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
    const int i_xe = ad_diff.XeIndex();
    const Real xe = (i_xe >= 0) ? 0.5 * (prim(b, i_xe, k, j - 1, i) + prim(b, i_xe, k, j, i))
                                : -1.0;
    eta = ad_diff.Get(bmag, 0.5 * (prim(b, IDN, k, j - 1, i) + prim(b, IDN, k, j, i)),
                      0.5 * (prim(b, IPR, k, j - 1, i) + prim(b, IPR, k, j, i)) /
                          (0.5 * (prim(b, IDN, k, j - 1, i) + prim(b, IDN, k, j, i))),
                      xe);
  }
  ADPerpEMF(eta, j1, j2, j3, b1, b2, b3, e1, e2, e3);
}

//! Ambipolar perp-EMF at the k-1/2 (X3) face; GLM AmbipolarDiffFluxIsoFixed X3 stencil (3D).
template <class Prim, class Eta, class Coords, class Diff>
KOKKOS_INLINE_FUNCTION void
AmbiFaceEMF_X3(const Prim &prim, const Eta &eta_pack, const Coords &c, const Diff &ad_diff,
               const bool use_cache, const int ndim, const int b, const int k, const int j,
               const int i, Real &e1, Real &e2, Real &e3) {
  const Real inv2dx = 0.5 / c.template Dxc<X1DIR>(k, j, i);
  const Real inv2dy = 0.5 / c.template Dxc<X2DIR>(k, j, i);
  const Real invdz = 1.0 / c.template Dxc<X3DIR>(k, j, i);
  const auto d3B1 = (prim(b, IB1, k, j, i) - prim(b, IB1, k - 1, j, i)) * invdz;
  const auto d3B2 = (prim(b, IB2, k, j, i) - prim(b, IB2, k - 1, j, i)) * invdz;
  const auto d2B3 = (0.5 * (prim(b, IB3, k - 1, j + 1, i) + prim(b, IB3, k, j + 1, i)) -
                     0.5 * (prim(b, IB3, k - 1, j - 1, i) + prim(b, IB3, k, j - 1, i))) *
                    inv2dy;
  const auto d2B1 = (0.5 * (prim(b, IB1, k - 1, j + 1, i) + prim(b, IB1, k, j + 1, i)) -
                     0.5 * (prim(b, IB1, k - 1, j - 1, i) + prim(b, IB1, k, j - 1, i))) *
                    inv2dy;
  const auto d1B3 = (0.5 * (prim(b, IB3, k - 1, j, i + 1) + prim(b, IB3, k, j, i + 1)) -
                     0.5 * (prim(b, IB3, k - 1, j, i - 1) + prim(b, IB3, k, j, i - 1))) *
                    inv2dx;
  const auto d1B2 = (0.5 * (prim(b, IB2, k - 1, j, i + 1) + prim(b, IB2, k, j, i + 1)) -
                     0.5 * (prim(b, IB2, k - 1, j, i - 1) + prim(b, IB2, k, j, i - 1))) *
                    inv2dx;
  const auto j1 = d2B3 - d3B2;
  const auto j2 = d3B1 - d1B3;
  const auto j3 = d1B2 - d2B1;
  const Real b1 = 0.5 * (prim(b, IB1, k - 1, j, i) + prim(b, IB1, k, j, i));
  const Real b2 = 0.5 * (prim(b, IB2, k - 1, j, i) + prim(b, IB2, k, j, i));
  const Real b3 = 0.5 * (prim(b, IB3, k - 1, j, i) + prim(b, IB3, k, j, i));
  Real eta;
  if (use_cache) {
    eta = 0.5 * (eta_pack(b, NonidealEtaIdx::A, k - 1, j, i) +
                 eta_pack(b, NonidealEtaIdx::A, k, j, i));
  } else {
    const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
    const int i_xe = ad_diff.XeIndex();
    const Real xe = (i_xe >= 0) ? 0.5 * (prim(b, i_xe, k - 1, j, i) + prim(b, i_xe, k, j, i))
                                : -1.0;
    eta = ad_diff.Get(bmag, 0.5 * (prim(b, IDN, k - 1, j, i) + prim(b, IDN, k, j, i)),
                      0.5 * (prim(b, IPR, k - 1, j, i) + prim(b, IPR, k, j, i)) /
                          (0.5 * (prim(b, IDN, k - 1, j, i) + prim(b, IDN, k, j, i))),
                      xe);
  }
  ADPerpEMF(eta, j1, j2, j3, b1, b2, b3, e1, e2, e3);
}

//----------------------------------------------------------------------------------------
//! DIRECT edge-centred ambipolar perp-EMF from cell-centred prim B -- the "third way" that is
//! simultaneously (i) GRID-SCALE AWARE: the two IN-PLANE derivatives (transverse to the edge)
//! use TIGHT 1-cell differences (averaged over the one out-of-edge-plane index that the edge
//! straddles), so a single-cell current sheet is fully resolved -- unlike the flux-CT face
//! average, which blends a tight and a wide estimate of the SAME derivative and smooths it;
//! (ii) STABLE + dissipative: the full J and B are CO-LOCATED at the edge in ONE PerpCurrentEMF
//! evaluation (E.J=eta|J_perp|^2>=0 pointwise) -- unlike the compact scheme, whose transverse J
//! came from OFFSET neighbour edges and decohered J -> a grid-scale dispersive instability;
//! (iii) DIV-FREE: one value per edge -> div(curl E)=0 to round-off. The one OUT-OF-PLANE
//! derivative (along the edge) is a wide 2-cell central difference (2*Dxc, exact on the
//! uniform-within-block grid), averaged over the 4 cells sharing the edge -- the same benign
//! wide-transverse treatment GLM uses at faces. This targets the confirmed defect (per-level
//! B.M: flux-CT matched GLM on the resolved levels L2-L5 but under-dissipated 39x at the
//! finest level L7 and went anti-dissipative at L6, purely from the face-average smoothing of
//! the grid-scale core-edge sheet). eta uses the shared 4-cell NonidealEdgeEta. See DEV_LOG
//! 2026-07-29. (The unused AmbiFaceEMF_X* / compact helpers above are kept for reference.)

//! e3 at the E3 edge (i-1/2,j-1/2,k). In-plane = xy (tight); along-edge = z (wide, 3D only).
template <class Prim, class Eta, class Coords, class Diff>
KOKKOS_INLINE_FUNCTION Real
AmbiEdgeEMF_E3(const Prim &prim, const Eta &eta_pack, const Coords &c, const Diff &ad_diff,
               const bool use_cache, const int ndim, const int b, const int k, const int j,
               const int i) {
  const Real hix = 0.5 / c.template Dxc<X1DIR>(k, j, i);
  const Real hiy = 0.5 / c.template Dxc<X2DIR>(k, j, i);
  // TIGHT in-plane: 1-cell x/y diffs averaged over the other in-plane index (j-1..j / i-1..i).
  const Real dBy_dx = hix * ((prim(b, IB2, k, j, i) - prim(b, IB2, k, j, i - 1)) +
                             (prim(b, IB2, k, j - 1, i) - prim(b, IB2, k, j - 1, i - 1)));
  const Real dBx_dy = hiy * ((prim(b, IB1, k, j, i) - prim(b, IB1, k, j - 1, i)) +
                             (prim(b, IB1, k, j, i - 1) - prim(b, IB1, k, j - 1, i - 1)));
  const Real dBz_dx = hix * ((prim(b, IB3, k, j, i) - prim(b, IB3, k, j, i - 1)) +
                             (prim(b, IB3, k, j - 1, i) - prim(b, IB3, k, j - 1, i - 1)));
  const Real dBz_dy = hiy * ((prim(b, IB3, k, j, i) - prim(b, IB3, k, j - 1, i)) +
                             (prim(b, IB3, k, j, i - 1) - prim(b, IB3, k, j - 1, i - 1)));
  Real dBx_dz = 0.0, dBy_dz = 0.0;
  if (ndim > 2) {
    const Real qiz = 0.25 * 0.5 / c.template Dxc<X3DIR>(k, j, i);
    dBx_dz = qiz * ((prim(b, IB1, k + 1, j, i) - prim(b, IB1, k - 1, j, i)) +
                    (prim(b, IB1, k + 1, j, i - 1) - prim(b, IB1, k - 1, j, i - 1)) +
                    (prim(b, IB1, k + 1, j - 1, i) - prim(b, IB1, k - 1, j - 1, i)) +
                    (prim(b, IB1, k + 1, j - 1, i - 1) - prim(b, IB1, k - 1, j - 1, i - 1)));
    dBy_dz = qiz * ((prim(b, IB2, k + 1, j, i) - prim(b, IB2, k - 1, j, i)) +
                    (prim(b, IB2, k + 1, j, i - 1) - prim(b, IB2, k - 1, j, i - 1)) +
                    (prim(b, IB2, k + 1, j - 1, i) - prim(b, IB2, k - 1, j - 1, i)) +
                    (prim(b, IB2, k + 1, j - 1, i - 1) - prim(b, IB2, k - 1, j - 1, i - 1)));
  }
  const Real Jx = dBz_dy - dBy_dz;
  const Real Jy = dBx_dz - dBz_dx;
  const Real Jz = dBy_dx - dBx_dy;
  const Real Bx = 0.25 * (prim(b, IB1, k, j, i) + prim(b, IB1, k, j, i - 1) +
                          prim(b, IB1, k, j - 1, i) + prim(b, IB1, k, j - 1, i - 1));
  const Real By = 0.25 * (prim(b, IB2, k, j, i) + prim(b, IB2, k, j, i - 1) +
                          prim(b, IB2, k, j - 1, i) + prim(b, IB2, k, j - 1, i - 1));
  const Real Bz = 0.25 * (prim(b, IB3, k, j, i) + prim(b, IB3, k, j, i - 1) +
                          prim(b, IB3, k, j - 1, i) + prim(b, IB3, k, j - 1, i - 1));
  const Real bmag = std::sqrt(SQR(Bx) + SQR(By) + SQR(Bz));
  const Real eta = NonidealEdgeEta<NonidealEtaIdx::A>(
      ad_diff, prim, eta_pack, use_cache, bmag, b, k, j, i, k, j, i - 1, k, j - 1, i, k,
      j - 1, i - 1);
  Real e1, e2, e3;
  ADPerpEMF(eta, Jx, Jy, Jz, Bx, By, Bz, e1, e2, e3);
  return e3;
}

//! e1 at the E1 edge (i,j-1/2,k-1/2). In-plane = yz (tight); along-edge = x (wide). 3D only.
template <class Prim, class Eta, class Coords, class Diff>
KOKKOS_INLINE_FUNCTION Real
AmbiEdgeEMF_E1(const Prim &prim, const Eta &eta_pack, const Coords &c, const Diff &ad_diff,
               const bool use_cache, const int b, const int k, const int j, const int i) {
  const Real hiy = 0.5 / c.template Dxc<X2DIR>(k, j, i);
  const Real hiz = 0.5 / c.template Dxc<X3DIR>(k, j, i);
  const Real dBz_dy = hiy * ((prim(b, IB3, k, j, i) - prim(b, IB3, k, j - 1, i)) +
                             (prim(b, IB3, k - 1, j, i) - prim(b, IB3, k - 1, j - 1, i)));
  const Real dBy_dz = hiz * ((prim(b, IB2, k, j, i) - prim(b, IB2, k - 1, j, i)) +
                             (prim(b, IB2, k, j - 1, i) - prim(b, IB2, k - 1, j - 1, i)));
  const Real dBx_dy = hiy * ((prim(b, IB1, k, j, i) - prim(b, IB1, k, j - 1, i)) +
                             (prim(b, IB1, k - 1, j, i) - prim(b, IB1, k - 1, j - 1, i)));
  const Real dBx_dz = hiz * ((prim(b, IB1, k, j, i) - prim(b, IB1, k - 1, j, i)) +
                             (prim(b, IB1, k, j - 1, i) - prim(b, IB1, k - 1, j - 1, i)));
  const Real qix = 0.25 * 0.5 / c.template Dxc<X1DIR>(k, j, i);
  const Real dBy_dx = qix * ((prim(b, IB2, k, j, i + 1) - prim(b, IB2, k, j, i - 1)) +
                             (prim(b, IB2, k, j - 1, i + 1) - prim(b, IB2, k, j - 1, i - 1)) +
                             (prim(b, IB2, k - 1, j, i + 1) - prim(b, IB2, k - 1, j, i - 1)) +
                             (prim(b, IB2, k - 1, j - 1, i + 1) -
                              prim(b, IB2, k - 1, j - 1, i - 1)));
  const Real dBz_dx = qix * ((prim(b, IB3, k, j, i + 1) - prim(b, IB3, k, j, i - 1)) +
                             (prim(b, IB3, k, j - 1, i + 1) - prim(b, IB3, k, j - 1, i - 1)) +
                             (prim(b, IB3, k - 1, j, i + 1) - prim(b, IB3, k - 1, j, i - 1)) +
                             (prim(b, IB3, k - 1, j - 1, i + 1) -
                              prim(b, IB3, k - 1, j - 1, i - 1)));
  const Real Jx = dBz_dy - dBy_dz;
  const Real Jy = dBx_dz - dBz_dx;
  const Real Jz = dBy_dx - dBx_dy;
  const Real Bx = 0.25 * (prim(b, IB1, k, j, i) + prim(b, IB1, k, j - 1, i) +
                          prim(b, IB1, k - 1, j, i) + prim(b, IB1, k - 1, j - 1, i));
  const Real By = 0.25 * (prim(b, IB2, k, j, i) + prim(b, IB2, k, j - 1, i) +
                          prim(b, IB2, k - 1, j, i) + prim(b, IB2, k - 1, j - 1, i));
  const Real Bz = 0.25 * (prim(b, IB3, k, j, i) + prim(b, IB3, k, j - 1, i) +
                          prim(b, IB3, k - 1, j, i) + prim(b, IB3, k - 1, j - 1, i));
  const Real bmag = std::sqrt(SQR(Bx) + SQR(By) + SQR(Bz));
  const Real eta = NonidealEdgeEta<NonidealEtaIdx::A>(
      ad_diff, prim, eta_pack, use_cache, bmag, b, k, j, i, k, j - 1, i, k - 1, j, i, k - 1,
      j - 1, i);
  Real e1, e2, e3;
  ADPerpEMF(eta, Jx, Jy, Jz, Bx, By, Bz, e1, e2, e3);
  return e1;
}

//! e2 at the E2 edge (i-1/2,j,k-1/2). In-plane = xz (tight); along-edge = y (wide). 3D only.
template <class Prim, class Eta, class Coords, class Diff>
KOKKOS_INLINE_FUNCTION Real
AmbiEdgeEMF_E2(const Prim &prim, const Eta &eta_pack, const Coords &c, const Diff &ad_diff,
               const bool use_cache, const int b, const int k, const int j, const int i) {
  const Real hix = 0.5 / c.template Dxc<X1DIR>(k, j, i);
  const Real hiz = 0.5 / c.template Dxc<X3DIR>(k, j, i);
  const Real dBz_dx = hix * ((prim(b, IB3, k, j, i) - prim(b, IB3, k, j, i - 1)) +
                             (prim(b, IB3, k - 1, j, i) - prim(b, IB3, k - 1, j, i - 1)));
  const Real dBx_dz = hiz * ((prim(b, IB1, k, j, i) - prim(b, IB1, k - 1, j, i)) +
                             (prim(b, IB1, k, j, i - 1) - prim(b, IB1, k - 1, j, i - 1)));
  const Real dBy_dx = hix * ((prim(b, IB2, k, j, i) - prim(b, IB2, k, j, i - 1)) +
                             (prim(b, IB2, k - 1, j, i) - prim(b, IB2, k - 1, j, i - 1)));
  const Real dBy_dz = hiz * ((prim(b, IB2, k, j, i) - prim(b, IB2, k - 1, j, i)) +
                             (prim(b, IB2, k, j, i - 1) - prim(b, IB2, k - 1, j, i - 1)));
  const Real qiy = 0.25 * 0.5 / c.template Dxc<X2DIR>(k, j, i);
  const Real dBx_dy = qiy * ((prim(b, IB1, k, j + 1, i) - prim(b, IB1, k, j - 1, i)) +
                             (prim(b, IB1, k, j + 1, i - 1) - prim(b, IB1, k, j - 1, i - 1)) +
                             (prim(b, IB1, k - 1, j + 1, i) - prim(b, IB1, k - 1, j - 1, i)) +
                             (prim(b, IB1, k - 1, j + 1, i - 1) -
                              prim(b, IB1, k - 1, j - 1, i - 1)));
  const Real dBz_dy = qiy * ((prim(b, IB3, k, j + 1, i) - prim(b, IB3, k, j - 1, i)) +
                             (prim(b, IB3, k, j + 1, i - 1) - prim(b, IB3, k, j - 1, i - 1)) +
                             (prim(b, IB3, k - 1, j + 1, i) - prim(b, IB3, k - 1, j - 1, i)) +
                             (prim(b, IB3, k - 1, j + 1, i - 1) -
                              prim(b, IB3, k - 1, j - 1, i - 1)));
  const Real Jx = dBz_dy - dBy_dz;
  const Real Jy = dBx_dz - dBz_dx;
  const Real Jz = dBy_dx - dBx_dy;
  const Real Bx = 0.25 * (prim(b, IB1, k, j, i) + prim(b, IB1, k, j, i - 1) +
                          prim(b, IB1, k - 1, j, i) + prim(b, IB1, k - 1, j, i - 1));
  const Real By = 0.25 * (prim(b, IB2, k, j, i) + prim(b, IB2, k, j, i - 1) +
                          prim(b, IB2, k - 1, j, i) + prim(b, IB2, k - 1, j, i - 1));
  const Real Bz = 0.25 * (prim(b, IB3, k, j, i) + prim(b, IB3, k, j, i - 1) +
                          prim(b, IB3, k - 1, j, i) + prim(b, IB3, k - 1, j, i - 1));
  const Real bmag = std::sqrt(SQR(Bx) + SQR(By) + SQR(Bz));
  const Real eta = NonidealEdgeEta<NonidealEtaIdx::A>(
      ad_diff, prim, eta_pack, use_cache, bmag, b, k, j, i, k, j, i - 1, k - 1, j, i, k - 1,
      j, i - 1);
  Real e1, e2, e3;
  ADPerpEMF(eta, Jx, Jy, Jz, Bx, By, Bz, e1, e2, e3);
  return e2;
}

//----------------------------------------------------------------------------------------
//! Non-ideal (ambipolar) edge EMF for CT -- DIRECT co-located edge stencil (AmbiEdgeEMF_E*).
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

  // ---- E3 edge (i-1/2,j-1/2,k): Jz tight (Bf curl); Jx,Jy interpolated from E1/E2. ----
  {
    IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E3);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_AmbiEMF_E3", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(b);
          // E3 edge: DIRECT co-located tight-in-plane edge EMF (no face averaging).
          pack.flux(b, TE::E3, Bf(), k, j, i) +=
              AmbiEdgeEMF_E3(prim, eta_pack, c, ad_diff, use_cache, ndim, b, k, j, i);
        });
  }

  if (three_d) {
    // ---- E1 edge (i,j-1/2,k-1/2): Jx tight; Jy,Jz interpolated from E2/E3. ----
    {
      IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E1);
      IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E1);
      IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E1);
      parthenon::par_for(
          DEFAULT_LOOP_PATTERN, "CT_AmbiEMF_E1", parthenon::DevExecSpace(), 0, nb - 1,
          kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
          KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
            const auto &c = pack.GetCoordinates(b);
            // E1 edge: DIRECT co-located tight-in-plane edge EMF (no face averaging).
            pack.flux(b, TE::E1, Bf(), k, j, i) +=
                AmbiEdgeEMF_E1(prim, eta_pack, c, ad_diff, use_cache, b, k, j, i);
          });
    }
    // ---- E2 edge (i-1/2,j,k-1/2): Jy tight; Jx,Jz interpolated from E1/E3. ----
    {
      IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E2);
      IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E2);
      IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E2);
      parthenon::par_for(
          DEFAULT_LOOP_PATTERN, "CT_AmbiEMF_E2", parthenon::DevExecSpace(), 0, nb - 1,
          kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
          KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
            const auto &c = pack.GetCoordinates(b);
            // E2 edge: DIRECT co-located tight-in-plane edge EMF (no face averaging).
            pack.flux(b, TE::E2, Bf(), k, j, i) +=
                AmbiEdgeEMF_E2(prim, eta_pack, c, ad_diff, use_cache, b, k, j, i);
          });
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Non-ideal (Hall) edge EMF for CT -- COMPACT edge-current formulation (see the helper
//! comment above). The Hall EMF E_H = eta_H (J x B)/|B| (+ optional Ohmic floor eta_O J) is
//! formed ONCE per edge from the tight edge currents and the edge-interpolated field, which
//! preserves the whistler dispersion (the face-EMF-averaged version was ~25-90% too slow).
//! Hall is dispersive and unsplit-only under CT. The matching cons.flux(IBn) deposit in
//! HallDiffFluxIsoFixed is gated off under CT; its cons.flux(IEN) Poynting term stays FV.
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
  // B11: per-cell Ohmic stabilizer = max(absolute floor, ratio*|eta_H|). CT forbids RKL2, so
  // eta_h is always evaluated on these edges and the effective floor is available directly.
  // ratio = 0 reproduces the constant floor exactly.
  const Real hall_floor_ratio = hall_diff.GetOhmicFloorRatio();
  const bool use_cache = hydro_pkg->Param<bool>("nonideal_eta_cache");

  auto prim = md->PackVariables(std::vector<std::string>{"prim"});
  const auto eta_pack =
      md->PackVariables(std::vector<std::string>{use_cache ? "nonideal_eta" : "prim"});
  static auto desc =
      parthenon::MakePackDescriptor<Bf>(md, {}, {parthenon::PDOpt::WithFluxes});
  auto pack = desc.GetPack(md);
  const int nb = pack.GetNBlocks();

  // ---- E3 edge (i-1/2,j-1/2,k): e3 = eta_H (Jx By - Jy Bx)/|B| + eta_floor Jz. 2D & 3D. ----
  {
    IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E3);
    IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E3);
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "CT_HallEMF_E3", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
        kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(b);
          const Real Jz = HallJzE3(pack, c, b, k, j, i);
          // J_x from E1 edges: interp x (i,i-1), z (k,k+1) -> preserves any dBz/dy stencil.
          const Real Jx =
              three_d ? 0.25 * (HallJxE1(pack, c, b, k, j, i, three_d) +
                                HallJxE1(pack, c, b, k, j, i - 1, three_d) +
                                HallJxE1(pack, c, b, k + 1, j, i, three_d) +
                                HallJxE1(pack, c, b, k + 1, j, i - 1, three_d))
                      : 0.5 * (HallJxE1(pack, c, b, k, j, i, three_d) +
                               HallJxE1(pack, c, b, k, j, i - 1, three_d));
          // J_y from E2 edges: interp y (j,j-1), z (k,k+1) -> preserves the dBz/dx stencil.
          const Real Jy =
              three_d ? 0.25 * (HallJyE2(pack, c, b, k, j, i, three_d) +
                                HallJyE2(pack, c, b, k, j - 1, i, three_d) +
                                HallJyE2(pack, c, b, k + 1, j, i, three_d) +
                                HallJyE2(pack, c, b, k + 1, j - 1, i, three_d))
                      : 0.5 * (HallJyE2(pack, c, b, k, j, i, three_d) +
                               HallJyE2(pack, c, b, k, j - 1, i, three_d));
          const Real Bx = 0.5 * (pack(b, TE::F1, Bf(), k, j, i) +
                                 pack(b, TE::F1, Bf(), k, j - 1, i));
          const Real By = 0.5 * (pack(b, TE::F2, Bf(), k, j, i) +
                                 pack(b, TE::F2, Bf(), k, j, i - 1));
          Real Bz = 0.25 * (pack(b, TE::F3, Bf(), k, j, i) +
                            pack(b, TE::F3, Bf(), k, j - 1, i) +
                            pack(b, TE::F3, Bf(), k, j, i - 1) +
                            pack(b, TE::F3, Bf(), k, j - 1, i - 1));
          if (three_d) {
            Bz = 0.5 * Bz + 0.125 * (pack(b, TE::F3, Bf(), k + 1, j, i) +
                                     pack(b, TE::F3, Bf(), k + 1, j - 1, i) +
                                     pack(b, TE::F3, Bf(), k + 1, j, i - 1) +
                                     pack(b, TE::F3, Bf(), k + 1, j - 1, i - 1));
          }
          const Real bmag = std::sqrt(SQR(Bx) + SQR(By) + SQR(Bz));
          const Real eta_h = HallEdgeEta(hall_diff, prim, eta_pack, use_cache, bmag, b, k, j,
                                         i, k, j, i - 1, k, j - 1, i, k, j - 1, i - 1);
          Real e1, e2, e3;
          HallEMFLocal(eta_h,
                       (hall_floor_ratio > 0.0) ? hall_diff.EffectiveOhmicFloor(eta_h)
                                                : eta_floor,
                       bmag, Jx, Jy, Jz, Bx, By, Bz, e1, e2, e3);
          pack.flux(b, TE::E3, Bf(), k, j, i) += e3;
        });
  }

  if (three_d) {
    // ---- E1 edge (i,j-1/2,k-1/2): e1 = eta_H (Jy Bz - Jz By)/|B| + eta_floor Jx. ----
    {
      IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E1);
      IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E1);
      IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E1);
      parthenon::par_for(
          DEFAULT_LOOP_PATTERN, "CT_HallEMF_E1", parthenon::DevExecSpace(), 0, nb - 1,
          kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
          KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
            const auto &c = pack.GetCoordinates(b);
            const Real Jx = HallJxE1(pack, c, b, k, j, i, true);
            const Real Jy = 0.25 * (HallJyE2(pack, c, b, k, j, i, true) +
                                    HallJyE2(pack, c, b, k, j - 1, i, true) +
                                    HallJyE2(pack, c, b, k, j, i + 1, true) +
                                    HallJyE2(pack, c, b, k, j - 1, i + 1, true));
            const Real Jz = 0.25 * (HallJzE3(pack, c, b, k, j, i) +
                                    HallJzE3(pack, c, b, k, j, i + 1) +
                                    HallJzE3(pack, c, b, k - 1, j, i) +
                                    HallJzE3(pack, c, b, k - 1, j, i + 1));
            const Real Bx = 0.125 * (pack(b, TE::F1, Bf(), k, j, i) +
                                     pack(b, TE::F1, Bf(), k, j, i + 1) +
                                     pack(b, TE::F1, Bf(), k, j - 1, i) +
                                     pack(b, TE::F1, Bf(), k, j - 1, i + 1) +
                                     pack(b, TE::F1, Bf(), k - 1, j, i) +
                                     pack(b, TE::F1, Bf(), k - 1, j, i + 1) +
                                     pack(b, TE::F1, Bf(), k - 1, j - 1, i) +
                                     pack(b, TE::F1, Bf(), k - 1, j - 1, i + 1));
            const Real By = 0.5 * (pack(b, TE::F2, Bf(), k, j, i) +
                                   pack(b, TE::F2, Bf(), k - 1, j, i));
            const Real Bz = 0.5 * (pack(b, TE::F3, Bf(), k, j, i) +
                                   pack(b, TE::F3, Bf(), k, j - 1, i));
            const Real bmag = std::sqrt(SQR(Bx) + SQR(By) + SQR(Bz));
            const Real eta_h =
                HallEdgeEta(hall_diff, prim, eta_pack, use_cache, bmag, b, k, j, i, k,
                            j - 1, i, k - 1, j, i, k - 1, j - 1, i);
            Real e1, e2, e3;
            HallEMFLocal(eta_h,
                       (hall_floor_ratio > 0.0) ? hall_diff.EffectiveOhmicFloor(eta_h)
                                                : eta_floor,
                       bmag, Jx, Jy, Jz, Bx, By, Bz, e1, e2, e3);
            pack.flux(b, TE::E1, Bf(), k, j, i) += e1;
          });
    }
    // ---- E2 edge (i-1/2,j,k-1/2): e2 = eta_H (Jz Bx - Jx Bz)/|B| + eta_floor Jy. ----
    {
      IndexRange ib = md->GetBoundsI(CellLevel::same, IndexDomain::interior, TE::E2);
      IndexRange jb = md->GetBoundsJ(CellLevel::same, IndexDomain::interior, TE::E2);
      IndexRange kb = md->GetBoundsK(CellLevel::same, IndexDomain::interior, TE::E2);
      parthenon::par_for(
          DEFAULT_LOOP_PATTERN, "CT_HallEMF_E2", parthenon::DevExecSpace(), 0, nb - 1,
          kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
          KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
            const auto &c = pack.GetCoordinates(b);
            const Real Jy = HallJyE2(pack, c, b, k, j, i, true);
            const Real Jx = 0.25 * (HallJxE1(pack, c, b, k, j, i, true) +
                                    HallJxE1(pack, c, b, k, j + 1, i, true) +
                                    HallJxE1(pack, c, b, k, j, i - 1, true) +
                                    HallJxE1(pack, c, b, k, j + 1, i - 1, true));
            const Real Jz = 0.25 * (HallJzE3(pack, c, b, k, j, i) +
                                    HallJzE3(pack, c, b, k, j + 1, i) +
                                    HallJzE3(pack, c, b, k - 1, j, i) +
                                    HallJzE3(pack, c, b, k - 1, j + 1, i));
            const Real Bx = 0.5 * (pack(b, TE::F1, Bf(), k, j, i) +
                                   pack(b, TE::F1, Bf(), k - 1, j, i));
            const Real By = 0.125 * (pack(b, TE::F2, Bf(), k, j, i) +
                                     pack(b, TE::F2, Bf(), k, j + 1, i) +
                                     pack(b, TE::F2, Bf(), k, j, i - 1) +
                                     pack(b, TE::F2, Bf(), k, j + 1, i - 1) +
                                     pack(b, TE::F2, Bf(), k - 1, j, i) +
                                     pack(b, TE::F2, Bf(), k - 1, j + 1, i) +
                                     pack(b, TE::F2, Bf(), k - 1, j, i - 1) +
                                     pack(b, TE::F2, Bf(), k - 1, j + 1, i - 1));
            const Real Bz = 0.5 * (pack(b, TE::F3, Bf(), k, j, i) +
                                   pack(b, TE::F3, Bf(), k, j, i - 1));
            const Real bmag = std::sqrt(SQR(Bx) + SQR(By) + SQR(Bz));
            const Real eta_h =
                HallEdgeEta(hall_diff, prim, eta_pack, use_cache, bmag, b, k, j, i, k, j,
                            i - 1, k - 1, j, i, k - 1, j, i - 1);
            Real e1, e2, e3;
            HallEMFLocal(eta_h,
                       (hall_floor_ratio > 0.0) ? hall_diff.EffectiveOhmicFloor(eta_h)
                                                : eta_floor,
                       bmag, Jx, Jy, Jz, Bx, By, Bz, e1, e2, e3);
            pack.flux(b, TE::E2, Bf(), k, j, i) += e2;
          });
    }
  }
  return TaskStatus::complete;
}

} // namespace CT
} // namespace Hydro
