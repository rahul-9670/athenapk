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

} // namespace CT
} // namespace Hydro
