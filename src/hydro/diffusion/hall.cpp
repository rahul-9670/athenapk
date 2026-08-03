//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2024, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hall.cpp
//! \brief Hall effect. The Hall diffusivity follows the constant-coefficient model of
//!        Athena++ (eta_H = coeff * B / rho), but unlike Athena++ -- where the Hall EMF is
//!        commented out and never applied -- the EMF E_H = eta_H (J x B) / |B| is actually
//!        added to the induction equation here. The Hall term is dispersive (whistler
//!        waves) and is therefore never super-time-stepped: under diffusion/integrator=
//!        rkl2 the Hall EMF is applied unsplit (added to the hyperbolic fluxes each stage)
//!        with its whistler dt as a strict constraint, while the parabolic terms go into
//!        RKL2 (mixed mode, wired in hydro.cpp / hydro_driver.cpp). An optional Ohmic
//!        resistivity floor (eta_O J) can be added inside the kernel for numerical
//!        stabilization on the cell-centered grid; the eta_h_on/floor_on arguments select
//!        which of the two parts a given call applies.

// Parthenon headers
#include <cmath>
#include <parthenon/package.hpp>

// AthenaPK headers
#include "../../main.hpp"
#include "config.hpp"
#include "diffusion.hpp"
#include "kokkos_abstraction.hpp"
#include "utils/error_checking.hpp"

using namespace parthenon::package::prelude;

// NOTE: HallDiffusivity::Get bodies moved to diffusion.hpp -- they are
// KOKKOS_INLINE_FUNCTION and are called from multiple translation units
// (flux kernels, dt estimators, PrecomputeNonidealEta).

//----------------------------------------------------------------------------------------
//! Hall EMF E = eta_H (J x B)/|B| + eta_floor J (the latter is an optional Ohmic
//! stabilizer).
KOKKOS_INLINE_FUNCTION
void HallEMF(const Real eta_h, const Real eta_floor, const Real bmag, const Real j1,
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

Real EstimateHallTimestep(MeshData<Real> *md, const bool whistler_on,
                          const bool floor_on) {
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real min_dt_hall = std::numeric_limits<Real>::max();
  const auto ndim = prim_pack.GetNdim();

  // Whistler timestep factors (cf. Athena++ NewDiffusionDt): 0.5 in 1D, 1.0 in 2D/3D.
  const Real fac_whistler = (ndim == 1) ? 0.5 : 1.0;
  // The optional Ohmic stabilizer eta_floor*J applied inside HallEMF is a real parabolic
  // diffusion and needs the (stricter) parabolic factor: 0.5 (1D), 0.25 (2D), 1/6 (3D).
  Real fac_par = 0.5;
  if (ndim == 2) {
    fac_par = 0.25;
  } else if (ndim == 3) {
    fac_par = 1.0 / 6.0;
  }

  const auto &hall_diff = hydro_pkg->Param<HallDiffusivity>("hall_diff");
  const Real eta_floor = hall_diff.GetOhmicFloor();
  // B11: with a per-cell floor the stability limit dx^2/eta_floor is also per-cell, and it
  // is the LARGEST floor in the mesh that binds -- so eta_H has to be evaluated here too,
  // even when the whistler part is switched off (mixed rkl2 mode). Gated on the ratio so
  // the ratio=0 path does no extra work and is bit-identical.
  const Real floor_ratio = hall_diff.GetOhmicFloorRatio();
  const bool need_eta_h = whistler_on || (floor_on && floor_ratio > 0.0);

  // The two constraints carry different stability factors, so fold them into the
  // reduction (only cfl_diff is applied outside).
  Kokkos::parallel_reduce(
      "EstimateHallTimestep",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExecSpace(), {0, kb.s, jb.s, ib.s},
          {prim_pack.GetDim(5), kb.e + 1, jb.e + 1, ib.e + 1},
          {1, 1, 1, ib.e + 1 - ib.s}),
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &min_dt) {
        const auto &coords = prim_pack.GetCoords(b);
        const auto &prim = prim_pack(b);
        Real mindx2 = SQR(coords.Dxc<1>(k, j, i));
        if (ndim >= 2) {
          mindx2 = fmin(mindx2, SQR(coords.Dxc<2>(k, j, i)));
        }
        if (ndim >= 3) {
          mindx2 = fmin(mindx2, SQR(coords.Dxc<3>(k, j, i)));
        }
        Real eta_h_cell = 0.0;
        if (need_eta_h) {
          const auto bmag = std::sqrt(SQR(prim(IB1, k, j, i)) + SQR(prim(IB2, k, j, i)) +
                                      SQR(prim(IB3, k, j, i)));
          const auto temp = prim(IPR, k, j, i) / prim(IDN, k, j, i); // code T (p = rho*T)
          const int i_xe = hall_diff.XeIndex();
          const Real xe = (i_xe >= 0) ? prim(i_xe, k, j, i) : -1.0;
          eta_h_cell = hall_diff.Get(bmag, prim(IDN, k, j, i), temp, xe);
        }
        if (whistler_on) {
          const auto eta = std::abs(eta_h_cell);
          min_dt = fmin(min_dt, fac_whistler * mindx2 / (eta + TINY_NUMBER));
        }
        if (floor_on) {
          // B11: per-cell floor => per-cell parabolic constraint.
          const Real eta_f = (floor_ratio > 0.0)
                                 ? hall_diff.EffectiveOhmicFloor(eta_h_cell)
                                 : eta_floor;
          if (eta_f > 0.0) {
            min_dt = fmin(min_dt, fac_par * mindx2 / (eta_f + TINY_NUMBER));
          }
        }
      },
      Kokkos::Min<Real>(min_dt_hall));

  const auto &cfl_diff = hydro_pkg->Param<Real>("cfl_diff");
  return cfl_diff * min_dt_hall;
}

//---------------------------------------------------------------------------------------
//! Calculate the Hall EMF (dispersive part gated by eta_h_on, Ohmic-floor stabilizer by
//! floor_on; see the declaration in diffusion.hpp)

void HallDiffFluxIsoFixed(MeshData<Real> *md, const bool eta_h_on, const bool floor_on) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  // Pack "cons" by NAME (not {Independent}): pins the fixed IB*/IEN flux indices to cons
  // regardless of registration order (v6-STS-bug trap class). Bit-identical today.
  const std::vector<std::string> cons_names{"cons"};
  auto cons_pack = md->PackVariablesAndFluxes(cons_names, cons_names);
  auto hydro_pkg = pmb->packages.Get("Hydro");

  // Under Constrained Transport (divergence_control=ct) the Hall induction is advanced from
  // an edge-centered curl assembled in CT_AddHallEMF (ct.cpp), so the cell-centered
  // cons.flux(IB*) induction deposits below are gated OFF (GS05 would double-count them).
  // The Hall cons.flux(IEN) Poynting term stays on the FV energy flux. GLM path bit-identical.
  // CT forbids RKL2, so under CT this kernel is only ever called unsplit (eta_h_on=floor_on).
  const bool ct_induction = hydro_pkg->Param<bool>("use_ct");

  auto const &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});

  const int ndim = pmb->pmy_mesh->ndim;

  const auto &hall_diff = hydro_pkg->Param<HallDiffusivity>("hall_diff");
  // A disabled part contributes exactly zero to the EMF: the floor via eta_floor = 0,
  // the dispersive part by skipping the eta_H evaluation (eta_h = 0) in the kernels.
  const auto eta_floor_abs = floor_on ? hall_diff.GetOhmicFloor() : 0.0;
  // B11: with a ratio set, the floor is per-cell (max of the absolute floor and
  // ratio*|eta_H|), so eta_H must be evaluated for the FLOOR as well -- including in the
  // mixed rkl2 mode, where this kernel is called once with (eta_h_on=true, floor_on=false)
  // for the dispersive part and again with (eta_h_on=false, floor_on=true) for the floor.
  // Without `need_eta_h` that second call would see eta_h = 0 and silently fall back to the
  // absolute floor -- i.e. the fix would be inert in exactly the production integrator.
  const Real floor_ratio = floor_on ? hall_diff.GetOhmicFloorRatio() : 0.0;
  const bool need_eta_h = eta_h_on || (floor_ratio > 0.0);

  // Cell-centered eta cache (PrecomputeNonidealEta): face eta_H = arithmetic average of
  // the two adjacent cached values. When the cache is off, "prim" is packed as a dummy
  // to keep the lambda capture valid; it is never indexed on that branch.
  const bool use_cache = hydro_pkg->Param<bool>("nonideal_eta_cache");
  const auto eta_pack = md->PackVariables(
      std::vector<std::string>{use_cache ? "nonideal_eta" : "prim"});

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Hall X1 fluxes", DevExecSpace(), 0, cons_pack.GetDim(5) - 1,
      kb.s, kb.e, jb.s, jb.e, ib.s, ib.e + 1,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = prim_pack.GetCoords(b);
        auto &cons = cons_pack(b);
        const auto &prim = prim_pack(b);

        const auto d3B1 =
            ndim > 2 ? (0.5 * (prim(IB1, k + 1, j, i - 1) + prim(IB1, k + 1, j, i)) -
                        0.5 * (prim(IB1, k - 1, j, i - 1) + prim(IB1, k - 1, j, i))) /
                           (coords.Xf<3, 1>(k + 1, j, i) - coords.Xf<3, 1>(k - 1, j, i))
                     : 0.0;
        const auto d1B3 =
            (prim(IB3, k, j, i) - prim(IB3, k, j, i - 1)) / coords.Dxc<1>(k, j, i);
        const auto j2 = d3B1 - d1B3;

        const auto d1B2 =
            (prim(IB2, k, j, i) - prim(IB2, k, j, i - 1)) / coords.Dxc<1>(k, j, i);
        const auto d2B1 =
            ndim > 1 ? (0.5 * (prim(IB1, k, j + 1, i - 1) + prim(IB1, k, j + 1, i)) -
                        0.5 * (prim(IB1, k, j - 1, i - 1) + prim(IB1, k, j - 1, i))) /
                           (coords.Xf<2, 1>(k, j + 1, i) - coords.Xf<2, 1>(k, j - 1, i))
                     : 0.0;
        const auto j3 = d1B2 - d2B1;

        const auto d2B3 =
            ndim > 1 ? (0.5 * (prim(IB3, k, j + 1, i - 1) + prim(IB3, k, j + 1, i)) -
                        0.5 * (prim(IB3, k, j - 1, i - 1) + prim(IB3, k, j - 1, i))) /
                           (coords.Xf<2, 1>(k, j + 1, i) - coords.Xf<2, 1>(k, j - 1, i))
                     : 0.0;
        const auto d3B2 =
            ndim > 2 ? (0.5 * (prim(IB2, k + 1, j, i - 1) + prim(IB2, k + 1, j, i)) -
                        0.5 * (prim(IB2, k - 1, j, i - 1) + prim(IB2, k - 1, j, i))) /
                           (coords.Xf<3, 1>(k + 1, j, i) - coords.Xf<3, 1>(k - 1, j, i))
                     : 0.0;
        const auto j1 = d2B3 - d3B2;

        const Real b1 = 0.5 * (prim(IB1, k, j, i - 1) + prim(IB1, k, j, i));
        const Real b2 = 0.5 * (prim(IB2, k, j, i - 1) + prim(IB2, k, j, i));
        const Real b3 = 0.5 * (prim(IB3, k, j, i - 1) + prim(IB3, k, j, i));
        const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
        // B11: `need_eta_h` (not `eta_h_on`) gates the evaluation, so the per-cell floor
        // gets a real eta_H even on the floor-only call of the mixed rkl2 mode. The
        // DISPERSIVE coefficient is still gated on eta_h_on, so that mode's split is
        // preserved exactly.
        Real eta_h = 0.0;
        if (need_eta_h && use_cache) {
          eta_h = 0.5 * (eta_pack(b, NonidealEtaIdx::H, k, j, i - 1) +
                         eta_pack(b, NonidealEtaIdx::H, k, j, i));
        } else if (need_eta_h) {
          const Real rho = 0.5 * (prim(IDN, k, j, i - 1) + prim(IDN, k, j, i));
          const Real prs = 0.5 * (prim(IPR, k, j, i - 1) + prim(IPR, k, j, i));
          const int i_xe = hall_diff.XeIndex();
          const Real xe = (i_xe >= 0)
                              ? 0.5 * (prim(i_xe, k, j, i - 1) + prim(i_xe, k, j, i))
                              : -1.0;
          eta_h = hall_diff.Get(bmag, rho, prs / rho, xe);
        }
        const Real eta_floor_c =
            (floor_ratio > 0.0) ? hall_diff.EffectiveOhmicFloor(eta_h) : eta_floor_abs;

        Real e1, e2, e3;
        HallEMF(eta_h_on ? eta_h : 0.0, eta_floor_c, bmag, j1, j2, j3, b1, b2, b3, e1, e2,
                e3);

        if (!ct_induction) {
          cons.flux(X1DIR, IB2, k, j, i) += -e3;
          cons.flux(X1DIR, IB3, k, j, i) += e2;
        }
        cons.flux(X1DIR, IEN, k, j, i) += e2 * b3 - e3 * b2;
      });

  if (ndim < 2) {
    return;
  }

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Hall X2 fluxes", DevExecSpace(), 0, cons_pack.GetDim(5) - 1,
      kb.s, kb.e, jb.s, jb.e + 1, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = prim_pack.GetCoords(b);
        auto &cons = cons_pack(b);
        const auto &prim = prim_pack(b);

        const auto d1B2 = (0.5 * (prim(IB2, k, j - 1, i + 1) + prim(IB2, k, j, i + 1)) -
                           0.5 * (prim(IB2, k, j - 1, i - 1) + prim(IB2, k, j, i - 1))) /
                          (coords.Xf<1, 2>(k, j, i + 1) - coords.Xf<1, 2>(k, j, i - 1));
        const auto d2B1 =
            (prim(IB1, k, j, i) - prim(IB1, k, j - 1, i)) / coords.Dxc<2>(k, j, i);
        const auto j3 = d1B2 - d2B1;

        const auto d2B3 =
            (prim(IB3, k, j, i) - prim(IB3, k, j - 1, i)) / coords.Dxc<2>(k, j, i);
        const auto d3B2 =
            ndim > 2 ? (0.5 * (prim(IB2, k + 1, j - 1, i) + prim(IB2, k + 1, j, i)) -
                        0.5 * (prim(IB2, k - 1, j - 1, i) + prim(IB2, k - 1, j, i))) /
                           (coords.Xf<3, 2>(k + 1, j, i) - coords.Xf<3, 2>(k - 1, j, i))
                     : 0.0;
        const auto j1 = d2B3 - d3B2;

        const auto d3B1 =
            ndim > 2 ? (0.5 * (prim(IB1, k + 1, j - 1, i) + prim(IB1, k + 1, j, i)) -
                        0.5 * (prim(IB1, k - 1, j - 1, i) + prim(IB1, k - 1, j, i))) /
                           (coords.Xf<3, 2>(k + 1, j, i) - coords.Xf<3, 2>(k - 1, j, i))
                     : 0.0;
        const auto d1B3 = (0.5 * (prim(IB3, k, j - 1, i + 1) + prim(IB3, k, j, i + 1)) -
                           0.5 * (prim(IB3, k, j - 1, i - 1) + prim(IB3, k, j, i - 1))) /
                          (coords.Xf<1, 2>(k, j, i + 1) - coords.Xf<1, 2>(k, j, i - 1));
        const auto j2 = d3B1 - d1B3;

        const Real b1 = 0.5 * (prim(IB1, k, j - 1, i) + prim(IB1, k, j, i));
        const Real b2 = 0.5 * (prim(IB2, k, j - 1, i) + prim(IB2, k, j, i));
        const Real b3 = 0.5 * (prim(IB3, k, j - 1, i) + prim(IB3, k, j, i));
        const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
        // B11: `need_eta_h` (not `eta_h_on`) gates the evaluation, so the per-cell floor
        // gets a real eta_H even on the floor-only call of the mixed rkl2 mode. The
        // DISPERSIVE coefficient is still gated on eta_h_on, so that mode's split is
        // preserved exactly.
        Real eta_h = 0.0;
        if (need_eta_h && use_cache) {
          eta_h = 0.5 * (eta_pack(b, NonidealEtaIdx::H, k, j - 1, i) +
                         eta_pack(b, NonidealEtaIdx::H, k, j, i));
        } else if (need_eta_h) {
          const Real rho = 0.5 * (prim(IDN, k, j - 1, i) + prim(IDN, k, j, i));
          const Real prs = 0.5 * (prim(IPR, k, j - 1, i) + prim(IPR, k, j, i));
          const int i_xe = hall_diff.XeIndex();
          const Real xe = (i_xe >= 0)
                              ? 0.5 * (prim(i_xe, k, j - 1, i) + prim(i_xe, k, j, i))
                              : -1.0;
          eta_h = hall_diff.Get(bmag, rho, prs / rho, xe);
        }
        const Real eta_floor_c =
            (floor_ratio > 0.0) ? hall_diff.EffectiveOhmicFloor(eta_h) : eta_floor_abs;

        Real e1, e2, e3;
        HallEMF(eta_h_on ? eta_h : 0.0, eta_floor_c, bmag, j1, j2, j3, b1, b2, b3, e1, e2,
                e3);

        if (!ct_induction) {
          cons.flux(X2DIR, IB1, k, j, i) += e3;
          cons.flux(X2DIR, IB3, k, j, i) += -e1;
        }
        cons.flux(X2DIR, IEN, k, j, i) += e3 * b1 - e1 * b3;
      });

  if (ndim < 3) {
    return;
  }

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Hall X3 fluxes", DevExecSpace(), 0, cons_pack.GetDim(5) - 1,
      kb.s, kb.e + 1, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = prim_pack.GetCoords(b);
        auto &cons = cons_pack(b);
        const auto &prim = prim_pack(b);

        const auto d2B3 = (0.5 * (prim(IB3, k - 1, j + 1, i) + prim(IB3, k, j + 1, i)) -
                           0.5 * (prim(IB3, k - 1, j - 1, i) + prim(IB3, k, j - 1, i))) /
                          (coords.Xf<2, 3>(k, j + 1, i) - coords.Xf<2, 3>(k, j - 1, i));
        const auto d3B2 =
            (prim(IB2, k, j, i) - prim(IB2, k - 1, j, i)) / coords.Dxc<3>(k, j, i);
        const auto j1 = d2B3 - d3B2;

        const auto d3B1 =
            (prim(IB1, k, j, i) - prim(IB1, k - 1, j, i)) / coords.Dxc<3>(k, j, i);
        const auto d1B3 = (0.5 * (prim(IB3, k - 1, j, i + 1) + prim(IB3, k, j, i + 1)) -
                           0.5 * (prim(IB3, k - 1, j, i - 1) + prim(IB3, k, j, i - 1))) /
                          (coords.Xf<1, 3>(k, j, i + 1) - coords.Xf<1, 3>(k, j, i - 1));
        const auto j2 = d3B1 - d1B3;

        const auto d1B2 = (0.5 * (prim(IB2, k - 1, j, i + 1) + prim(IB2, k, j, i + 1)) -
                           0.5 * (prim(IB2, k - 1, j, i - 1) + prim(IB2, k, j, i - 1))) /
                          (coords.Xf<1, 3>(k, j, i + 1) - coords.Xf<1, 3>(k, j, i - 1));
        const auto d2B1 = (0.5 * (prim(IB1, k - 1, j + 1, i) + prim(IB1, k, j + 1, i)) -
                           0.5 * (prim(IB1, k - 1, j - 1, i) + prim(IB1, k, j - 1, i))) /
                          (coords.Xf<2, 3>(k, j + 1, i) - coords.Xf<2, 3>(k, j - 1, i));
        const auto j3 = d1B2 - d2B1;

        const Real b1 = 0.5 * (prim(IB1, k - 1, j, i) + prim(IB1, k, j, i));
        const Real b2 = 0.5 * (prim(IB2, k - 1, j, i) + prim(IB2, k, j, i));
        const Real b3 = 0.5 * (prim(IB3, k - 1, j, i) + prim(IB3, k, j, i));
        const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
        // B11: `need_eta_h` (not `eta_h_on`) gates the evaluation, so the per-cell floor
        // gets a real eta_H even on the floor-only call of the mixed rkl2 mode. The
        // DISPERSIVE coefficient is still gated on eta_h_on, so that mode's split is
        // preserved exactly.
        Real eta_h = 0.0;
        if (need_eta_h && use_cache) {
          eta_h = 0.5 * (eta_pack(b, NonidealEtaIdx::H, k - 1, j, i) +
                         eta_pack(b, NonidealEtaIdx::H, k, j, i));
        } else if (need_eta_h) {
          const Real rho = 0.5 * (prim(IDN, k - 1, j, i) + prim(IDN, k, j, i));
          const Real prs = 0.5 * (prim(IPR, k - 1, j, i) + prim(IPR, k, j, i));
          const int i_xe = hall_diff.XeIndex();
          const Real xe = (i_xe >= 0)
                              ? 0.5 * (prim(i_xe, k - 1, j, i) + prim(i_xe, k, j, i))
                              : -1.0;
          eta_h = hall_diff.Get(bmag, rho, prs / rho, xe);
        }
        const Real eta_floor_c =
            (floor_ratio > 0.0) ? hall_diff.EffectiveOhmicFloor(eta_h) : eta_floor_abs;

        Real e1, e2, e3;
        HallEMF(eta_h_on ? eta_h : 0.0, eta_floor_c, bmag, j1, j2, j3, b1, b2, b3, e1, e2,
                e3);

        if (!ct_induction) {
          cons.flux(X3DIR, IB1, k, j, i) += -e2;
          cons.flux(X3DIR, IB2, k, j, i) += e1;
        }
        cons.flux(X3DIR, IEN, k, j, i) += e1 * b2 - e2 * b1;
      });
}
