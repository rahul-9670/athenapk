//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2024, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ambipolar.cpp
//! \brief Ambipolar diffusion, following the constant-coefficient model of Athena++
//!        (field_diffusion). The ambipolar diffusivity is eta_A = coeff * B^2 and the EMF
//!        is the perpendicular current E_A = eta_A (J - (J.bhat) bhat). The term is added
//!        to the cell-centered (GLM-MHD) magnetic field fluxes together with the
//!        corresponding Poynting flux in the energy equation, exactly as in
//!        OhmicDiffFluxIsoFixed (resistivity.cpp).

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

// NOTE: AmbipolarDiffusivity::Get bodies moved to diffusion.hpp -- they are
// KOKKOS_INLINE_FUNCTION and are called from multiple translation units
// (flux kernels, dt estimators, PrecomputeNonidealEta).

//----------------------------------------------------------------------------------------
//! Perpendicular-current EMF used for ambipolar diffusion:
//!   E = eta * (J - (J.bhat) bhat)
KOKKOS_INLINE_FUNCTION
void PerpCurrentEMF(const Real eta, const Real j1, const Real j2, const Real j3,
                    const Real b1, const Real b2, const Real b3, Real &e1, Real &e2,
                    Real &e3) {
  const Real bsq = SQR(b1) + SQR(b2) + SQR(b3) + TINY_NUMBER;
  const Real jdotb = j1 * b1 + j2 * b2 + j3 * b3;
  e1 = eta * (j1 - jdotb * b1 / bsq);
  e2 = eta * (j2 - jdotb * b2 / bsq);
  e3 = eta * (j3 - jdotb * b3 / bsq);
}

Real EstimateAmbipolarTimestep(MeshData<Real> *md) {
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real min_dt_ad = std::numeric_limits<Real>::max();
  const auto ndim = prim_pack.GetNdim();

  Real fac = 0.5;
  if (ndim == 2) {
    fac = 0.25;
  } else if (ndim == 3) {
    fac = 1.0 / 6.0;
  }

  const auto &ad_diff = hydro_pkg->Param<AmbipolarDiffusivity>("ad_diff");

  Kokkos::parallel_reduce(
      "EstimateAmbipolarTimestep",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExecSpace(), {0, kb.s, jb.s, ib.s},
          {prim_pack.GetDim(5), kb.e + 1, jb.e + 1, ib.e + 1},
          {1, 1, 1, ib.e + 1 - ib.s}),
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &min_dt) {
        const auto &coords = prim_pack.GetCoords(b);
        const auto &prim = prim_pack(b);
        const auto bmag = std::sqrt(SQR(prim(IB1, k, j, i)) + SQR(prim(IB2, k, j, i)) +
                                    SQR(prim(IB3, k, j, i)));
        const auto temp = prim(IPR, k, j, i) / prim(IDN, k, j, i); // code T (p = rho*T)
        const int i_xe = ad_diff.XeIndex();
        const Real xe = (i_xe >= 0) ? prim(i_xe, k, j, i) : -1.0;
        const auto eta = ad_diff.Get(bmag, prim(IDN, k, j, i), temp, xe);
        min_dt = fmin(min_dt, SQR(coords.Dxc<1>(k, j, i)) / (eta + TINY_NUMBER));
        if (ndim >= 2) {
          min_dt = fmin(min_dt, SQR(coords.Dxc<2>(k, j, i)) / (eta + TINY_NUMBER));
        }
        if (ndim >= 3) {
          min_dt = fmin(min_dt, SQR(coords.Dxc<3>(k, j, i)) / (eta + TINY_NUMBER));
        }
      },
      Kokkos::Min<Real>(min_dt_ad));

  const auto &cfl_diff = hydro_pkg->Param<Real>("cfl_diff");
  return cfl_diff * fac * min_dt_ad;
}

//---------------------------------------------------------------------------------------
//! Calculate ambipolar diffusion with fixed coefficient

void AmbipolarDiffFluxIsoFixed(MeshData<Real> *md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  // Pack "cons" by NAME (not {Independent}): pins the fixed IB*/IEN flux indices to cons
  // regardless of registration order (v6-STS-bug trap class). Bit-identical today.
  const std::vector<std::string> cons_names{"cons"};
  auto cons_pack = md->PackVariablesAndFluxes(cons_names, cons_names);
  auto hydro_pkg = pmb->packages.Get("Hydro");

  // Under Constrained Transport (divergence_control=ct) the magnetic *induction* is
  // advanced from an edge-centered curl of the ambipolar perp-current EMF assembled in
  // CT_AddAmbipolarEMF (ct.cpp), so we must NOT also deposit the AD contribution into the
  // cell-centered cons.flux(IB*) here (GS05 reads those and would double-count it). The
  // ambipolar *energy* (Poynting) term in cons.flux(IEN) stays on the finite-volume energy
  // flux. GLM path (use_ct=false): deposits unchanged, bit-identical. Mirrors resistivity.cpp.
  const bool ct_induction = hydro_pkg->Param<bool>("use_ct");

  auto const &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});

  const int ndim = pmb->pmy_mesh->ndim;

  const auto &ad_diff = hydro_pkg->Param<AmbipolarDiffusivity>("ad_diff");

  // Cell-centered eta cache (PrecomputeNonidealEta): face eta = arithmetic average of
  // the two adjacent cached values. When the cache is off, "prim" is packed as a dummy
  // to keep the lambda capture valid; it is never indexed on that branch.
  const bool use_cache = hydro_pkg->Param<bool>("nonideal_eta_cache");
  const auto eta_pack = md->PackVariables(
      std::vector<std::string>{use_cache ? "nonideal_eta" : "prim"});

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Ambipolar X1 fluxes", DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e + 1,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = prim_pack.GetCoords(b);
        auto &cons = cons_pack(b);
        const auto &prim = prim_pack(b);

        // Face-centered current densities at the i-1/2 face.
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

        // Face-averaged magnetic field.
        const Real b1 = 0.5 * (prim(IB1, k, j, i - 1) + prim(IB1, k, j, i));
        const Real b2 = 0.5 * (prim(IB2, k, j, i - 1) + prim(IB2, k, j, i));
        const Real b3 = 0.5 * (prim(IB3, k, j, i - 1) + prim(IB3, k, j, i));
        Real eta;
        if (use_cache) {
          eta = 0.5 * (eta_pack(b, NonidealEtaIdx::A, k, j, i - 1) +
                       eta_pack(b, NonidealEtaIdx::A, k, j, i));
        } else {
          const Real rho = 0.5 * (prim(IDN, k, j, i - 1) + prim(IDN, k, j, i));
          const Real prs = 0.5 * (prim(IPR, k, j, i - 1) + prim(IPR, k, j, i));
          const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
          const int i_xe = ad_diff.XeIndex();
          const Real xe = (i_xe >= 0)
                              ? 0.5 * (prim(i_xe, k, j, i - 1) + prim(i_xe, k, j, i))
                              : -1.0;
          eta = ad_diff.Get(bmag, rho, prs / rho, xe);
        }

        Real e1, e2, e3;
        PerpCurrentEMF(eta, j1, j2, j3, b1, b2, b3, e1, e2, e3);

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
      DEFAULT_LOOP_PATTERN, "Ambipolar X2 fluxes", DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e + 1, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = prim_pack.GetCoords(b);
        auto &cons = cons_pack(b);
        const auto &prim = prim_pack(b);

        // Face-centered current densities at the j-1/2 face.
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

        // Face-averaged magnetic field.
        const Real b1 = 0.5 * (prim(IB1, k, j - 1, i) + prim(IB1, k, j, i));
        const Real b2 = 0.5 * (prim(IB2, k, j - 1, i) + prim(IB2, k, j, i));
        const Real b3 = 0.5 * (prim(IB3, k, j - 1, i) + prim(IB3, k, j, i));
        Real eta;
        if (use_cache) {
          eta = 0.5 * (eta_pack(b, NonidealEtaIdx::A, k, j - 1, i) +
                       eta_pack(b, NonidealEtaIdx::A, k, j, i));
        } else {
          const Real rho = 0.5 * (prim(IDN, k, j - 1, i) + prim(IDN, k, j, i));
          const Real prs = 0.5 * (prim(IPR, k, j - 1, i) + prim(IPR, k, j, i));
          const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
          const int i_xe = ad_diff.XeIndex();
          const Real xe = (i_xe >= 0)
                              ? 0.5 * (prim(i_xe, k, j - 1, i) + prim(i_xe, k, j, i))
                              : -1.0;
          eta = ad_diff.Get(bmag, rho, prs / rho, xe);
        }

        Real e1, e2, e3;
        PerpCurrentEMF(eta, j1, j2, j3, b1, b2, b3, e1, e2, e3);

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
      DEFAULT_LOOP_PATTERN, "Ambipolar X3 fluxes", DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, kb.s, kb.e + 1, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = prim_pack.GetCoords(b);
        auto &cons = cons_pack(b);
        const auto &prim = prim_pack(b);

        // Face-centered current densities at the k-1/2 face.
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

        // Face-averaged magnetic field.
        const Real b1 = 0.5 * (prim(IB1, k - 1, j, i) + prim(IB1, k, j, i));
        const Real b2 = 0.5 * (prim(IB2, k - 1, j, i) + prim(IB2, k, j, i));
        const Real b3 = 0.5 * (prim(IB3, k - 1, j, i) + prim(IB3, k, j, i));
        Real eta;
        if (use_cache) {
          eta = 0.5 * (eta_pack(b, NonidealEtaIdx::A, k - 1, j, i) +
                       eta_pack(b, NonidealEtaIdx::A, k, j, i));
        } else {
          const Real rho = 0.5 * (prim(IDN, k - 1, j, i) + prim(IDN, k, j, i));
          const Real prs = 0.5 * (prim(IPR, k - 1, j, i) + prim(IPR, k, j, i));
          const Real bmag = std::sqrt(SQR(b1) + SQR(b2) + SQR(b3));
          const int i_xe = ad_diff.XeIndex();
          const Real xe = (i_xe >= 0)
                              ? 0.5 * (prim(i_xe, k - 1, j, i) + prim(i_xe, k, j, i))
                              : -1.0;
          eta = ad_diff.Get(bmag, rho, prs / rho, xe);
        }

        Real e1, e2, e3;
        PerpCurrentEMF(eta, j1, j2, j3, b1, b2, b3, e1, e2, e3);

        if (!ct_induction) {
          cons.flux(X3DIR, IB1, k, j, i) += -e2;
          cons.flux(X3DIR, IB2, k, j, i) += e1;
        }
        cons.flux(X3DIR, IEN, k, j, i) += e1 * b2 - e2 * b1;
      });
}
