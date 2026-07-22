//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2024, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file resistivity.cpp
//! \brief

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

// NOTE: OhmicDiffusivity::Get bodies moved to diffusion.hpp -- they are
// KOKKOS_INLINE_FUNCTION and are called from multiple translation units
// (flux kernels, dt estimators, PrecomputeNonidealEta).

Real EstimateResistivityTimestep(MeshData<Real> *md) {
  // get to package via first block in Meshdata (which exists by construction)
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real min_dt_resist = std::numeric_limits<Real>::max();
  const auto ndim = prim_pack.GetNdim();

  Real fac = 0.5;
  if (ndim == 2) {
    fac = 0.25;
  } else if (ndim == 3) {
    fac = 1.0 / 6.0;
  }

  const auto &ohm_diff = hydro_pkg->Param<OhmicDiffusivity>("ohm_diff");

  // Per-cell estimate so it works for both fixed and (spatially varying) ionization
  // coefficients. eta_O is B-independent, but Get() takes bmag/temp for the unified API.
  Kokkos::parallel_reduce(
      "EstimateResistivityTimestep",
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
        const int i_xe = ohm_diff.XeIndex();
        const Real xe = (i_xe >= 0) ? prim(i_xe, k, j, i) : -1.0;
        const auto eta = ohm_diff.Get(bmag, prim(IDN, k, j, i), temp, xe);
        min_dt = fmin(min_dt, SQR(coords.Dxc<1>(k, j, i)) / (eta + TINY_NUMBER));
        if (ndim >= 2) {
          min_dt = fmin(min_dt, SQR(coords.Dxc<2>(k, j, i)) / (eta + TINY_NUMBER));
        }
        if (ndim >= 3) {
          min_dt = fmin(min_dt, SQR(coords.Dxc<3>(k, j, i)) / (eta + TINY_NUMBER));
        }
      },
      Kokkos::Min<Real>(min_dt_resist));

  const auto &cfl_diff = hydro_pkg->Param<Real>("cfl_diff");
  return cfl_diff * fac * min_dt_resist;
}

//---------------------------------------------------------------------------------------
//! Calculate isotropic resistivity with fixed coefficient

void OhmicDiffFluxIsoFixed(MeshData<Real> *md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  // Pack "cons" by NAME: the {Independent} flag also matches rad.Er/Fr and grav.phi, so
  // the fixed IB*/IEN flux indices below would silently target the wrong component if a
  // non-cons Independent field ever sorted ahead of cons (same trap class as the v6 STS
  // bug). The named pack pins the component indices to cons exactly; bit-identical today.
  const std::vector<std::string> cons_names{"cons"};
  auto cons_pack = md->PackVariablesAndFluxes(cons_names, cons_names);
  auto hydro_pkg = pmb->packages.Get("Hydro");

  // Under Constrained Transport (divergence_control=ct) the magnetic *induction* is
  // advanced from an edge-centered curl of eta*J assembled in CT_AddOhmicEMF (ct.cpp),
  // so we must NOT also deposit the resistive contribution into the cell-centered
  // cons.flux(IB*) here -- GS05 reads those and would double-count it. The resistive
  // *energy* (Poynting/heating) term in cons.flux(IEN) is unaffected and stays on the
  // finite-volume energy flux. GLM path (use_ct=false): deposits unchanged, bit-identical.
  const bool ct_induction = hydro_pkg->Param<bool>("use_ct");

  auto const &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});

  const int ndim = pmb->pmy_mesh->ndim;

  const auto &ohm_diff = hydro_pkg->Param<OhmicDiffusivity>("ohm_diff");
  // eta is computed per-face inside the kernels: for ResistivityCoeff::fixed it is the
  // (uniform) coefficient, for ResistivityCoeff::ionization it varies with the local
  // density/temperature via the self-consistent ionization model.

  // Cell-centered eta cache (PrecomputeNonidealEta): face eta = arithmetic average of
  // the two adjacent cached values. When the cache is off, "prim" is packed as a dummy
  // to keep the lambda capture valid; it is never indexed on that branch.
  const bool use_cache = hydro_pkg->Param<bool>("nonideal_eta_cache");
  const auto eta_pack = md->PackVariables(
      std::vector<std::string>{use_cache ? "nonideal_eta" : "prim"});

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Resist. X1 fluxes (ohmic)", DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e + 1,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = prim_pack.GetCoords(b);
        auto &cons = cons_pack(b);
        const auto &prim = prim_pack(b);

        // Face centered current densities
        // j2 = d3B1 - d1B3
        const auto d3B1 =
            ndim > 2 ? (0.5 * (prim(IB1, k + 1, j, i - 1) + prim(IB1, k + 1, j, i)) -
                        0.5 * (prim(IB1, k - 1, j, i - 1) + prim(IB1, k - 1, j, i))) /
                           (coords.Xf<3, 1>(k + 1, j, i) - coords.Xf<3, 1>(k - 1, j, i))
                     : 0.0;

        const auto d1B3 =
            (prim(IB3, k, j, i) - prim(IB3, k, j, i - 1)) / coords.Dxc<1>(k, j, i);

        const auto j2 = d3B1 - d1B3;

        // j3 = d1B2 - d2B1
        const auto d1B2 =
            (prim(IB2, k, j, i) - prim(IB2, k, j, i - 1)) / coords.Dxc<1>(k, j, i);

        const auto d2B1 =
            ndim > 1 ? (0.5 * (prim(IB1, k, j + 1, i - 1) + prim(IB1, k, j + 1, i)) -
                        0.5 * (prim(IB1, k, j - 1, i - 1) + prim(IB1, k, j - 1, i))) /
                           (coords.Xf<2, 1>(k, j + 1, i) - coords.Xf<2, 1>(k, j - 1, i))
                     : 0.0;

        const auto j3 = d1B2 - d2B1;

        // Diffusivity at the i-1/2 face: cached cell-centered average, or evaluated
        // from the face-averaged state.
        Real eta;
        if (use_cache) {
          eta = 0.5 * (eta_pack(b, NonidealEtaIdx::O, k, j, i - 1) +
                       eta_pack(b, NonidealEtaIdx::O, k, j, i));
        } else {
          const Real bm =
              std::sqrt(SQR(0.5 * (prim(IB1, k, j, i - 1) + prim(IB1, k, j, i))) +
                        SQR(0.5 * (prim(IB2, k, j, i - 1) + prim(IB2, k, j, i))) +
                        SQR(0.5 * (prim(IB3, k, j, i - 1) + prim(IB3, k, j, i))));
          const Real rho = 0.5 * (prim(IDN, k, j, i - 1) + prim(IDN, k, j, i));
          const Real prs = 0.5 * (prim(IPR, k, j, i - 1) + prim(IPR, k, j, i));
          const int i_xe = ohm_diff.XeIndex();
          const Real xe = (i_xe >= 0)
                              ? 0.5 * (prim(i_xe, k, j, i - 1) + prim(i_xe, k, j, i))
                              : -1.0;
          eta = ohm_diff.Get(bm, rho, prs / rho, xe);
        }

        if (!ct_induction) {
          cons.flux(X1DIR, IB2, k, j, i) += -eta * j3;
          cons.flux(X1DIR, IB3, k, j, i) += eta * j2;
        }
        cons.flux(X1DIR, IEN, k, j, i) +=
            0.5 * eta *
            ((prim(IB3, k, j, i - 1) + prim(IB3, k, j, i)) * j2 -
             (prim(IB2, k, j, i - 1) + prim(IB2, k, j, i)) * j3);
      });

  if (ndim < 2) {
    return;
  }

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Resist. X2 fluxes (ohmic)", parthenon::DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e + 1, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = prim_pack.GetCoords(b);
        auto &cons = cons_pack(b);
        const auto &prim = prim_pack(b);

        // Face centered current densities
        // j3 = d1B2 - d2B1
        const auto d1B2 = (0.5 * (prim(IB2, k, j - 1, i + 1) + prim(IB2, k, j, i + 1)) -
                           0.5 * (prim(IB2, k, j - 1, i - 1) + prim(IB2, k, j, i - 1))) /
                          (coords.Xf<1, 2>(k, j, i + 1) - coords.Xf<1, 2>(k, j, i - 1));

        const auto d2B1 =
            (prim(IB1, k, j, i) - prim(IB1, k, j - 1, i)) / coords.Dxc<2>(k, j, i);

        const auto j3 = d1B2 - d2B1;

        // j1 = d2B3 - d3B2
        const auto d2B3 =
            (prim(IB3, k, j, i) - prim(IB3, k, j - 1, i)) / coords.Dxc<2>(k, j, i);

        const auto d3B2 =
            ndim > 2 ? (0.5 * (prim(IB2, k + 1, j - 1, i) + prim(IB2, k + 1, j, i)) -
                        0.5 * (prim(IB2, k - 1, j - 1, i) + prim(IB2, k - 1, j, i))) /
                           (coords.Xf<3, 2>(k + 1, j, i) - coords.Xf<3, 2>(k - 1, j, i))
                     : 0.0;

        const auto j1 = d2B3 - d3B2;

        // Diffusivity at the j-1/2 face: cached cell-centered average, or evaluated
        // from the face-averaged state.
        Real eta;
        if (use_cache) {
          eta = 0.5 * (eta_pack(b, NonidealEtaIdx::O, k, j - 1, i) +
                       eta_pack(b, NonidealEtaIdx::O, k, j, i));
        } else {
          const Real bm =
              std::sqrt(SQR(0.5 * (prim(IB1, k, j - 1, i) + prim(IB1, k, j, i))) +
                        SQR(0.5 * (prim(IB2, k, j - 1, i) + prim(IB2, k, j, i))) +
                        SQR(0.5 * (prim(IB3, k, j - 1, i) + prim(IB3, k, j, i))));
          const Real rho = 0.5 * (prim(IDN, k, j - 1, i) + prim(IDN, k, j, i));
          const Real prs = 0.5 * (prim(IPR, k, j - 1, i) + prim(IPR, k, j, i));
          const int i_xe = ohm_diff.XeIndex();
          const Real xe = (i_xe >= 0)
                              ? 0.5 * (prim(i_xe, k, j - 1, i) + prim(i_xe, k, j, i))
                              : -1.0;
          eta = ohm_diff.Get(bm, rho, prs / rho, xe);
        }

        if (!ct_induction) {
          cons.flux(X2DIR, IB1, k, j, i) += eta * j3;
          cons.flux(X2DIR, IB3, k, j, i) += -eta * j1;
        }
        cons.flux(X2DIR, IEN, k, j, i) +=
            0.5 * eta *
            ((prim(IB1, k, j - 1, i) + prim(IB1, k, j, i)) * j3 -
             (prim(IB3, k, j - 1, i) + prim(IB3, k, j, i)) * j1);
      });

  if (ndim < 3) {
    return;
  }

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Resist. X3 fluxes (ohmic)", parthenon::DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, kb.s, kb.e + 1, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = prim_pack.GetCoords(b);
        auto &cons = cons_pack(b);
        const auto &prim = prim_pack(b);

        // Face centered current densities
        // j1 = d2B3 - d3B2
        const auto d2B3 = (0.5 * (prim(IB3, k - 1, j + 1, i) + prim(IB3, k, j + 1, i)) -
                           0.5 * (prim(IB3, k - 1, j - 1, i) + prim(IB3, k, j - 1, i))) /
                          (coords.Xf<2, 3>(k, j + 1, i) - coords.Xf<2, 3>(k, j - 1, i));

        const auto d3B2 =
            (prim(IB2, k, j, i) - prim(IB2, k - 1, j, i)) / coords.Dxc<3>(k, j, i);

        const auto j1 = d2B3 - d3B2;

        // j2 = d3B1 - d1B3
        const auto d3B1 =
            (prim(IB1, k, j, i) - prim(IB1, k - 1, j, i)) / coords.Dxc<3>(k, j, i);

        const auto d1B3 = (0.5 * (prim(IB3, k - 1, j, i + 1) + prim(IB3, k, j, i + 1)) -
                           0.5 * (prim(IB3, k - 1, j, i - 1) + prim(IB3, k, j, i - 1))) /
                          (coords.Xf<1, 3>(k, j, i + 1) - coords.Xf<1, 3>(k, j, i - 1));

        const auto j2 = d3B1 - d1B3;

        // Diffusivity at the k-1/2 face: cached cell-centered average, or evaluated
        // from the face-averaged state.
        Real eta;
        if (use_cache) {
          eta = 0.5 * (eta_pack(b, NonidealEtaIdx::O, k - 1, j, i) +
                       eta_pack(b, NonidealEtaIdx::O, k, j, i));
        } else {
          const Real bm =
              std::sqrt(SQR(0.5 * (prim(IB1, k - 1, j, i) + prim(IB1, k, j, i))) +
                        SQR(0.5 * (prim(IB2, k - 1, j, i) + prim(IB2, k, j, i))) +
                        SQR(0.5 * (prim(IB3, k - 1, j, i) + prim(IB3, k, j, i))));
          const Real rho = 0.5 * (prim(IDN, k - 1, j, i) + prim(IDN, k, j, i));
          const Real prs = 0.5 * (prim(IPR, k - 1, j, i) + prim(IPR, k, j, i));
          const int i_xe = ohm_diff.XeIndex();
          const Real xe = (i_xe >= 0)
                              ? 0.5 * (prim(i_xe, k - 1, j, i) + prim(i_xe, k, j, i))
                              : -1.0;
          eta = ohm_diff.Get(bm, rho, prs / rho, xe);
        }

        if (!ct_induction) {
          cons.flux(X3DIR, IB1, k, j, i) += -eta * j2;
          cons.flux(X3DIR, IB2, k, j, i) += eta * j1;
        }
        cons.flux(X3DIR, IEN, k, j, i) +=
            0.5 * eta *
            ((prim(IB2, k - 1, j, i) + prim(IB2, k, j, i)) * j1 -
             (prim(IB1, k - 1, j, i) + prim(IB1, k, j, i)) * j2);
      });
}

//---------------------------------------------------------------------------------------
//! TODO(pgrete) Calculate Ohmic diffusion, general case, e.g., with varying (Spitzer)
//! coefficient

void OhmicDiffFluxGeneral(MeshData<Real> *md) { PARTHENON_THROW("Needs impl."); }