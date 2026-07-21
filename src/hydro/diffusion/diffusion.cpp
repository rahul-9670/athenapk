//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2021-2023, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file diffusion.cpp
//! \brief

// C++ headers
#include <cmath>
#include <limits>

// Parthenon headers
#include <parthenon/package.hpp>

// AthenaPK headers
#include "../../eos/adiabatic_glmmhd.hpp" // table T for the ionization model under eos=hydrogen
#include "../../main.hpp"
#include "diffusion.hpp"

using namespace parthenon::package::prelude;

namespace {
//! Assemble the per-cell fused non-ideal evaluator (see FusedNonidealEval in
//! diffusion.hpp) from the package's active diffusivity objects, and report whether the
//! coefficient mix is fusable at all: every ACTIVE term on an ionization-family
//! coefficient (plain "ionization" or chem-coupled "ionization_chem"; NOT "fixed").
//! The ionization model is shared by all active terms (identical copies, including the
//! AD closure); the x_e prim index comes from any chem-coupled term (all are built from
//! the same diffusion/xe_scalar_index).
bool MakeFusedNonidealEval(parthenon::StateDescriptor *pkg, FusedNonidealEval &fe) {
  const bool have_ohm = pkg->Param<Resistivity>("resistivity") != Resistivity::none;
  const bool have_ad = pkg->Param<Ambipolar>("ambipolar") != Ambipolar::none;
  const bool have_hall = pkg->Param<Hall>("hall") != Hall::none;
  fe = FusedNonidealEval();
  fe.ohm_on = have_ohm;
  fe.ad_on = have_ad;
  fe.hall_on = have_hall;
  bool fusable = have_ohm || have_ad || have_hall;
  if (have_ohm) {
    const auto &d = pkg->Param<OhmicDiffusivity>("ohm_diff");
    const auto ct = d.GetCoeffType();
    fe.ohm_chem = (ct == ResistivityCoeff::ionization_chem);
    if (ct != ResistivityCoeff::ionization && !fe.ohm_chem) {
      fusable = false;
    } else {
      fe.ion = d.GetIonModel();
      fe.ohm_cap = d.EtaCap();
    }
    if (fe.ohm_chem) fe.i_xe = d.XeIndex();
  }
  if (have_ad) {
    const auto &d = pkg->Param<AmbipolarDiffusivity>("ad_diff");
    const auto ct = d.GetCoeffType();
    fe.ad_chem = (ct == AmbipolarCoeff::ionization_chem);
    if (ct != AmbipolarCoeff::ionization && !fe.ad_chem) {
      fusable = false;
    } else {
      fe.ion = d.GetIonModel();
    }
    if (fe.ad_chem && fe.i_xe < 0) fe.i_xe = d.XeIndex();
  }
  if (have_hall) {
    const auto &d = pkg->Param<HallDiffusivity>("hall_diff");
    const auto ct = d.GetCoeffType();
    fe.hall_chem = (ct == HallCoeff::ionization_chem);
    if (ct != HallCoeff::ionization && !fe.hall_chem) {
      fusable = false;
    } else {
      fe.ion = d.GetIonModel();
      fe.hall_cap = d.EtaCap();
    }
    if (fe.hall_chem && fe.i_xe < 0) fe.i_xe = d.XeIndex();
  }
  return fusable;
}
} // namespace

//----------------------------------------------------------------------------------------
//! Fused non-ideal diffusive timestep (ionization/ionization_chem coeff family). See
//! diffusion.hpp. Evaluates the minimum number of Wardle-tensor solves per cell
//! (FusedNonidealEval) and reduces the diffusive dt over the active Ohmic / ambipolar /
//! Hall terms, replacing one full tensor solve per term with (usually) one total.
Real EstimateNonidealTimestepIonizationFused(MeshData<Real> *md) {
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  const auto ndim = prim_pack.GetNdim();

  const bool have_ohm = hydro_pkg->Param<Resistivity>("resistivity") != Resistivity::none;
  const bool have_ad = hydro_pkg->Param<Ambipolar>("ambipolar") != Ambipolar::none;
  const bool have_hall = hydro_pkg->Param<Hall>("hall") != Hall::none;

  // Per-term CFL factors, matching the individual estimators exactly:
  //   Ohmic/ambipolar (parabolic): 0.5 (1D), 0.25 (2D), 1/6 (3D)
  //   Hall (whistler, dispersive): 0.5 (1D), 1.0 (2D/3D)
  Real fac_par = 0.5;
  if (ndim == 2) {
    fac_par = 0.25;
  } else if (ndim == 3) {
    fac_par = 1.0 / 6.0;
  }
  const Real fac_hall = (ndim == 1) ? 0.5 : 1.0;
  const Real cfl_diff = hydro_pkg->Param<Real>("cfl_diff");
  // Hall's optional Ohmic stabilizer eta_floor*J (applied inside HallEMF) is a real
  // parabolic diffusion that adds to the Ohmic flux, so it enters the parabolic dt
  // constraint summed with eta_O.
  const Real eta_floor_hall =
      have_hall ? hydro_pkg->Param<HallDiffusivity>("hall_diff").GetOhmicFloor() : 0.0;

  // Fused per-cell evaluator (shared ionization model + coefficient-family flags); the
  // caller's "nonideal_dt_fused" gate guarantees the mix is fusable.
  FusedNonidealEval feval;
  MakeFusedNonidealEval(hydro_pkg.get(), feval);
  // Under eos=hydrogen the ionization model's temperature must come from the EOS (mu varies
  // with dissociation/ionization, so T_code != p/rho); ideal path uses p/rho (bit-identical).
  const auto eos_nid = hydro_pkg->Param<AdiabaticGLMMHDEOS>("eos");
  const auto eos_tab_nid = eos_nid.GetEosTable();
  const bool useh2_nid = eos_nid.UseH2Diss();

  Real min_dt = std::numeric_limits<Real>::max();
  Kokkos::parallel_reduce(
      "EstimateNonidealTimestepIonizationFused",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExecSpace(), {0, kb.s, jb.s, ib.s},
          {prim_pack.GetDim(5), kb.e + 1, jb.e + 1, ib.e + 1},
          {1, 1, 1, ib.e + 1 - ib.s}),
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &ldt) {
        const auto &coords = prim_pack.GetCoords(b);
        const auto &prim = prim_pack(b);
        const Real bmag = std::sqrt(SQR(prim(IB1, k, j, i)) + SQR(prim(IB2, k, j, i)) +
                                    SQR(prim(IB3, k, j, i)));
        const Real rho = prim(IDN, k, j, i);
        const Real temp = useh2_nid ? eos_tab_nid.TemperatureKFromPres(
                                          rho, prim(IPR, k, j, i)) /
                                          feval.ion.T_unit
                                    : prim(IPR, k, j, i) / rho; // code T (p = rho*T)
        const Real xe = (feval.i_xe >= 0) ? prim(feval.i_xe, k, j, i) : -1.0;

        // Minimum tensor solves -> all three diffusivities (matches per-term Get()).
        Real eta_O, eta_H, eta_A;
        feval.Eta(rho, temp, bmag, xe, eta_O, eta_H, eta_A);

        // eta is a scalar per cell, so min_d(dx_d^2/eta) = (min_d dx_d^2)/eta exactly.
        Real mindx2 = SQR(coords.Dxc<1>(k, j, i));
        if (ndim >= 2) {
          mindx2 = fmin(mindx2, SQR(coords.Dxc<2>(k, j, i)));
        }
        if (ndim >= 3) {
          mindx2 = fmin(mindx2, SQR(coords.Dxc<3>(k, j, i)));
        }
        // Total Ohmic-like diffusivity: the Ohmic term (if active) plus Hall's Ohmic
        // floor -- the flux kernels apply them additively, so the stability limit is
        // set by their sum.
        const Real eta_O_tot = (have_ohm ? eta_O : 0.0) + eta_floor_hall;
        if (eta_O_tot > 0.0) {
          ldt = fmin(ldt, cfl_diff * fac_par * mindx2 / (eta_O_tot + TINY_NUMBER));
        }
        if (have_ad) {
          ldt = fmin(ldt, cfl_diff * fac_par * mindx2 / (eta_A + TINY_NUMBER));
        }
        if (have_hall) {
          ldt = fmin(ldt, cfl_diff * fac_hall * mindx2 / (std::abs(eta_H) + TINY_NUMBER));
        }
      },
      Kokkos::Min<Real>(min_dt));

  return min_dt;
}

//----------------------------------------------------------------------------------------
//! Two-reducer variant of the fused estimator for the mixed rkl2+Hall mode: one Wardle
//! tensor evaluation per cell feeding TWO separate minima -- dt_par (parabolic: Ohmic,
//! ambipolar, and Hall's Ohmic floor when it is RKL2-integrated) and dt_strict (the
//! dispersive whistler limit, plus the floor's parabolic limit when it stays unsplit).
//! Formulas match the per-term estimators exactly (same eta, same CFL factors), so the
//! resulting dts agree to machine precision.
void EstimateNonidealTimestepIonizationFusedMixed(MeshData<Real> *md,
                                                  const bool floor_strict, Real &dt_par,
                                                  Real &dt_strict) {
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  const auto ndim = prim_pack.GetNdim();

  const bool have_ohm = hydro_pkg->Param<Resistivity>("resistivity") != Resistivity::none;
  const bool have_ad = hydro_pkg->Param<Ambipolar>("ambipolar") != Ambipolar::none;
  const bool have_hall = hydro_pkg->Param<Hall>("hall") != Hall::none;

  Real fac_par = 0.5;
  if (ndim == 2) {
    fac_par = 0.25;
  } else if (ndim == 3) {
    fac_par = 1.0 / 6.0;
  }
  const Real fac_hall = (ndim == 1) ? 0.5 : 1.0;
  const Real cfl_diff = hydro_pkg->Param<Real>("cfl_diff");
  const Real eta_floor_hall =
      have_hall ? hydro_pkg->Param<HallDiffusivity>("hall_diff").GetOhmicFloor() : 0.0;
  // Floor placement (see hall_floor_integrator): strict -> it rides with the unsplit
  // Hall EMF and constrains dt_strict; otherwise it is RKL2-integrated and adds to the
  // parabolic Ohmic diffusivity in dt_par.
  const Real eta_floor_par = floor_strict ? 0.0 : eta_floor_hall;
  const Real eta_floor_strict = floor_strict ? eta_floor_hall : 0.0;

  // Fused per-cell evaluator (shared ionization model + coefficient-family flags); the
  // caller's "nonideal_dt_fused" gate guarantees the mix is fusable.
  FusedNonidealEval feval;
  MakeFusedNonidealEval(hydro_pkg.get(), feval);
  const auto eos_nid = hydro_pkg->Param<AdiabaticGLMMHDEOS>("eos");
  const auto eos_tab_nid = eos_nid.GetEosTable();
  const bool useh2_nid = eos_nid.UseH2Diss();

  Real min_dt_par = std::numeric_limits<Real>::max();
  Real min_dt_strict = std::numeric_limits<Real>::max();
  Kokkos::parallel_reduce(
      "EstimateNonidealTimestepIonizationFusedMixed",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExecSpace(), {0, kb.s, jb.s, ib.s},
          {prim_pack.GetDim(5), kb.e + 1, jb.e + 1, ib.e + 1},
          {1, 1, 1, ib.e + 1 - ib.s}),
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &ldt_par,
                    Real &ldt_strict) {
        const auto &coords = prim_pack.GetCoords(b);
        const auto &prim = prim_pack(b);
        const Real bmag = std::sqrt(SQR(prim(IB1, k, j, i)) + SQR(prim(IB2, k, j, i)) +
                                    SQR(prim(IB3, k, j, i)));
        const Real rho = prim(IDN, k, j, i);
        const Real temp = useh2_nid ? eos_tab_nid.TemperatureKFromPres(
                                          rho, prim(IPR, k, j, i)) /
                                          feval.ion.T_unit
                                    : prim(IPR, k, j, i) / rho; // code T (p = rho*T)
        const Real xe = (feval.i_xe >= 0) ? prim(feval.i_xe, k, j, i) : -1.0;

        Real eta_O, eta_H, eta_A;
        feval.Eta(rho, temp, bmag, xe, eta_O, eta_H, eta_A);

        Real mindx2 = SQR(coords.Dxc<1>(k, j, i));
        if (ndim >= 2) {
          mindx2 = fmin(mindx2, SQR(coords.Dxc<2>(k, j, i)));
        }
        if (ndim >= 3) {
          mindx2 = fmin(mindx2, SQR(coords.Dxc<3>(k, j, i)));
        }
        // Separate per-term minima (NOT the eta_O+floor sum of the unsplit fused
        // estimator), matching the per-term path this replaces --
        // EstimateResistivityTimestep + EstimateHallTimestep(false,true) +
        // EstimateHallTimestep(true,!floor_rkl2). Bit-identical on the ideal-gas path;
        // under eos=hydrogen the per-term estimators evaluate the tensor at temp=p/rho
        // while this kernel (like PrecomputeNonidealEta and the flux kernels) uses the
        // EOS-table temperature -- the consistent choice, but a (small) dt difference
        // to quantify in A/B before switching a production run over.
        if (have_ohm) {
          ldt_par = fmin(ldt_par, cfl_diff * fac_par * mindx2 / (eta_O + TINY_NUMBER));
        }
        if (eta_floor_par > 0.0) {
          ldt_par =
              fmin(ldt_par, cfl_diff * fac_par * mindx2 / (eta_floor_par + TINY_NUMBER));
        }
        if (have_ad) {
          ldt_par = fmin(ldt_par, cfl_diff * fac_par * mindx2 / (eta_A + TINY_NUMBER));
        }
        if (have_hall) {
          ldt_strict =
              fmin(ldt_strict, cfl_diff * fac_hall * mindx2 / (std::abs(eta_H) + TINY_NUMBER));
        }
        if (eta_floor_strict > 0.0) {
          ldt_strict =
              fmin(ldt_strict, cfl_diff * fac_par * mindx2 / (eta_floor_strict + TINY_NUMBER));
        }
      },
      Kokkos::Min<Real>(min_dt_par), Kokkos::Min<Real>(min_dt_strict));

  dt_par = min_dt_par;
  dt_strict = min_dt_strict;
}

//----------------------------------------------------------------------------------------
//! Fill the cell-centered non-ideal diffusivity cache (see diffusion.hpp). Bounds are the
//! interior expanded by one layer in each existing dimension -- exactly the cells whose
//! values the face-averaging in the flux kernels reads. Ghost prim is valid there (the
//! per-face path already read the same cells).
TaskStatus PrecomputeNonidealEta(StateDescriptor *hydro_pkg, MeshData<Real> *md) {
  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});
  auto eta_pack = md->PackVariables(std::vector<std::string>{"nonideal_eta"});

  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  const auto ndim = prim_pack.GetNdim();
  const int il = ib.s - 1, iu = ib.e + 1;
  const int jl = (ndim > 1) ? jb.s - 1 : jb.s, ju = (ndim > 1) ? jb.e + 1 : jb.e;
  const int kl = (ndim > 2) ? kb.s - 1 : kb.s, ku = (ndim > 2) ? kb.e + 1 : kb.e;

  const bool have_ohm = hydro_pkg->Param<Resistivity>("resistivity") != Resistivity::none;
  const bool have_ad = hydro_pkg->Param<Ambipolar>("ambipolar") != Ambipolar::none;
  const bool have_hall = hydro_pkg->Param<Hall>("hall") != Hall::none;

  // Inactive terms get dummy-constructed diffusivities (never called; the value-capture
  // into the device lambda needs *some* object).
  const auto ohm_diff =
      have_ohm ? hydro_pkg->Param<OhmicDiffusivity>("ohm_diff")
               : OhmicDiffusivity(Resistivity::none, ResistivityCoeff::none, 0.0, 0.0,
                                  0.0, 0.0);
  const auto ad_diff =
      have_ad ? hydro_pkg->Param<AmbipolarDiffusivity>("ad_diff")
              : AmbipolarDiffusivity(Ambipolar::none, AmbipolarCoeff::none, 0.0, 0.0, 0.0,
                                     0.0);
  const auto hall_diff = have_hall
                             ? hydro_pkg->Param<HallDiffusivity>("hall_diff")
                             : HallDiffusivity(Hall::none, HallCoeff::none, 0.0, 0.0, 0.0,
                                               0.0, 0.0);

  // Ionization-family fast path: FusedNonidealEval fills all three components from the
  // minimum number of tensor solves per cell -- one equilibrium solve when all terms are
  // plain "ionization" (bit-identical to the per-term Get calls, which each compute the
  // full tensor anyway), plus only the strictly needed chem-x_e evaluations for
  // "ionization_chem" terms (production tier-4 mix: Ohm/Hall equilibrium + chem-capped
  // AD = ONE tensor solve instead of three).
  FusedNonidealEval feval;
  const bool fused = MakeFusedNonidealEval(hydro_pkg, feval);

  const auto eos_nid = hydro_pkg->Param<AdiabaticGLMMHDEOS>("eos");
  const auto eos_tab_nid = eos_nid.GetEosTable();
  const bool useh2_nid = eos_nid.UseH2Diss();
  const auto T_unit_nid = feval.ion.T_unit;

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "PrecomputeNonidealEta", DevExecSpace(), 0,
      prim_pack.GetDim(5) - 1, kl, ku, jl, ju, il, iu,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &prim = prim_pack(b);
        const Real rho = prim(IDN, k, j, i);
        const Real temp = useh2_nid ? eos_tab_nid.TemperatureKFromPres(
                                          rho, prim(IPR, k, j, i)) /
                                          T_unit_nid
                                    : prim(IPR, k, j, i) / rho; // code T (p = rho*T)
        const Real bmag = std::sqrt(SQR(prim(IB1, k, j, i)) + SQR(prim(IB2, k, j, i)) +
                                    SQR(prim(IB3, k, j, i)));
        if (fused) {
          const Real xe = (feval.i_xe >= 0) ? prim(feval.i_xe, k, j, i) : -1.0;
          Real eO, eH, eA;
          feval.Eta(rho, temp, bmag, xe, eO, eH, eA);
          eta_pack(b, NonidealEtaIdx::O, k, j, i) = eO;
          eta_pack(b, NonidealEtaIdx::H, k, j, i) = eH;
          eta_pack(b, NonidealEtaIdx::A, k, j, i) = eA;
        } else {
          Real eO = 0.0, eH = 0.0, eA = 0.0;
          if (have_ohm) {
            const int ix = ohm_diff.XeIndex();
            const Real xe = (ix >= 0) ? prim(ix, k, j, i) : -1.0;
            eO = ohm_diff.Get(bmag, rho, temp, xe);
          }
          if (have_hall) {
            const int ix = hall_diff.XeIndex();
            const Real xe = (ix >= 0) ? prim(ix, k, j, i) : -1.0;
            eH = hall_diff.Get(bmag, rho, temp, xe);
          }
          if (have_ad) {
            const int ix = ad_diff.XeIndex();
            const Real xe = (ix >= 0) ? prim(ix, k, j, i) : -1.0;
            eA = ad_diff.Get(bmag, rho, temp, xe);
          }
          eta_pack(b, NonidealEtaIdx::O, k, j, i) = eO;
          eta_pack(b, NonidealEtaIdx::H, k, j, i) = eH;
          eta_pack(b, NonidealEtaIdx::A, k, j, i) = eA;
        }
      });
  return TaskStatus::complete;
}

TaskStatus CalcDiffFluxes(StateDescriptor *hydro_pkg, MeshData<Real> *md,
                          DiffTermSet term_set, bool refill_eta) {
  // Term-set selection (see DiffTermSet in diffusion.hpp): "all" is the unsplit
  // integrator; the mixed rkl2+Hall mode calls "parabolic" from the RKL2 stages and
  // "dispersive" from the per-stage hyperbolic flux calculation. Hall's Ohmic floor is
  // parabolic physics inside the Hall kernel: hall_floor_int_rkl2 decides which of the
  // two calls applies it.
  const bool do_parabolic = (term_set != DiffTermSet::dispersive);
  const auto &hall = hydro_pkg->Param<Hall>("hall");
  const bool floor_rkl2 =
      (hall != Hall::none) && hydro_pkg->Param<bool>("hall_floor_int_rkl2");
  // Hall term (either part) requested in this set?
  bool do_hall = false;
  bool hall_eta_h_on = false, hall_floor_on = false;
  if (hall != Hall::none) {
    if (term_set == DiffTermSet::all) {
      do_hall = true;
      hall_eta_h_on = true;
      hall_floor_on = true;
    } else if (term_set == DiffTermSet::dispersive) {
      do_hall = true;
      hall_eta_h_on = true;
      hall_floor_on = !floor_rkl2;
    } else { // parabolic: only the floor, and only if it is RKL2-integrated
      do_hall = floor_rkl2;
      hall_floor_on = true;
    }
  }

  if (do_parabolic) {
    const auto &conduction = hydro_pkg->Param<Conduction>("conduction");
    if (conduction != Conduction::none) {
      const auto &thermal_diff = hydro_pkg->Param<ThermalDiffusivity>("thermal_diff");

      if (conduction == Conduction::isotropic &&
          thermal_diff.GetCoeffType() == ConductionCoeff::fixed) {
        ThermalFluxIsoFixed(md);
      } else {
        ThermalFluxGeneral(md);
      }
    }
    const auto &viscosity = hydro_pkg->Param<Viscosity>("viscosity");
    if (viscosity != Viscosity::none) {
      const auto &mom_diff = hydro_pkg->Param<MomentumDiffusivity>("mom_diff");

      if (viscosity == Viscosity::isotropic &&
          mom_diff.GetCoeffType() == ViscosityCoeff::fixed) {
        MomentumDiffFluxIsoFixed(md);
      } else {
        MomentumDiffFluxGeneral(md);
      }
    }
  }
  const auto &resistivity = hydro_pkg->Param<Resistivity>("resistivity");
  const auto &ambipolar = hydro_pkg->Param<Ambipolar>("ambipolar");
  const bool do_ohm = do_parabolic && (resistivity != Resistivity::none);
  const bool do_ad = do_parabolic && (ambipolar != Ambipolar::none);
  // Fill the cell-centered non-ideal diffusivity cache once per call (per stage); the
  // Ohmic/ambipolar/Hall flux kernels below then face-average the cached values instead
  // of re-running the ionization solve per face. Skipped when this term set runs no
  // eta-consuming kernel (e.g. the parabolic set with an unsplit Hall floor), or when
  // the caller freezes eta across RKL2 stages (refill_eta=false; cache still holds the
  // values from the first stage of this Strang half).
  if (refill_eta && hydro_pkg->Param<bool>("nonideal_eta_cache") &&
      (do_ohm || do_ad || (do_hall && hall_eta_h_on))) {
    PrecomputeNonidealEta(hydro_pkg, md);
  }
  if (do_ohm) {
    const auto &ohm_diff = hydro_pkg->Param<OhmicDiffusivity>("ohm_diff");

    // OhmicDiffFluxIsoFixed evaluates the diffusivity per-face via OhmicDiffusivity::Get,
    // so it handles both the fixed and the (spatially varying) ionization coefficient.
    if (resistivity == Resistivity::ohmic &&
        (ohm_diff.GetCoeffType() == ResistivityCoeff::fixed ||
         ohm_diff.GetCoeffType() == ResistivityCoeff::ionization ||
         ohm_diff.GetCoeffType() == ResistivityCoeff::ionization_chem)) {
      OhmicDiffFluxIsoFixed(md);
    } else {
      OhmicDiffFluxGeneral(md);
    }
  }
  if (do_ad) {
    AmbipolarDiffFluxIsoFixed(md);
  }
  if (do_hall) {
    HallDiffFluxIsoFixed(md, hall_eta_h_on, hall_floor_on);
  }
  return TaskStatus::complete;
}
