//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2021, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file diffusion.hpp
//! \brief

#ifndef HYDRO_DIFFUSION_DIFFUSION_HPP_
#define HYDRO_DIFFUSION_DIFFUSION_HPP_

// Parthenon headers
#include <limits>

#include <parthenon/package.hpp>
#include <utils/error_checking.hpp>

// AthenaPK headers
#include "../../main.hpp"
#include "ionization.hpp"

using namespace parthenon::package::prelude;

namespace limiters {
/*----------------------------------------------------------------------------*/
/* vanleer: van Leer slope limiter
 */

KOKKOS_INLINE_FUNCTION Real vanleer(const Real A, const Real B) {
  if (A * B > 0) {
    return 2.0 * A * B / (A + B);
  } else {
    return 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/* minmod: minmod slope limiter
 */

KOKKOS_INLINE_FUNCTION Real minmod(const Real A, const Real B) {
  if (A * B > 0) {
    if (A > 0) {
      return std::min(A, B);
    } else {
      return std::max(A, B);
    }
  } else {
    return 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/* mc: monotonized central slope limiter
 */

KOKKOS_INLINE_FUNCTION Real mc(const Real A, const Real B) {
  return minmod(2.0 * minmod(A, B), (A + B) / 2.0);
}
/*----------------------------------------------------------------------------*/
/* limiter2 and limiter4: call slope limiters to preserve monotonicity
 */

KOKKOS_INLINE_FUNCTION Real lim2(const Real A, const Real B) {
  /* slope limiter */
  return mc(A, B);
}

KOKKOS_INLINE_FUNCTION Real lim4(const Real A, const Real B, const Real C, const Real D) {
  return lim2(lim2(A, B), lim2(C, D));
}
} // namespace limiters

struct ThermalDiffusivity {
 private:
  Real mbar_, me_, kb_;
  Conduction conduction_;
  ConductionCoeff conduction_coeff_type_;
  // "free" coefficient/prefactor. Value depends on conduction is set in the constructor.
  Real coeff_;

 public:
  KOKKOS_INLINE_FUNCTION
  ThermalDiffusivity(Conduction conduction, ConductionCoeff conduction_coeff_type,
                     Real coeff, Real mbar, Real me, Real kb)
      : conduction_(conduction), conduction_coeff_type_(conduction_coeff_type),
        coeff_(coeff), mbar_(mbar), me_(me), kb_(kb) {}

  KOKKOS_INLINE_FUNCTION
  Real Get(const Real pres, const Real rho) const;

  KOKKOS_INLINE_FUNCTION
  Conduction GetType() const { return conduction_; }

  KOKKOS_INLINE_FUNCTION
  ConductionCoeff GetCoeffType() const { return conduction_coeff_type_; }
};

Real EstimateConductionTimestep(MeshData<Real> *md);

//! Calculate isotropic thermal conduction with fixed coefficient
void ThermalFluxIsoFixed(MeshData<Real> *md);
//! Calculate thermal conduction (general case incl. anisotropic and saturated)
void ThermalFluxGeneral(MeshData<Real> *md);

struct MomentumDiffusivity {
 private:
  Real mbar_, me_, kb_;
  Viscosity viscosity_;
  ViscosityCoeff viscosity_coeff_type_;
  // "free" coefficient/prefactor. Value depends on viscosity set in the constructor.
  Real coeff_;

 public:
  KOKKOS_INLINE_FUNCTION
  MomentumDiffusivity(Viscosity viscosity, ViscosityCoeff viscosity_coeff_type,
                      Real coeff, Real mbar, Real me, Real kb)
      : viscosity_(viscosity), viscosity_coeff_type_(viscosity_coeff_type), coeff_(coeff),
        mbar_(mbar), me_(me), kb_(kb) {}

  KOKKOS_INLINE_FUNCTION
  Real Get(const Real pres, const Real rho) const;

  KOKKOS_INLINE_FUNCTION
  Viscosity GetType() const { return viscosity_; }

  KOKKOS_INLINE_FUNCTION
  ViscosityCoeff GetCoeffType() const { return viscosity_coeff_type_; }
};

Real EstimateViscosityTimestep(MeshData<Real> *md);

//! Calculate isotropic viscosity with fixed coefficient
void MomentumDiffFluxIsoFixed(MeshData<Real> *md);
//! Calculate viscosity (general case incl. anisotropic)
void MomentumDiffFluxGeneral(MeshData<Real> *md);

struct OhmicDiffusivity {
 private:
  Real mbar_, me_, kb_;
  Resistivity resistivity_;
  ResistivityCoeff resistivity_coeff_type_;
  // "free" coefficient/prefactor. Value depends on resistivity set in the constructor.
  Real coeff_;
  // Reduced ionization model used when resistivity_coeff_type_ == ionization[_chem].
  Ionization::IonizationModel ion_;
  // prim-pack component index of the chemistry x_e scalar (ionization_chem only; -1 else).
  int i_xe_;
  // Ceiling on the ionization-model eta_O (code units), applied identically in the flux
  // kernels (via the eta cache), the dt estimators, and the fused evaluator. In the
  // first-core interior eta_O grows without bound as x_e collapses while the field is
  // already fully decoupled there; the cap bounds the parabolic dt cost of resolving a
  // diffusion the dynamics no longer depends on. Disabled by default (huge value).
  Real eta_cap_;

 public:
  KOKKOS_INLINE_FUNCTION
  OhmicDiffusivity(Resistivity resistivity, ResistivityCoeff resistivity_coeff_type,
                   Real coeff, Real mbar, Real me, Real kb,
                   Ionization::IonizationModel ion = Ionization::IonizationModel(),
                   int i_xe = -1, Real eta_cap = std::numeric_limits<Real>::max())
      : resistivity_(resistivity), resistivity_coeff_type_(resistivity_coeff_type),
        coeff_(coeff), mbar_(mbar), me_(me), kb_(kb), ion_(ion), i_xe_(i_xe),
        eta_cap_(eta_cap) {}

  // bmag, temp (= code-unit T = p/rho) are used only for the ionization model;
  // ignored for ResistivityCoeff::fixed.
  KOKKOS_INLINE_FUNCTION
  Real Get(const Real bmag, const Real rho, const Real temp) const;

  // x_e-aware overload: for ResistivityCoeff::ionization_chem uses the supplied chemistry
  // electron abundance in the Wardle tensor; otherwise ignores xe and dispatches as 3-arg.
  KOKKOS_INLINE_FUNCTION
  Real Get(const Real bmag, const Real rho, const Real temp, const Real xe) const;

  KOKKOS_INLINE_FUNCTION
  Resistivity GetType() const { return resistivity_; }

  KOKKOS_INLINE_FUNCTION
  ResistivityCoeff GetCoeffType() const { return resistivity_coeff_type_; }

  // prim-pack index of the x_e scalar (-1 if the chem coupling is not active).
  KOKKOS_INLINE_FUNCTION
  int XeIndex() const { return i_xe_; }

  // The reduced ionization model (shared by all non-ideal terms). Used by the fused
  // non-ideal timestep estimator to evaluate all three eta's from ONE Diffusivities call.
  KOKKOS_INLINE_FUNCTION
  const Ionization::IonizationModel &GetIonModel() const { return ion_; }

  // Ceiling on the ionization-model eta_O (code units; huge = disabled).
  KOKKOS_INLINE_FUNCTION
  Real EtaCap() const { return eta_cap_; }
};

// NOTE: the Get() bodies live here in the header (not in resistivity.cpp) because they
// are KOKKOS_INLINE_FUNCTION (implicitly inline): every translation unit that calls them
// (flux kernels, dt estimators, PrecomputeNonidealEta) must see the definition.
KOKKOS_INLINE_FUNCTION
Real OhmicDiffusivity::Get(const Real bmag, const Real rho, const Real temp) const {
  if (resistivity_coeff_type_ == ResistivityCoeff::fixed) {
    return coeff_;
  } else if (resistivity_coeff_type_ == ResistivityCoeff::ionization ||
             resistivity_coeff_type_ == ResistivityCoeff::ionization_chem) {
    // Self-consistent eta_O from the reduced CR+thermal+grain ionization model
    // (Ohmic = c^2/(4 pi sigma_O)). ionization_chem without a supplied x_e falls back to
    // the equilibrium charge solve (the 4-arg Get supplies the chemistry x_e).
    Real eta_O, eta_H, eta_A;
    Ionization::Diffusivities(ion_, rho, temp, bmag, eta_O, eta_H, eta_A);
    return (eta_O < eta_cap_) ? eta_O : eta_cap_;
  } else if (resistivity_coeff_type_ == ResistivityCoeff::spitzer) {
    PARTHENON_FAIL("needs impl");
  } else {
    PARTHENON_FAIL("Unknown Resistivity coeff");
  }
}

KOKKOS_INLINE_FUNCTION
Real OhmicDiffusivity::Get(const Real bmag, const Real rho, const Real temp,
                           const Real xe) const {
  if (resistivity_coeff_type_ == ResistivityCoeff::ionization_chem) {
    // eta_O from the Wardle tensor built on the chemistry-evolved electron abundance
    // (smooth/monotonic; avoids the high-density non-convergence of the equilibrium solve).
    Real eta_O, eta_H, eta_A;
    Ionization::DiffusivitiesFromXe(ion_, rho, temp, bmag, xe, eta_O, eta_H, eta_A);
    return (eta_O < eta_cap_) ? eta_O : eta_cap_;
  }
  return Get(bmag, rho, temp);
}

Real EstimateResistivityTimestep(MeshData<Real> *md);

//! Calculate isotropic resistivity with fixed coefficient
void OhmicDiffFluxIsoFixed(MeshData<Real> *md);

//! Calculate resistivity (general case incl. Spitzer)
void OhmicDiffFluxGeneral(MeshData<Real> *md);

//----------------------------------------------------------------------------------------
//! Ambipolar diffusion. Following Athena++ (ConstDiffusivity), the ambipolar diffusivity
//! is eta_A = coeff * B^2 and the EMF is the perpendicular current,
//! E_A = eta_A (J - (J.bhat) bhat). Parabolic, so it is compatible with both the unsplit
//! and the RKL2 super-time-stepping integrators.
struct AmbipolarDiffusivity {
 private:
  Real mbar_, me_, kb_;
  Ambipolar ambipolar_;
  AmbipolarCoeff ambipolar_coeff_type_;
  // "free" coefficient/prefactor Q_A such that eta_A = Q_A * B^2 (AmbipolarCoeff::fixed).
  Real coeff_;
  // Reduced ionization model used when ambipolar_coeff_type_ == ionization[_chem].
  Ionization::IonizationModel ion_;
  // prim-pack component index of the chemistry x_e scalar (ionization_chem only; -1 else).
  int i_xe_;

 public:
  KOKKOS_INLINE_FUNCTION
  AmbipolarDiffusivity(Ambipolar ambipolar, AmbipolarCoeff ambipolar_coeff_type, Real coeff,
                       Real mbar, Real me, Real kb,
                       Ionization::IonizationModel ion = Ionization::IonizationModel(),
                       int i_xe = -1)
      : ambipolar_(ambipolar), ambipolar_coeff_type_(ambipolar_coeff_type), coeff_(coeff),
        mbar_(mbar), me_(me), kb_(kb), ion_(ion), i_xe_(i_xe) {}

  // temp = code-unit temperature (p = rho*T); ignored for AmbipolarCoeff::fixed.
  KOKKOS_INLINE_FUNCTION
  Real Get(const Real bmag, const Real rho, const Real temp) const;

  // x_e-aware overload: for AmbipolarCoeff::ionization_chem uses the supplied chemistry
  // electron abundance; otherwise ignores xe and dispatches as the 3-arg Get.
  KOKKOS_INLINE_FUNCTION
  Real Get(const Real bmag, const Real rho, const Real temp, const Real xe) const;

  KOKKOS_INLINE_FUNCTION
  Ambipolar GetType() const { return ambipolar_; }

  KOKKOS_INLINE_FUNCTION
  AmbipolarCoeff GetCoeffType() const { return ambipolar_coeff_type_; }

  // prim-pack index of the x_e scalar (-1 if the chem coupling is not active).
  KOKKOS_INLINE_FUNCTION
  int XeIndex() const { return i_xe_; }

  KOKKOS_INLINE_FUNCTION
  const Ionization::IonizationModel &GetIonModel() const { return ion_; }
};

// Get() bodies in the header for cross-TU inlining (see note above OhmicDiffusivity::Get).
KOKKOS_INLINE_FUNCTION
Real AmbipolarDiffusivity::Get(const Real bmag, const Real rho, const Real temp) const {
  if (ambipolar_coeff_type_ == AmbipolarCoeff::fixed) {
    // Athena++ convention: eta_A = coeff * B^2
    return coeff_ * SQR(bmag);
  } else if (ambipolar_coeff_type_ == AmbipolarCoeff::ionization) {
    // Self-consistent eta_A from the reduced CR+thermal+grain ionization model.
    Real eta_O, eta_H, eta_A;
    Ionization::Diffusivities(ion_, rho, temp, bmag, eta_O, eta_H, eta_A);
    return eta_A;
  } else if (ambipolar_coeff_type_ == AmbipolarCoeff::ionization_chem) {
    // chem x_e is supplied via the 4-arg Get; without it fall back to equilibrium.
    Real eta_O, eta_H, eta_A;
    Ionization::Diffusivities(ion_, rho, temp, bmag, eta_O, eta_H, eta_A);
    return eta_A;
  } else {
    PARTHENON_FAIL("Unknown Ambipolar coeff");
  }
}

KOKKOS_INLINE_FUNCTION
Real AmbipolarDiffusivity::Get(const Real bmag, const Real rho, const Real temp,
                               const Real xe) const {
  if (ambipolar_coeff_type_ == AmbipolarCoeff::ionization_chem) {
    // Single-fluid eta_A using the chemistry-evolved electron abundance for rho_i.
    const Real eta_chem = Ionization::AmbipolarEtaFromXe(ion_, rho, bmag, xe);
    // Cap at the equilibrium NICIL eta_A. In dense gas the chemistry x_e collapses
    // (C->CO ionization minimum, x_e ~ 1e-13) which makes the single-fluid
    // eta_A = B^2/(4 pi gamma rho_i rho) blow up and drives the explicit dt to ~0
    // (the dt-collapse that froze mhd4_ad_ohm_rt_chem). The equilibrium tensor model
    // (which ran stably in mhd4_full_apk at dt~1e-6) is the validated diffusivity
    // scale, so we use it as a self-calibrating ceiling: chemistry sets eta_A where
    // it is smaller/comparable, the equilibrium value bounds the dense-gas runaway.
    Real eta_O, eta_H, eta_eq;
    Ionization::Diffusivities(ion_, rho, temp, bmag, eta_O, eta_H, eta_eq);
    return (eta_chem < eta_eq) ? eta_chem : eta_eq;
  }
  return Get(bmag, rho, temp);
}

Real EstimateAmbipolarTimestep(MeshData<Real> *md);

//! Calculate ambipolar diffusion with fixed coefficient
void AmbipolarDiffFluxIsoFixed(MeshData<Real> *md);

//----------------------------------------------------------------------------------------
//! Hall effect. Following Athena++ (ConstDiffusivity), the Hall diffusivity is
//! eta_H = coeff * B / rho and the (non-dissipative) Hall EMF is E_H = eta_H (J x B) / |B|.
//! Hall is dispersive (whistler waves) and is therefore never super-time-stepped: with
//! diffusion/integrator=unsplit everything is unsplit; with rkl2 the parabolic terms go
//! into RKL2 while the Hall EMF stays unsplit under its own (whistler) strict dt (mixed
//! mode). An optional Ohmic resistivity floor can be added inside the kernel for
//! numerical stabilization; in mixed mode diffusion/hall_floor_integrator selects whether
//! the floor stays unsplit with Hall ("unsplit", default) or joins the RKL2 parabolic set
//! ("rkl2", lifting its dx^2 dt ceiling).
struct HallDiffusivity {
 private:
  Real mbar_, me_, kb_;
  Hall hall_;
  HallCoeff hall_coeff_type_;
  // "free" coefficient/prefactor Q_H such that eta_H = Q_H * B / rho.
  Real coeff_;
  // Optional Ohmic resistivity floor (code units) added for stabilization.
  Real ohmic_floor_;
  // Reduced ionization model used when hall_coeff_type_ == ionization[_chem].
  Ionization::IonizationModel ion_;
  // prim-pack component index of the chemistry x_e scalar (ionization_chem only; -1 else).
  int i_xe_;
  // Ceiling on |eta_H| for the IONIZATION-family coefficients (code units), the Hall
  // analogue of OhmicDiffusivity::eta_cap_ and applied identically everywhere eta_H is
  // consumed (flux kernels via the eta cache, dt estimators, fused evaluator). Where
  // eta_O sits at ITS cap the region is Ohmic-dominated and whistlers are diffusively
  // damped faster than Hall acts, while the uncapped |eta_H| ~ B/n_e runaway drives the
  // strict whistler dt (~dx^2/|eta_H|) to zero. eta_H is SIGNED: the cap clamps the
  // magnitude and preserves the sign. Disabled by default (huge value).
  Real eta_cap_;

 public:
  KOKKOS_INLINE_FUNCTION
  HallDiffusivity(Hall hall, HallCoeff hall_coeff_type, Real coeff, Real ohmic_floor,
                  Real mbar, Real me, Real kb,
                  Ionization::IonizationModel ion = Ionization::IonizationModel(),
                  int i_xe = -1, Real eta_cap = std::numeric_limits<Real>::max())
      : hall_(hall), hall_coeff_type_(hall_coeff_type), coeff_(coeff),
        ohmic_floor_(ohmic_floor), mbar_(mbar), me_(me), kb_(kb), ion_(ion), i_xe_(i_xe),
        eta_cap_(eta_cap) {}

  // temp (= code-unit T = p/rho) is used only for the ionization model; ignored for
  // HallCoeff::fixed.
  KOKKOS_INLINE_FUNCTION
  Real Get(const Real bmag, const Real rho, const Real temp) const;

  // x_e-aware overload: for HallCoeff::ionization_chem uses the supplied chemistry electron
  // abundance in the (signed) Wardle-tensor Hall; otherwise dispatches as the 3-arg Get.
  KOKKOS_INLINE_FUNCTION
  Real Get(const Real bmag, const Real rho, const Real temp, const Real xe) const;

  KOKKOS_INLINE_FUNCTION
  Real GetOhmicFloor() const { return ohmic_floor_; }

  KOKKOS_INLINE_FUNCTION
  Hall GetType() const { return hall_; }

  KOKKOS_INLINE_FUNCTION
  HallCoeff GetCoeffType() const { return hall_coeff_type_; }

  // prim-pack index of the x_e scalar (-1 if the chem coupling is not active).
  KOKKOS_INLINE_FUNCTION
  int XeIndex() const { return i_xe_; }

  KOKKOS_INLINE_FUNCTION
  const Ionization::IonizationModel &GetIonModel() const { return ion_; }

  // Ceiling on |eta_H| for the ionization-model coefficients (code units; huge = disabled).
  KOKKOS_INLINE_FUNCTION
  Real EtaCap() const { return eta_cap_; }
};

// Get() bodies in the header for cross-TU inlining (see note above OhmicDiffusivity::Get).
KOKKOS_INLINE_FUNCTION
Real HallDiffusivity::Get(const Real bmag, const Real rho, const Real temp) const {
  if (hall_coeff_type_ == HallCoeff::fixed) {
    // Athena++ convention: eta_H = coeff * B / rho
    return coeff_ * bmag / (rho + TINY_NUMBER);
  } else if (hall_coeff_type_ == HallCoeff::ionization ||
             hall_coeff_type_ == HallCoeff::ionization_chem) {
    // Self-consistent (signed) eta_H from the reduced CR+thermal+grain ionization model;
    // the conductivity-tensor form caps the runaway of the single-fluid eta_H ~ B/n_e at
    // the ionization minimum and changes sign where grains/ions carry the current.
    // ionization_chem without a supplied x_e falls back to the equilibrium charge solve.
    Real eta_O, eta_H, eta_A;
    Ionization::Diffusivities(ion_, rho, temp, bmag, eta_O, eta_H, eta_A);
    // clamp |eta_H| at the cap, preserving the sign
    return (eta_H > eta_cap_) ? eta_cap_ : ((eta_H < -eta_cap_) ? -eta_cap_ : eta_H);
  } else {
    PARTHENON_FAIL("Unknown Hall coeff");
  }
}

KOKKOS_INLINE_FUNCTION
Real HallDiffusivity::Get(const Real bmag, const Real rho, const Real temp,
                          const Real xe) const {
  if (hall_coeff_type_ == HallCoeff::ionization_chem) {
    // Signed eta_H from the Wardle tensor built on the chemistry-evolved x_e.
    Real eta_O, eta_H, eta_A;
    Ionization::DiffusivitiesFromXe(ion_, rho, temp, bmag, xe, eta_O, eta_H, eta_A);
    // clamp |eta_H| at the cap, preserving the sign
    return (eta_H > eta_cap_) ? eta_cap_ : ((eta_H < -eta_cap_) ? -eta_cap_ : eta_H);
  }
  return Get(bmag, rho, temp);
}

//! Hall dt constraint. whistler_on selects the dispersive whistler limit
//! (fac_whistler*dx^2/|eta_H|), floor_on the parabolic limit of the Ohmic stabilizer
//! (fac_par*dx^2/eta_floor). The mixed rkl2+Hall mode needs them separately: the whistler
//! limit stays a strict step constraint while the floor (if integrated by RKL2) only
//! enters the parabolic aggregate that sets the STS stage count.
Real EstimateHallTimestep(MeshData<Real> *md, bool whistler_on = true,
                          bool floor_on = true);

//----------------------------------------------------------------------------------------
//! Per-cell fused evaluation of all three non-ideal diffusivities for the ionization /
//! ionization_chem coefficient family. Each term's per-term Get() recomputes the FULL
//! Pandey-Wardle conductivity tensor (the dominant GPU cost) but keeps only one eta; this
//! helper reproduces the per-term Get() results from the MINIMUM number of tensor solves:
//! one equilibrium Ionization::Diffusivities when any plain-"ionization" term (or the
//! chem-AD equilibrium cap) needs it, plus one Ionization::DiffusivitiesFromXe only when
//! Ohm or Hall are chem-coupled. The chem-coupled AD itself uses the cheap single-fluid
//! AmbipolarEtaFromXe capped at the equilibrium eta_A (matching the 4-arg
//! AmbipolarDiffusivity::Get). For the all-plain-"ionization" configuration this reduces
//! to the single equilibrium solve of the original fused path (bit-identical).
struct FusedNonidealEval {
  Ionization::IonizationModel ion;
  bool ohm_on = false, ad_on = false, hall_on = false;
  bool ohm_chem = false, ad_chem = false, hall_chem = false;
  int i_xe = -1; // prim-pack index of the chemistry x_e scalar (-1 if no chem term)
  // eta_O ceiling mirrored from OhmicDiffusivity::EtaCap() (huge = disabled), so the
  // fused path clamps exactly where the per-term Get() does.
  Real ohm_cap = std::numeric_limits<Real>::max();
  // |eta_H| ceiling mirrored from HallDiffusivity::EtaCap() (signed clamp; huge = disabled).
  Real hall_cap = std::numeric_limits<Real>::max();

  KOKKOS_INLINE_FUNCTION
  void Eta(const Real rho, const Real temp, const Real bmag, const Real xe, Real &eta_O,
           Real &eta_H, Real &eta_A) const {
    Real eO_eq = 0.0, eH_eq = 0.0, eA_eq = 0.0;
    if ((ohm_on && !ohm_chem) || (hall_on && !hall_chem) || ad_on) {
      Ionization::Diffusivities(ion, rho, temp, bmag, eO_eq, eH_eq, eA_eq);
    }
    Real eO_xe = 0.0, eH_xe = 0.0, eA_xe = 0.0;
    if ((ohm_on && ohm_chem) || (hall_on && hall_chem)) {
      Ionization::DiffusivitiesFromXe(ion, rho, temp, bmag, xe, eO_xe, eH_xe, eA_xe);
    }
    eta_O = ohm_on ? (ohm_chem ? eO_xe : eO_eq) : 0.0;
    if (eta_O > ohm_cap) eta_O = ohm_cap;
    eta_H = hall_on ? (hall_chem ? eH_xe : eH_eq) : 0.0;
    if (eta_H > hall_cap) eta_H = hall_cap;
    if (eta_H < -hall_cap) eta_H = -hall_cap;
    if (ad_on) {
      if (ad_chem) {
        // Chemistry x_e single-fluid eta_A capped at the equilibrium value (the
        // dense-gas runaway ceiling; see AmbipolarDiffusivity::Get 4-arg overload).
        const Real eta_chem = Ionization::AmbipolarEtaFromXe(ion, rho, bmag, xe);
        eta_A = (eta_chem < eA_eq) ? eta_chem : eA_eq;
      } else {
        eta_A = eA_eq;
      }
    } else {
      eta_A = 0.0;
    }
  }
};

//----------------------------------------------------------------------------------------
//! Fused non-ideal (Ohmic + ambipolar + Hall) diffusive timestep for the IONIZATION /
//! IONIZATION_CHEM coefficient family. The three separate EstimateX Timestep reductions
//! each recompute the full conductivity tensor but keep only one eta; this routine uses
//! FusedNonidealEval (minimum tensor solves per cell) and reduces the min over whichever
//! of the three terms are active, using each term's own CFL factor. The eta's match the
//! per-term Get() exactly. The parabolic constraint uses the SUM
//! eta_O + hall_ohmic_floor (the flux kernels apply them additively), so with a Hall floor
//! this is slightly stricter (safer) than the separate per-term estimators.
//! Valid when every ACTIVE non-ideal term uses an "ionization"/"ionization_chem" coeff
//! (NOT "fixed"); the caller gates on the "nonideal_dt_fused" package flag.
Real EstimateNonidealTimestepIonizationFused(MeshData<Real> *md);

//! Two-reducer fused estimator for the mixed rkl2+Hall mode: one Wardle tensor
//! evaluation per cell reducing dt_par (parabolic: Ohmic, ambipolar, Hall floor if
//! RKL2-integrated) and dt_strict (whistler + floor if unsplit) separately.
//! floor_strict = (hall_floor_integrator == unsplit). Replaces up to four per-term
//! reduction sweeps per cycle. Uses the EOS-table temperature under eos=hydrogen
//! (consistent with the flux kernels; the per-term estimators use p/rho).
void EstimateNonidealTimestepIonizationFusedMixed(MeshData<Real> *md, bool floor_strict,
                                                  Real &dt_par, Real &dt_strict);

//----------------------------------------------------------------------------------------
//! Component indices of the cell-centered non-ideal diffusivity cache field
//! "nonideal_eta" (registered by the Hydro package when diffusion/eta_cache is on).
namespace NonidealEtaIdx {
constexpr int O = 0; // Ohmic eta_O
constexpr int H = 1; // Hall eta_H (signed)
constexpr int A = 2; // ambipolar eta_A
} // namespace NonidealEtaIdx

//! Fill the cell-centered (eta_O, eta_H, eta_A) cache from prim over the interior
//! expanded by one ghost layer (all cells the face-averaging in the flux kernels
//! touches). Uses the fused FusedNonidealEval (minimum tensor solves per cell) when
//! every active non-ideal term uses an "ionization"/"ionization_chem" coefficient;
//! per-term Get() otherwise (still once per cell instead of once per face per term).
//! Called from CalcDiffFluxes when the "nonideal_eta_cache" package flag is set.
TaskStatus PrecomputeNonidealEta(StateDescriptor *hydro_pkg, MeshData<Real> *md);

//! Calculate the Hall EMF. eta_h_on applies the dispersive term eta_H (JxB)/|B|,
//! floor_on the parabolic Ohmic stabilizer eta_floor*J; both share the same J stencil.
//! (true,true) is the classic unsplit Hall term; the mixed rkl2+Hall mode splits it into
//! an unsplit dispersive part and (optionally) an RKL2-integrated floor part.
void HallDiffFluxIsoFixed(MeshData<Real> *md, bool eta_h_on = true, bool floor_on = true);

//! Which subset of the active diffusion terms CalcDiffFluxes applies.
//!   all        -- every active term (the unsplit integrator).
//!   parabolic  -- conduction, viscosity, Ohmic, ambipolar, and Hall's Ohmic floor if
//!                 diffusion/hall_floor_integrator=rkl2 (the RKL2 stage operator).
//!   dispersive -- the Hall whistler term, plus the floor if it stays unsplit (added to
//!                 the hyperbolic fluxes every stage in the mixed rkl2+Hall mode).
enum class DiffTermSet { all, parabolic, dispersive };

// Calculate diffusion fluxes for the selected term set, i.e., update the .flux views in md.
// refill_eta=false skips the PrecomputeNonidealEta cache refresh so the flux kernels reuse
// the previously cached diffusivities ("diffusion/rkl2_freeze_eta": the RKL2 stage loop
// refreshes eta only at the first stage of each Strang half; RKL2's stability polynomial
// (Meyer+2012) assumes a fixed operator across the super-step, and the resulting
// coefficient lag is O(dt/2) on eta(rho,T,B), which evolves on the collapse timescale).
TaskStatus CalcDiffFluxes(StateDescriptor *hydro_pkg, MeshData<Real> *md,
                          DiffTermSet term_set = DiffTermSet::all, bool refill_eta = true);

#endif //  HYDRO_DIFFUSION_DIFFUSION_HPP_
