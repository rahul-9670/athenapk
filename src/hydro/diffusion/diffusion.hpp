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
  // Ceiling on eta_A for the IONIZATION-family coefficients (code units), the ambipolar
  // analogue of OhmicDiffusivity::eta_cap_ / HallDiffusivity::eta_cap_ and applied
  // identically everywhere eta_A is consumed (flux kernels via the eta cache, dt
  // estimators, fused evaluator). The single-fluid/tensor eta_A = B^2/(4 pi gamma rho_i rho)
  // runs away (~1/rho^2) in low-density magnetized cells and drives the explicit/RKL2
  // parabolic dt to ~0 there; this bounds that runaway in the diffusion-decoupled regime
  // exactly as eta_ohm_cap does for Ohm. Default (huge) = disabled (bit-identical).
  Real eta_cap_;

 public:
  KOKKOS_INLINE_FUNCTION
  AmbipolarDiffusivity(Ambipolar ambipolar, AmbipolarCoeff ambipolar_coeff_type, Real coeff,
                       Real mbar, Real me, Real kb,
                       Ionization::IonizationModel ion = Ionization::IonizationModel(),
                       int i_xe = -1, Real eta_cap = std::numeric_limits<Real>::max())
      : ambipolar_(ambipolar), ambipolar_coeff_type_(ambipolar_coeff_type), coeff_(coeff),
        mbar_(mbar), me_(me), kb_(kb), ion_(ion), i_xe_(i_xe), eta_cap_(eta_cap) {}

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

  KOKKOS_INLINE_FUNCTION
  Real EtaCap() const { return eta_cap_; }

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
    // Athena++ convention: eta_A = coeff * B^2 (uncapped: fixed-coeff test/idealised path)
    return coeff_ * SQR(bmag);
  } else if (ambipolar_coeff_type_ == AmbipolarCoeff::ionization) {
    // Self-consistent eta_A from the reduced CR+thermal+grain ionization model.
    Real eta_O, eta_H, eta_A;
    Ionization::Diffusivities(ion_, rho, temp, bmag, eta_O, eta_H, eta_A);
    return (eta_A < eta_cap_) ? eta_A : eta_cap_;
  } else if (ambipolar_coeff_type_ == AmbipolarCoeff::ionization_chem) {
    // chem x_e is supplied via the 4-arg Get; without it fall back to equilibrium.
    Real eta_O, eta_H, eta_A;
    Ionization::Diffusivities(ion_, rho, temp, bmag, eta_O, eta_H, eta_A);
    return (eta_A < eta_cap_) ? eta_A : eta_cap_;
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
    const Real eta_A = (eta_chem < eta_eq) ? eta_chem : eta_eq;
    // Absolute ceiling (diffusion-decoupled runaway guard) ON TOP of the equilibrium
    // ceiling: even eta_eq blows up (~1/rho^2) in low-density magnetized cells.
    return (eta_A < eta_cap_) ? eta_A : eta_cap_;
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
  // B11 (WP-16 part 3): the Ohmic stabilizer expressed as a FRACTION of the local |eta_H|,
  // so the floor tracks the term it stabilizes instead of being a fixed number that merely
  // happens to be large enough. 0 = disabled => EffectiveOhmicFloor() returns ohmic_floor_
  // exactly, and every consumer is bit-identical to the pre-B11 code.
  Real ohmic_floor_ratio_;

 public:
  KOKKOS_INLINE_FUNCTION
  HallDiffusivity(Hall hall, HallCoeff hall_coeff_type, Real coeff, Real ohmic_floor,
                  Real mbar, Real me, Real kb,
                  Ionization::IonizationModel ion = Ionization::IonizationModel(),
                  int i_xe = -1, Real eta_cap = std::numeric_limits<Real>::max(),
                  Real ohmic_floor_ratio = 0.0)
      : hall_(hall), hall_coeff_type_(hall_coeff_type), coeff_(coeff),
        ohmic_floor_(ohmic_floor), mbar_(mbar), me_(me), kb_(kb), ion_(ion), i_xe_(i_xe),
        eta_cap_(eta_cap), ohmic_floor_ratio_(ohmic_floor_ratio) {}

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
  Real GetOhmicFloorRatio() const { return ohmic_floor_ratio_; }

  //! B11: the Ohmic stabilizer ACTUALLY applied in a cell whose Hall diffusivity is eta_h.
  //!
  //! WHY A RATIO. `hall_ohmic_floor_code` alone is an absolute number applied regardless of
  //! the local eta_H, which is safe only while eta_H happens to stay small. Measured on the
  //! analytic Hall eigenmode (src/pgen/diffusion.cpp iprob=60 -- a circularly polarized mode
  //! the non-dissipative Hall term must conserve EXACTLY): in 3D the scheme AMPLIFIES rather
  //! than damps once eta_floor/|eta_H| falls below ~0.11. The onset is sharp -- crossing it
  //! takes the amplitude error from +1.6 % to +2.65e5, and 256^3 crashes outright -- and it
  //! is a RATIO, not a resolution-dependent number: above threshold, 64^3 and 128^3 agree to
  //! ~5 %. That is what theory predicts, since Hall (eta_H k^2) and Ohmic (eta_O k^2) damping
  //! carry the same power of k. 1D is unaffected at every resolution to 1024, which is why
  //! this went unnoticed: the historical Hall validation was 1D-only.
  //!
  //! max(), NOT replace: the absolute floor still applies where it is the larger of the two,
  //! so raising the ratio can only ever ADD dissipation, never remove it. That one-sided
  //! property is what makes the switch safe to enable on a run whose absolute floor already
  //! dominates -- it is then a provable no-op (see docs/validation/WP16b_hall_3d_instability.md).
  //!
  //! NOT FREE: the stabilizer is a real resistivity in the dispersion relation, so a larger
  //! ratio buys stability at the cost of Hall fidelity (measured omega error 6.6e-3 at ratio
  //! 0.125 rising to 4.7e-2 at 1.0). Sit just above threshold, not far above it.
  KOKKOS_INLINE_FUNCTION
  Real EffectiveOhmicFloor(const Real eta_h) const {
    const Real from_ratio = ohmic_floor_ratio_ * (eta_h < 0.0 ? -eta_h : eta_h);
    return (from_ratio > ohmic_floor_) ? from_ratio : ohmic_floor_;
  }

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
  // eta_A ceiling mirrored from AmbipolarDiffusivity::EtaCap() (huge = disabled). Bounds the
  // low-density single-fluid/tensor eta_A runaway that otherwise sets the parabolic dt.
  Real ad_cap = std::numeric_limits<Real>::max();

  // --- WS-4 DUST COUPLING (flagship audit item 5) ---------------------------------------
  // prim-pack indices of the evolved dust scalars (f_dg, a_c) and the reference grain
  // radius the STATIC MRN bins correspond to. dust_on=false (default) => the ionization
  // model keeps its frozen ISM MRN population and every result is bit-identical.
  int i_fdg = -1;      // prim index of the dust-to-gas mass ratio
  int i_ac = -1;       // prim index of the characteristic grain radius a_c [cm]
  Real dust_a_ref = 1.0e-5; // [cm] a_c the static bins represent (DustModel::a_ref)
  Real dust_fdg_ref = 0.01; // dust-to-gas ratio the static bins represent
  bool dust_on = false;

  //! Per-cell ionization model: the static one, or a copy rescaled to the evolved grain
  //! population. The copy is only made when dust coupling is ON, so uncoupled runs pay
  //! nothing (the model is a ~40-Real POD; copying it per cell is not free on GPU).
  KOKKOS_INLINE_FUNCTION
  Ionization::IonizationModel CellIon(const Real f_dg, const Real a_c) const {
    Ionization::IonizationModel mc = ion;
    // Guard degenerate/uninitialized scalars: fall back to the static population rather
    // than divide by zero or hand the charge solve a zero-radius grain.
    //
    // HISTORY (2026-07-30): collapse_be used to leave these scalars at 0.0 even with
    // <dust> evolve=true (its species zero-fill ran past the chemistry block), so both
    // guards fired and this coupling was a silent no-op. FIXED in collapse_be.cpp, which
    // now seeds f_dg = <dust> f_dg_ref and a_c = <dust> a_ref_cm. The guards stay as a
    // backstop for other pgens: falling back to the static population is always safer than
    // handing the charge solve a zero-radius grain. NOTE that a_scale == 1 at t=0 BY
    // DESIGN (the static MRN bins represent exactly the reference population), so this
    // coupling only bites once the dust package moves (f_dg, a_c) away from the reference
    // -- changing <dust> a_ref_cm alone moves the IC and the reference together and does
    // nothing.
    mc.f_dg_scale = (f_dg > 0.0 && dust_fdg_ref > 0.0) ? f_dg / dust_fdg_ref : 1.0;
    mc.a_scale = (a_c > 0.0 && dust_a_ref > 0.0) ? a_c / dust_a_ref : 1.0;
    return mc;
  }

  //! UNCAPPED diffusivities: everything Eta() does except the three ceiling clamps.
  //! Eta() == EtaRaw() + clamps, so the physics path is unchanged by construction; this
  //! exists so the cap diagnostic (FillNonidealCapDiag) can report how much diffusivity
  //! the ceilings actually remove, without a second Wardle tensor solve.
  KOKKOS_INLINE_FUNCTION
  void EtaRaw(const Real rho, const Real temp, const Real bmag, const Real xe, Real &eta_O,
              Real &eta_H, Real &eta_A, const Real f_dg = -1.0,
              const Real a_c = -1.0) const {
    // Dust coupling: rescale the grain population to the evolved (f_dg, a_c) for THIS cell.
    // Pointer (not a ternary) on purpose: `dust_on ? CellIon(...) : ion` is a PRVALUE, so it
    // would copy the model even in the OFF branch -- the exact cost this is avoiding.
    Ionization::IonizationModel mcell;
    const Ionization::IonizationModel *pim = &ion;
    if (dust_on) {
      mcell = CellIon(f_dg, a_c);
      pim = &mcell;
    }
    const Ionization::IonizationModel &im = *pim;
    Real eO_eq = 0.0, eH_eq = 0.0, eA_eq = 0.0;
    if ((ohm_on && !ohm_chem) || (hall_on && !hall_chem) || ad_on) {
      Ionization::Diffusivities(im, rho, temp, bmag, eO_eq, eH_eq, eA_eq);
    }
    Real eO_xe = 0.0, eH_xe = 0.0, eA_xe = 0.0;
    if ((ohm_on && ohm_chem) || (hall_on && hall_chem)) {
      Ionization::DiffusivitiesFromXe(im, rho, temp, bmag, xe, eO_xe, eH_xe, eA_xe);
    }
    eta_O = ohm_on ? (ohm_chem ? eO_xe : eO_eq) : 0.0;
    eta_H = hall_on ? (hall_chem ? eH_xe : eH_eq) : 0.0;
    if (ad_on) {
      if (ad_chem) {
        // NOTE: the equilibrium ceiling min(eta_chem, eA_eq) is PHYSICS (the dense-gas
        // runaway closure), not a numerical ceiling, so it stays in the "raw" value; only
        // the three explicit eta_*_cap_code ceilings are treated as caps by the diagnostic.
        const Real eta_chem = Ionization::AmbipolarEtaFromXe(im, rho, bmag, xe);
        eta_A = (eta_chem < eA_eq) ? eta_chem : eA_eq;
      } else {
        eta_A = eA_eq;
      }
    } else {
      eta_A = 0.0;
    }
  }

  //! Apply the three explicit numerical ceilings in place. Split out of Eta() so the cap
  //! diagnostic applies EXACTLY the same clamps the physics path does (no drift between a
  //! diagnostic copy and the real thing). When a term is off its eta is 0 and its cap is
  //! positive, so clamping unconditionally is a no-op there.
  KOKKOS_INLINE_FUNCTION
  void Clamp(Real &eta_O, Real &eta_H, Real &eta_A) const {
    if (eta_O > ohm_cap) eta_O = ohm_cap;
    if (eta_H > hall_cap) eta_H = hall_cap;
    if (eta_H < -hall_cap) eta_H = -hall_cap;
    // ad_on guard so an off-AD eta_A stays exactly 0 even for a pathological ad_cap <= 0
    // (pre-refactor the eta_A clamp lived inside the ad_on branch).
    if (ad_on && eta_A > ad_cap) eta_A = ad_cap;
  }

  KOKKOS_INLINE_FUNCTION
  void Eta(const Real rho, const Real temp, const Real bmag, const Real xe, Real &eta_O,
           Real &eta_H, Real &eta_A, const Real f_dg = -1.0, const Real a_c = -1.0) const {
    EtaRaw(rho, temp, bmag, xe, eta_O, eta_H, eta_A, f_dg, a_c);
    Clamp(eta_O, eta_H, eta_A);
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

//----------------------------------------------------------------------------------------
//! Component indices of the cap-activation diagnostic field "diff.capdiag" (registered when
//! diffusion/cap_diag=true). FLAGSHIP AUDIT ITEM 1: eta_ohm_cap_code / eta_ad_cap_code /
//! eta_hall_cap_code and the Hall Ohmic floor are numerical stabilizers that MODIFY the
//! induction equation. Without this diagnostic there is no way to tell whether a flux-
//! retention result is physics or a stabilizer. Components 0-2 are a 0/1 activation flag
//! (a capped cell returns the ceiling exactly, so the test is exact); components 3-5 are
//! log10(eta_raw / eta_applied) -- how many DECADES of diffusivity the ceiling removed
//! (0 where inactive). Both are per-cell, so the phdf gives the spatial picture (where the
//! caps bite relative to the core/accretion shock) for free.
namespace CapDiagIdx {
constexpr int flagO = 0;
constexpr int flagH = 1;
constexpr int flagA = 2;
constexpr int decO = 3;
constexpr int decH = 4;
constexpr int decA = 5;
constexpr int NCOMP = 6;
} // namespace CapDiagIdx

//! Fill "diff.capdiag" from the same per-cell state PrecomputeNonidealEta uses. Called
//! from PrecomputeNonidealEta itself (gated on the "nonideal_cap_diag" package flag) so
//! there is no second Wardle tensor solve and no extra task in the driver.
//! Volume-weighted history reductions over "diff.capdiag" (see CapDiagIdx). All are
//! plain sums so they compose across meshblocks/ranks; form fractions in analysis by
//! dividing by "cap-Vtot" (the summed volume of the same cells).
Real NonidealCapHstVol(MeshData<Real> *md);                   // cap-Vtot  = sum dV
Real NonidealCapHstFlag(MeshData<Real> *md, int comp);        // sum flag*dV
Real NonidealCapHstDecade(MeshData<Real> *md, int comp);      // sum flag*decades*dV
Real NonidealCapHstMassFlag(MeshData<Real> *md, int comp);    // sum flag*rho*dV

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
