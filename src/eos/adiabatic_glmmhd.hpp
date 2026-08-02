#ifndef EOS_ADIABATIC_GLMMHD_HPP_
#define EOS_ADIABATIC_GLMMHD_HPP_
//! \file eos.hpp
//  \brief defines class EquationOfState
//  Contains data and functions that implement the equation of state

// C headers

// C++ headers
#include <limits> // std::numeric_limits<float>

// Parthenon headers
#include "mesh/mesh.hpp"

// Athena headers
#include "../main.hpp"
#include "eos.hpp"
#include "eos_table.hpp"

using parthenon::MeshBlock;
using parthenon::MeshBlockData;
using parthenon::MeshBlockVarPack;
using parthenon::Real;

class AdiabaticGLMMHDEOS : public EquationOfState {
 public:
  // Optional tabulated protostellar EOS (second-core physics) is selected by
  // use_h2diss=true; then pressure/sound-speed/internal-energy come from the interpolated
  // table (eos_table.hpp). When use_h2diss=false the class is bit-identical to the ideal
  // gamma-law EOS (all table branches compile out at runtime; defaults keep old constructor
  // call sites valid).
  AdiabaticGLMMHDEOS(Real pressure_floor, Real density_floor, Real internal_e_floor,
                     Real velocity_ceiling, Real internal_e_ceiling, Real gamma,
                     bool use_h2diss = false,
                     EOSTable::EosTable eos_tab = EOSTable::EosTable(),
                     Real boris_ca_max = 0.0, bool ct_glm_inert = false)
      : EquationOfState(pressure_floor, density_floor, internal_e_floor, velocity_ceiling,
                        internal_e_ceiling),
        gamma_{gamma}, use_h2diss_{use_h2diss}, eos_tab_{eos_tab},
        boris_inv_ca2_{boris_ca_max > 0.0 ? 1.0 / (boris_ca_max * boris_ca_max) : 0.0},
        ct_glm_inert_{ct_glm_inert} {}

  // --- Boris / semi-relativistic Alfven-speed limiter ---------------------------------
  // Optional cap c_b on the Alfven speed used in the FORCE (Maxwell stress) and SIGNAL
  // speeds, to relax the timestep / tame the vA runaway in evacuated magnetized cells
  // (e.g. the first-core magnetic wall). Disabled (bit-identical) when boris_inv_ca2_==0
  // (boris_ca_max<=0). Formulation: replace B^2 -> B_eff^2 = B^2 / (1 + B^2/(rho c_b^2))
  // everywhere B exerts force or sets a speed; then vA_eff^2 = B_eff^2/rho -> c_b^2 as
  // B->inf, and -> real vA when vA<<c_b. The INDUCTION equation keeps the real B (Boris
  // changes how B acts on the fluid, not how B evolves). boris_inv_ca2_ = 1/c_b^2.
  KOKKOS_FORCEINLINE_FUNCTION
  bool BorisOn() const { return boris_inv_ca2_ > 0.0; }
  // Boris factor lambda = 1/(1 + vA^2/c_b^2), vA^2 = b2/rho. lambda in (0,1]; ==1 when off.
  KOKKOS_FORCEINLINE_FUNCTION
  Real BorisLambda(const Real d, const Real b2) const {
    return boris_inv_ca2_ > 0.0 ? 1.0 / (1.0 + (b2 / d) * boris_inv_ca2_) : 1.0;
  }

  void ConservedToPrimitive(MeshData<Real> *md) const override;

  KOKKOS_INLINE_FUNCTION
  Real GetGamma() const { return gamma_; }

  KOKKOS_INLINE_FUNCTION
  bool UseH2Diss() const { return use_h2diss_; }

  // Internal-energy DENSITY from (rho, pressure). Ideal: P/(gamma-1); H2-diss: Saha
  // inversion. Used by the Riemann solvers to reconstruct the L/R total energy so the
  // energy flux stays consistent with the (general) EOS.
  KOKKOS_INLINE_FUNCTION
  Real InternalEnergyFromPres(const Real rho, const Real pres) const {
    if (use_h2diss_) return eos_tab_.EintFromRhoPres(rho, pres);
    return pres / (gamma_ - 1.0);
  }

  // Accessor so downstream packages (e.g. RT matter-coupling) can share the same table.
  KOKKOS_INLINE_FUNCTION
  const EOSTable::EosTable &GetEosTable() const { return eos_tab_; }

  //----------------------------------------------------------------------------------------
  // \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
  // \brief returns adiabatic sound speed given vector of primitive variables
  // TODO(pgrete): need to fix idx defs
  KOKKOS_INLINE_FUNCTION
  Real SoundSpeed(const Real prim[NHYDRO]) const {
    if (use_h2diss_)
      return std::sqrt(eos_tab_.AsqFromRhoPres(prim[IDN], prim[IPR]));
    return std::sqrt(gamma_ * prim[IPR] / prim[IDN]);
  }
  // fast magnetosonic speed function for adiabatic EOS
  KOKKOS_INLINE_FUNCTION
  Real FastMagnetosonicSpeed(const Real d, const Real p, const Real bx, const Real by,
                             const Real bz) const {
    // asq here is rho*c_s^2 (= gamma*p for an ideal gas); use the general Saha c_s^2 when
    // the H2-dissociation EOS is active.
    Real asq = use_h2diss_ ? d * eos_tab_.AsqFromRhoPres(d, p) : gamma_ * p;
    // Boris limiter: replace B^2 -> lambda*B^2 in the magnetic (speed) terms so the fast
    // speed used for BOTH the HLLD signal speeds and the CFL timestep is capped consistently
    // with the Maxwell-stress scaling in the flux. lambda==1 (bit-identical) when disabled.
    const Real lam = BorisLambda(d, bx * bx + by * by + bz * bz);
    const Real bx2 = lam * bx * bx;
    Real ct2 = lam * (by * by + bz * bz);
    Real qsq = bx2 + ct2 + asq;
    Real tmp = bx2 + ct2 - asq;
    return std::sqrt(0.5 * (qsq + std::sqrt(tmp * tmp + 4.0 * asq * ct2)) / d);
  }
  //

  //----------------------------------------------------------------------------------------
  // \!fn Real EquationOfState::ConsToPrim(View4D cons, View4D prim, const int& k, const
  // int& j, const int& i) \brief Fills an array of primitives given an array of
  // conserveds, potentially updating the conserved with floors
  template <typename View4D>
  KOKKOS_INLINE_FUNCTION void ConsToPrim(View4D cons, View4D prim, const int &nhydro,
                                         const int &nscalars, const int &k, const int &j,
                                         const int &i) const {
    auto gam = GetGamma();
    auto gm1 = gam - 1.0;
    auto density_floor_ = GetDensityFloor();
    auto pressure_floor_ = GetPressureFloor();
    auto e_floor_ = GetInternalEFloor();

    auto velocity_ceiling_ = GetVelocityCeiling();
    auto e_ceiling_ = GetInternalECeiling();

    Real &u_d = cons(IDN, k, j, i);
    Real &u_m1 = cons(IM1, k, j, i);
    Real &u_m2 = cons(IM2, k, j, i);
    Real &u_m3 = cons(IM3, k, j, i);
    Real &u_e = cons(IEN, k, j, i);
    Real &u_b1 = cons(IB1, k, j, i);
    Real &u_b2 = cons(IB2, k, j, i);
    Real &u_b3 = cons(IB3, k, j, i);
    Real &u_psi = cons(IPS, k, j, i);

    Real &w_d = prim(IDN, k, j, i);
    Real &w_vx = prim(IV1, k, j, i);
    Real &w_vy = prim(IV2, k, j, i);
    Real &w_vz = prim(IV3, k, j, i);
    Real &w_p = prim(IPR, k, j, i);
    Real &w_Bx = prim(IB1, k, j, i);
    Real &w_By = prim(IB2, k, j, i);
    Real &w_Bz = prim(IB3, k, j, i);
    Real &w_psi = prim(IPS, k, j, i);

    // Let's apply floors explicitly, i.e., by default floor will be disabled (<=0)
    // and the code will fail if a negative density is encountered.
    PARTHENON_REQUIRE(u_d > 0.0 || density_floor_ > 0.0,
                      "Got negative density. Consider enabling first-order flux "
                      "correction or setting a reasonble density floor.");
    // apply density floor, without changing momentum or energy
    u_d = (u_d > density_floor_) ? u_d : density_floor_;
    w_d = u_d;

    Real di = 1.0 / u_d;
    w_vx = u_m1 * di;
    w_vy = u_m2 * di;
    w_vz = u_m3 * di;

    w_Bx = u_b1;
    w_By = u_b2;
    w_Bz = u_b3;
    // Constrained transport: the GLM/Dedner scalar is INERT and must be held at exactly
    // zero. Under CT the face field is the primary magnetic variable and the cell-centered
    // B is overwritten by its projection every substage, which DISCARDS the psi->B half of
    // the Dedner cleaning loop. psi is still sourced by the Riemann solver
    // (flxi[IPS] = c_h^2 * bxi) but never relieved, so it grows without bound -- measured
    // |psi|_max = 4.8e2 on the CT flagship vs 1.9e0 on the matched GLM run (244x). That
    // runaway then feeds straight back into every Riemann solve through the decoupled
    // normal field bxi = 1/2(BL+BR) - (psiR-psiL)/(2 c_h) (glmmhd_hlld.hpp), corrupting the
    // momentum and energy fluxes. Zeroing the CONSERVED psi here (not just the primitive)
    // covers restarts and ghost cells: every cons->prim pass, i.e. before every flux
    // evaluation, leaves psi == 0, so bxi reduces to the plain average and flxi[IPS] is
    // re-zeroed each stage. GLM path (ct_glm_inert_ == false) is untouched/bit-identical.
    if (ct_glm_inert_) {
      u_psi = 0.0;
    }
    w_psi = u_psi;

    Real e_k = 0.5 * di * (SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
    Real e_B = 0.5 * (SQR(u_b1) + SQR(u_b2) + SQR(u_b3));
    const Real e_int = u_e - e_k - e_B;
    // Table lookup interpolates in log10(e_int/rho): for e_int <= 0 it would return NaN,
    // and NaN compares false against every floor below (silently poisoning the cell).
    // Hand non-positive e_int to the ideal-gamma expression instead so w_p goes negative
    // and the existing PARTHENON_REQUIRE / pressure-floor logic catches it.
    w_p = (use_h2diss_ && e_int > 0.0) ? eos_tab_.PresFromRhoEint(u_d, e_int)
                                       : gm1 * e_int;

    // apply velocity ceiling. By default ceiling is std::numeric_limits<Real>::infinity()
    const Real w_v2 = SQR(w_vx) + SQR(w_vy) + SQR(w_vz);
    if (w_v2 > SQR(velocity_ceiling_)) {
      const Real w_v = sqrt(w_v2);
      w_vx *= velocity_ceiling_ / w_v;
      w_vy *= velocity_ceiling_ / w_v;
      w_vz *= velocity_ceiling_ / w_v;

      u_m1 *= velocity_ceiling_ / w_v;
      u_m2 *= velocity_ceiling_ / w_v;
      u_m3 *= velocity_ceiling_ / w_v;

      Real e_k_new = 0.5 * u_d * SQR(velocity_ceiling_);
      u_e -= e_k - e_k_new;
      e_k = e_k_new;
    }

    // Let's apply floors explicitly, i.e., by default floor will be disabled (<=0)
    // and the code will fail if a negative pressure is encountered.
    PARTHENON_REQUIRE(w_p > 0.0 || pressure_floor_ > 0.0 || e_floor_ > 0.0,
                      "Got negative pressure. Consider enabling first-order flux "
                      "correction or setting a reasonble pressure or temperature floor.");

    // Pressure floor (if present) takes precedence over temperature floor
    if ((pressure_floor_ > 0.0) && (w_p < pressure_floor_)) {
      // apply pressure floor, correct total energy. InternalEnergyFromPres keeps the
      // e<->P mapping consistent with the active EOS (ideal: pressure_floor/gm1,
      // bit-identical to the old expression; table: Saha inversion).
      u_e = InternalEnergyFromPres(u_d, pressure_floor_) + e_k + e_B;
      w_p = pressure_floor_;
    }

    // temperature (internal energy) based pressure floor. Under the tabulated EOS the
    // pressure corresponding to the specific-internal-energy floor comes from the table
    // (the gm1 conversion would set w_p inconsistent with the u_e written below); the
    // ideal-gas branch is unchanged.
    const Real eff_pressure_floor = (use_h2diss_ && e_floor_ > 0.0)
                                        ? eos_tab_.PresFromRhoEint(u_d, u_d * e_floor_)
                                        : gm1 * u_d * e_floor_;
    if (w_p < eff_pressure_floor) {
      // apply temperature floor, correct total energy
      u_e = (u_d * e_floor_) + e_k + e_B;
      w_p = eff_pressure_floor;
    }

    // temperature (internal energy) based pressure ceiling (guard the table against the
    // default e_ceiling = infinity: log10(inf) would clamp to the table edge and
    // spuriously trigger the ceiling).
    const Real eff_pressure_ceiling =
        (use_h2diss_ && e_ceiling_ < std::numeric_limits<Real>::infinity())
            ? eos_tab_.PresFromRhoEint(u_d, u_d * e_ceiling_)
            : gm1 * u_d * e_ceiling_;
    if (w_p > eff_pressure_ceiling) {
      // apply temperature ceiling, correct total energy
      u_e = (u_d * e_ceiling_) + e_k + e_B;
      w_p = eff_pressure_ceiling;
    }

    // Convert passive scalars
    for (auto n = nhydro; n < nhydro + nscalars; ++n) {
      prim(n, k, j, i) = cons(n, k, j, i) * di;
    }
  }

 private:
  Real gamma_; // ratio of specific heats
  bool use_h2diss_;            // enable tabulated protostellar EOS (second-core physics)
  EOSTable::EosTable eos_tab_; // interpolated table (View handles; device-copyable)
  Real boris_inv_ca2_;        // 1/c_b^2 for the Boris Alfven-speed limiter (0 = disabled)
  bool ct_glm_inert_;         // CT: hold the GLM scalar psi identically zero (see ConsToPrim)
};

#endif // EOS_ADIABATIC_GLMMHD_HPP_
