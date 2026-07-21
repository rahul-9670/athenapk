//========================================================================================
// AthenaPK - Bonnor-Ebert sphere collapse problem generator
// Ported from Athena++'s collapse.cpp (Kengo Tomida).
// Adapted for AthenaPK's hydrodynamics interface and self-gravity package.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
//! \file collapse_be.cpp
//  \brief Bonnor-Ebert sphere collapse with zero-Dirichlet gravity BCs
//         (multipole BCs not yet supported in AthenaPK self-gravity).
//         Large-box setup: domain size ~4x sphere radius in each direction,
//         ambient density = BE profile value at r = rc (pressure-continuous).

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <parthenon/package.hpp>

#include "../main.hpp"
#include "../units.hpp" // canonical Heaviside-Lorentz magnetic-unit definition (audit #3)
#include "pgen.hpp"

namespace collapse_be {
using namespace parthenon::driver::prelude;

// ---------------------------------------------------------------------------
// Dimensionless constants from Tomida's normalization
// In code units: 4piG = 1, c_s(T=10K) = 1, central density of BE sphere = 1
// ---------------------------------------------------------------------------
constexpr Real four_pi_G_code = 1.0;
constexpr Real rc_code        = 6.45;        // BE radius in code units
constexpr Real rcsq_code      = 26.0 / 3.0;  // BE profile scale^2
constexpr Real bemass_code    = 197.561;     // total mass of critical BE sphere

// Dimensional constants (cgs)
constexpr Real cs10 = 1.9e4;        // cm/s, sound speed at 10 K
constexpr Real msun = 1.9891e33;    // g
constexpr Real au   = 1.4959787e13; // cm
constexpr Real yr   = 3.15569e7;    // s
constexpr Real G    = 6.67259e-8;   // cgs

// These hold the derived physical parameters for logging and the cooling fn.
// Parthenon doesn't give us a natural "namespace-scoped globals" pattern like
// Athena++, so we stash them in the Hydro package params after problem init.
struct CollapseParams {
  Real rho0;     // density normalization (cgs)
  Real rhocrit;  // barotropic transition density in CODE units
  Real f;        // BE density enhancement factor
  Real rc;       // sphere radius in code units
  Real omega;    // angular velocity in code units
  Real amp;      // m=2 perturbation amplitude
  Real gamma_eos;
};

// Approximated BE profile (Tomida 2011 PhD thesis)
KOKKOS_INLINE_FUNCTION
Real BEProfile(Real r) {
  return std::pow(1.0 + r * r / rcsq_code, -1.5);
}

// ---------------------------------------------------------------------------
// Problem generator: initialize cons variables on each block.
// No rotation-vs-ambient distinction — inside sphere has solid-body rotation,
// outside has zero velocity and ambient density = rho(rc).
// ---------------------------------------------------------------------------
void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
  // For the tabulated protostellar EOS, cold molecular H2 (10 K, rotation frozen) is
  // monatomic-like: e_int = 1.5 p (gamma=5/3), not p/(gamma-1). Seeding the ideal value
  // would make the table read the IC as hot. (At 10 K the table's pressure already equals
  // the isothermal p=rho since it uses the same mu=2.33; only the energy differs.)
  const bool use_h2diss = (pin->GetOrAddString("hydro", "eos", "adiabatic") == "hydrogen");
  const Real gam = pin->GetReal("hydro", "gamma");
  const Real gm1 = gam - 1.0;
  const Real igm1 = 1.0 / gm1;
  const Real e_int_over_p = use_h2diss ? 1.5 : igm1; // e_int = e_int_over_p * p

  // Read problem params
  const Real mass_msun = pin->GetReal("problem/collapse_be", "mass");
  const Real temp_K    = pin->GetReal("problem/collapse_be", "temperature");
  const Real f         = pin->GetReal("problem/collapse_be", "f");
  const Real amp       = pin->GetOrAddReal("problem/collapse_be", "amp", 0.0);
  const Real omegatff  = pin->GetOrAddReal("problem/collapse_be", "omegatff", 0.0);
  const Real rhocrit_cgs = pin->GetReal("problem/collapse_be", "rhocrit");

  // Derive normalization units
  const Real m0   = mass_msun * msun / (bemass_code * f);   // total mass = bemass*f in code units
  const Real v0   = cs10 * std::sqrt(temp_K / 10.0);         // code c_s = 1
  const Real rho0 = std::pow(v0, 6) / (m0 * m0) / (64.0 * M_PI * M_PI * M_PI * G * G * G);
  const Real t0   = 1.0 / std::sqrt(4.0 * M_PI * G * rho0);  // code t = 1
  const Real l0   = v0 * t0;

  // Angular velocity: omega * t_ff = omegatff
  const Real tff_code = std::sqrt(3.0 / (8.0 * f)) * M_PI;
  const Real omega_code = omegatff / tff_code;

  // rhocrit in code units
  const Real rhocrit_code = rhocrit_cgs / rho0;

  // Stash params in the hydro package for the cooling function to read.
  // (Only rank 0 prints; all ranks need the params.)
  auto hydro_pkg = pmb->packages.Get("Hydro");
  const bool mhd = (hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd);
  const Real B0z = mhd ? pin->GetOrAddReal("problem/collapse_be", "B0z", 0.0) : 0.0;
  if (!hydro_pkg->AllParams().hasKey("collapse_be_rhocrit")) {
    hydro_pkg->AddParam("collapse_be_rhocrit", rhocrit_code);
    hydro_pkg->AddParam("collapse_be_mhd", mhd);
    hydro_pkg->AddParam("collapse_be_rc", rc_code);
    hydro_pkg->AddParam("collapse_be_gamma", gam);

    // === Unit conversion factors to cgs ===
    // Multiply a code-unit value by these to obtain cgs.
    const Real code_length_cgs   = l0;
    const Real code_time_cgs     = t0;
    const Real code_mass_cgs     = m0;
    const Real code_density_cgs  = rho0;
    const Real code_velocity_cgs = l0 / t0;
    const Real code_energy_cgs   = m0 * (l0/t0) * (l0/t0);
    // Audit fix #3 / flagship Phase 1: the Heaviside-Lorentz magnetic unit
    // B_cgs = sqrt(4*pi*rho0)*v0 (P_mag = B^2/2) is now taken from the ONE definition in
    // units.hpp (Units::CodeMagneticCgs), the same formula PhysicalUnits uses -- so this
    // sidecar cannot drift from the code's internal B convention. The BE derivation gives
    // mass_unit = rho0*l0^3 exactly, so CodeMagneticCgs(m0,l0,t0) == sqrt(4*pi*rho0)*v0.
    // (Previously the sqrt(4*pi) was dropped, under-reporting physical B by ~3.545.)
    const Real code_bfield_cgs   = Units::CodeMagneticCgs(code_mass_cgs, code_length_cgs,
                                                          code_time_cgs);

    hydro_pkg->AddParam("units/code_length_cgs",   code_length_cgs);
    hydro_pkg->AddParam("units/code_time_cgs",     code_time_cgs);
    hydro_pkg->AddParam("units/code_mass_cgs",     code_mass_cgs);
    hydro_pkg->AddParam("units/code_density_cgs",  code_density_cgs);
    hydro_pkg->AddParam("units/code_velocity_cgs", code_velocity_cgs);
    hydro_pkg->AddParam("units/code_energy_cgs",   code_energy_cgs);
    hydro_pkg->AddParam("units/code_bfield_cgs",   code_bfield_cgs);

    // Sidecar JSON dump (rank 0 only; once per run).
    if (parthenon::Globals::my_rank == 0) {
      std::ofstream js("units.json");
      js.precision(12);
      js << "{\n"
         << "  \"code_length_cgs\":   " << code_length_cgs   << ",\n"
         << "  \"code_time_cgs\":     " << code_time_cgs     << ",\n"
         << "  \"code_mass_cgs\":     " << code_mass_cgs     << ",\n"
         << "  \"code_density_cgs\":  " << code_density_cgs  << ",\n"
         << "  \"code_velocity_cgs\": " << code_velocity_cgs << ",\n"
         << "  \"code_energy_cgs\":   " << code_energy_cgs   << ",\n"
         << "  \"code_bfield_cgs\":   " << code_bfield_cgs   << "\n"
         << "}\n";
    }
  }

  // Log (rank 0 only)
  if (parthenon::Globals::my_rank == 0 && pmb->gid == 0) {
    std::cout << "\n---  Dimensional parameters of the simulation  ---\n"
              << "Total mass          : " << mass_msun << " [Msun]\n"
              << "Initial temperature : " << temp_K    << " [K]\n"
              << "Sound speed         : " << v0        << " [cm/s]\n"
              << "Central density     : " << rho0 * f  << " [g/cm^3]\n"
              << "Cloud radius        : " << rc_code * l0 / au << " [au]\n"
              << "Free-fall time      : " << tff_code * t0 / yr << " [yr]\n"
              << "rhocrit (cgs)       : " << rhocrit_cgs << " [g/cm^3]\n"
              << "rhocrit (code)      : " << rhocrit_code << "\n"
              << "Density enhancement : " << f << "\n\n"
              << "---   Normalization Units    ---\n"
              << "Mass                : " << m0 / msun << " [Msun]\n"
              << "Length              : " << l0 / au   << " [au]\n"
              << "Time                : " << t0 / yr   << " [yr]\n"
              << "Density             : " << rho0      << " [g/cm^3]\n\n"
              << "--- Dimensionless parameters ---\n"
              << "Total mass          : " << bemass_code * f << "\n"
              << "Sound speed         : 1.0\n"
              << "Central density     : 1.0\n"
              << "Cloud radius        : " << rc_code << "\n"
              << "Free-fall time      : " << tff_code << "\n"
              << "m=2 perturbation    : " << amp << "\n"
              << "Omega * tff         : " << omegatff << std::endl;
  }

  auto &data = pmb->meshblock_data.Get();
  auto &u_dev = data->Get("cons").data;
  auto u = u_dev.GetHostMirrorAndCopy();

  parthenon::IndexRange ib = pmb->cellbounds.GetBoundsI(parthenon::IndexDomain::entire);
  parthenon::IndexRange jb = pmb->cellbounds.GetBoundsJ(parthenon::IndexDomain::entire);
  parthenon::IndexRange kb = pmb->cellbounds.GetBoundsK(parthenon::IndexDomain::entire);

  auto &coords = pmb->coords;
  // Ambient density: BE profile value at r = rc (pressure-continuous)
  const Real rho_amb = f * BEProfile(rc_code);

  // Chemistry species ICs (gated on passive scalars). Species ride on the passive
  // scalars as scalar_density_i = rho * x_i; seed all-atomic hydrogen (x_H=1, x_H2=0)
  // so a chemistry run starts from a physical, H-nucleus-conserving state. The first
  // scalar_density component sits just past the hydro vars in the cons pack.
  const int nscalars_ic = pin->GetOrAddInteger("hydro", "nscalars", 0);
  const int iscal0 = mhd ? (IPS + 1) : NHYDRO;

  for (int k = kb.s; k <= kb.e; ++k) {
    Real z = coords.Xc<3>(k);
    for (int j = jb.s; j <= jb.e; ++j) {
      Real y = coords.Xc<2>(j);
      for (int i = ib.s; i <= ib.e; ++i) {
        Real x = coords.Xc<1>(i);
        Real r = std::sqrt(x * x + y * y + z * z);

        Real rho, v1, v2, v3;
        if (r < rc_code) {
          // Inside sphere: BE profile + m=2 azimuthal perturbation + solid-body rotation
          const Real phi = std::atan2(y, x);
          rho = f * BEProfile(r) * (1.0 + amp * (r*r/(rc_code*rc_code)) * std::cos(2.0 * phi));
          v1 =  omega_code * y;
          v2 = -omega_code * x;
          v3 = 0.0;
        } else {
          // Outside: ambient, at rest
          rho = rho_amb;
          v1 = 0.0;
          v2 = 0.0;
          v3 = 0.0;
        }

        // Pressure: isothermal c_s = 1 (T=10K), so the BE-equilibrium pressure is
        // p = rho * c_s^2 = rho. This matches the barotropic floor the cooling source
        // enforces every step (e_th = rho/(gamma-1) for rho << rhocrit  =>  p = rho),
        // so there is no factor-gamma transient on the first step. NOTE: the initial
        // pressure is cosmetic — the barotropic cooling overwrites e_th on step 1, so
        // p=rho vs p=rho/gam give bit-identical evolution (verified A/B 2026-06-22).
        const Real p = rho;

        u(IDN, k, j, i) = rho;
        u(IM1, k, j, i) = rho * v1;
        u(IM2, k, j, i) = rho * v2;
        u(IM3, k, j, i) = rho * v3;
        u(IEN, k, j, i) = p * e_int_over_p + 0.5 * rho * (v1*v1 + v2*v2 + v3*v3);
        if (mhd) {
          u(IB1, k, j, i) = 0.0;
          u(IB2, k, j, i) = 0.0;
          u(IB3, k, j, i) = B0z;
          u(IPS, k, j, i) = 0.0;
          u(IEN, k, j, i) += 0.5 * B0z * B0z;
        }
        // Species ICs: scalar_density_i = rho * x_i. Layout depends on the chemistry
        // network (selected by the number of scalars):
        //   >=5 -> gow17_reduced {H2,H+,C+,CO,e-}: molecular-cloud core (H in H2, C in C+)
        //   ==2 -> H2 {H,H2}: all-atomic hydrogen
        if (nscalars_ic >= 5) {
          u(iscal0 + 0, k, j, i) = rho * 0.5;     // x_H2 = 0.5 (all H molecular)
          u(iscal0 + 1, k, j, i) = rho * 1.0e-8;  // x_H+
          u(iscal0 + 2, k, j, i) = rho * 1.6e-4;  // x_C+  = total gas-phase C
          u(iscal0 + 3, k, j, i) = rho * 1.0e-20; // x_CO
          u(iscal0 + 4, k, j, i) = rho * 1.6e-4;  // x_e   ~ x_C+ (+ x_H+)
          for (int s = 5; s < nscalars_ic; ++s) u(iscal0 + s, k, j, i) = 0.0;
        } else if (nscalars_ic >= 2) {
          u(iscal0 + 0, k, j, i) = rho; // x_H  = 1
          u(iscal0 + 1, k, j, i) = 0.0; // x_H2 = 0
          for (int s = 2; s < nscalars_ic; ++s) u(iscal0 + s, k, j, i) = 0.0;
        }
      }
    }
  }
  // ============================================================================
  // Initial-only turbulent velocity field (gated by turb_mach > 0).
  // Adds a random divergence-free + compressive velocity perturbation on top
  // of the deterministic IC above. The field is scaled to `turb_mach` via the
  // ANALYTIC ensemble RMS of the mode spectrum, and the ANALYTIC region-mean
  // velocity is subtracted so the cloud doesn't translate as a unit (both are
  // closed-form so every rank/block applies identical factors; see notes below).
  // ============================================================================
  const Real turb_mach = pin->GetOrAddReal("problem/collapse_be", "turb_mach", 0.0);
  if (turb_mach > 0.0) {
    const int   k_min   = pin->GetOrAddInteger("problem/collapse_be", "turb_kmin", 2);
    const int   k_max   = pin->GetOrAddInteger("problem/collapse_be", "turb_kmax", 6);
    const int   nmodes  = pin->GetOrAddInteger("problem/collapse_be", "turb_nmodes", 64);
    const Real  zeta    = pin->GetOrAddReal("problem/collapse_be", "turb_zeta", 0.5);
    const Real  alpha_s = pin->GetOrAddReal("problem/collapse_be", "turb_alpha", 2.0);
    const std::string region = pin->GetOrAddString(
        "problem/collapse_be", "turb_region", std::string("sphere"));
    const uint32_t seed = static_cast<uint32_t>(
        pin->GetOrAddInteger("problem/collapse_be", "turb_seed", 12345));

    const Real x1min_ = pin->GetReal("parthenon/mesh", "x1min");
    const Real x1max_ = pin->GetReal("parthenon/mesh", "x1max");
    const Real Lbox  = x1max_ - x1min_;
    const Real two_pi_over_L = 2.0 * M_PI / Lbox;

    // Build the SAME mode list on every rank (deterministic with seed).
    std::mt19937 mrng(seed);
    std::uniform_real_distribution<Real> uphase(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<Real> ucos(-1.0, 1.0);
    std::uniform_real_distribution<Real> uang(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<Real> ukmag(static_cast<Real>(k_min),
                                               static_cast<Real>(k_max));

    std::vector<std::array<Real,3>> k_vec(nmodes), phase(nmodes),
                                     eps_sol(nmodes), eps_comp(nmodes);
    std::vector<Real> t_amp(nmodes);

    for (int m = 0; m < nmodes; ++m) {
      const Real kmag = ukmag(mrng);
      const Real cos_t = ucos(mrng);
      const Real sin_t = std::sqrt(std::max(0.0, 1.0 - cos_t * cos_t));
      const Real phi_k = uang(mrng);
      const Real kxh = sin_t * std::cos(phi_k);
      const Real kyh = sin_t * std::sin(phi_k);
      const Real kzh = cos_t;

      k_vec[m] = {kmag * two_pi_over_L * kxh,
                  kmag * two_pi_over_L * kyh,
                  kmag * two_pi_over_L * kzh};
      t_amp[m] = std::pow(kmag, -0.5 * alpha_s);
      phase[m] = {uphase(mrng), uphase(mrng), uphase(mrng)};

      // Solenoidal direction: pick non-aligned axis, cross with k_hat.
      Real ax = 1.0, ay = 0.0, az = 0.0;
      if (std::fabs(kxh) > std::fabs(kyh) && std::fabs(kxh) > std::fabs(kzh)) {
        ax = 0.0; ay = 1.0; az = 0.0;
      }
      Real sx = ay * kzh - az * kyh;
      Real sy = az * kxh - ax * kzh;
      Real sz = ax * kyh - ay * kxh;
      Real snorm = std::sqrt(sx*sx + sy*sy + sz*sz);
      if (snorm < 1e-12) snorm = 1.0;
      eps_sol[m]  = {sx/snorm, sy/snorm, sz/snorm};
      eps_comp[m] = {kxh, kyh, kzh};
    }

    // Evaluate field on this block's local cells (host scratch).
    const int nx_loc = ib.e - ib.s + 1;
    const int ny_loc = jb.e - jb.s + 1;
    const int nz_loc = kb.e - kb.s + 1;
    std::vector<Real> v_turb(3 * nx_loc * ny_loc * nz_loc, 0.0);
    auto vt = [&](int c, int kk, int jj, int ii) -> Real & {
      return v_turb[((c * nz_loc + (kk - kb.s)) * ny_loc + (jj - jb.s)) * nx_loc
                    + (ii - ib.s)];
    };

    for (int k = kb.s; k <= kb.e; ++k) {
      const Real z = coords.Xc<3>(k);
      for (int j = jb.s; j <= jb.e; ++j) {
        const Real y = coords.Xc<2>(j);
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real x = coords.Xc<1>(i);
          Real vx_t = 0.0, vy_t = 0.0, vz_t = 0.0;
          for (int m = 0; m < nmodes; ++m) {
            const Real kdotr = k_vec[m][0] * x + k_vec[m][1] * y + k_vec[m][2] * z;
            const Real w_sol  = std::sqrt(zeta) * t_amp[m]
                              * std::cos(kdotr + phase[m][0]);
            const Real w_comp = std::sqrt(1.0 - zeta) * t_amp[m]
                              * std::cos(kdotr + phase[m][1]);
            vx_t += w_sol * eps_sol[m][0]  + w_comp * eps_comp[m][0];
            vy_t += w_sol * eps_sol[m][1]  + w_comp * eps_comp[m][1];
            vz_t += w_sol * eps_sol[m][2]  + w_comp * eps_comp[m][2];
          }
          vt(0, k, j, i) = vx_t;
          vt(1, k, j, i) = vy_t;
          vt(2, k, j, i) = vz_t;
        }
      }
    }

    // -----------------------------------------------------------------------
    // Apply analytical normalization to the unscaled mode-sum field.
    // We do NOT use a mass-weighted RMS reduction across the mesh — that
    // approach fails when ProblemGenerator is called per-block and block 0
    // (which prints diagnostics) sits entirely outside the sphere. With
    // 1 MPI rank the MPI_Allreduce returns this block's own data, not the
    // global mesh state, so per-block rescaling becomes inconsistent and
    // outside-sphere blocks silently skip turbulence entirely.
    //
    // Instead we compute the field's expected RMS analytically from the
    // mode spectrum and apply a fixed scale factor identical on every block.
    //
    // For an isotropic random field with t_amp[m] = kmag^(-alpha/2),
    // unit-variance cos terms, and orthogonal sol/comp projections each
    // unit-norm, the per-component variance is:
    //   <v_x^2> = sum_m [zeta * <cos^2> * eps_sol_x^2 + (1-zeta) * <cos^2> * eps_comp_x^2]
    //          * t_amp[m]^2
    // Summed over (x,y,z), the eps^2 components sum to 1 for each mode (since
    // both eps_sol and eps_comp are unit vectors), so:
    //   <|v|^2> = sum_m t_amp[m]^2 * <cos^2>(zeta + (1-zeta)) = 0.5 * sum_m t_amp[m]^2
    // -----------------------------------------------------------------------
    Real sum_amp_sq = 0.0;
    for (int m = 0; m < nmodes; ++m) {
      sum_amp_sq += t_amp[m] * t_amp[m];
    }
    const Real vrms_analytic = std::sqrt(0.5 * sum_amp_sq);
    const Real scale = (vrms_analytic > 0.0) ? (turb_mach / vrms_analytic) : 0.0;

    // -----------------------------------------------------------------------
    // Analytic mean of the mode-sum field over the application region, subtracted
    // before scaling so the cloud receives no net momentum kick (for ~64 modes the
    // sphere-restricted mean is O(a few %) of v_rms). Computed analytically -- so it
    // is deterministic and identical on every rank -- for the same reason as the RMS
    // normalization above: a mesh-wide reduction is not per-block consistent here.
    //   sphere (centered on origin): <cos(k.r + phi)> = cos(phi) * W(|k| rc),
    //       W(x) = 3 (sin x - x cos x) / x^3   (-> 1 as x -> 0)
    //   box:  <cos(k.r + phi)> = cos(k.r_c + phi) * prod_d sinc(k_d L/2)
    // -----------------------------------------------------------------------
    const bool in_sphere = (region == "sphere");
    Real vmean[3] = {0.0, 0.0, 0.0};
    for (int m = 0; m < nmodes; ++m) {
      Real wcos_sol, wcos_comp; // region-mean of the sol/comp cosine terms
      if (in_sphere) {
        const Real kmagp =
            std::sqrt(k_vec[m][0] * k_vec[m][0] + k_vec[m][1] * k_vec[m][1] +
                      k_vec[m][2] * k_vec[m][2]);
        const Real x = kmagp * rc_code;
        const Real W =
            (x < 1.0e-4) ? 1.0 : 3.0 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
        wcos_sol = W * std::cos(phase[m][0]);
        wcos_comp = W * std::cos(phase[m][1]);
      } else {
        // Cubic box; same center in each dimension (matches the k-space setup above).
        const Real xc = 0.5 * (x1min_ + x1max_);
        Real sincs = 1.0;
        Real kdotc = 0.0;
        for (int d = 0; d < 3; ++d) {
          const Real arg = 0.5 * k_vec[m][d] * Lbox;
          sincs *= (std::fabs(arg) < 1.0e-8) ? 1.0 : std::sin(arg) / arg;
          kdotc += k_vec[m][d] * xc;
        }
        wcos_sol = sincs * std::cos(kdotc + phase[m][0]);
        wcos_comp = sincs * std::cos(kdotc + phase[m][1]);
      }
      const Real a_sol = std::sqrt(zeta) * t_amp[m] * wcos_sol;
      const Real a_comp = std::sqrt(1.0 - zeta) * t_amp[m] * wcos_comp;
      for (int c = 0; c < 3; ++c) {
        vmean[c] += a_sol * eps_sol[m][c] + a_comp * eps_comp[m][c];
      }
    }

    if (parthenon::Globals::my_rank == 0 && pmb->gid == 0) {
      std::cout << "\n--- Initial turbulence ---\n"
                << "Region              : " << region << "\n"
                << "Target Mach (RMS)   : " << turb_mach << "\n"
                << "Modes               : " << nmodes
                << "  k in [" << k_min << ", " << k_max << "]"
                << "  zeta=" << zeta << "  alpha=" << alpha_s << "\n"
                << "Analytic vrms       : " << vrms_analytic << " (code v)\n"
                << "Scale factor        : " << scale << "\n"
                << "Mean subtracted     : (" << vmean[0] * scale << ", "
                << vmean[1] * scale << ", " << vmean[2] * scale << ") (code v)"
                << std::endl;
    }

    // Apply the scaled, mean-subtracted turbulent velocity to the region's cells.
    for (int k = kb.s; k <= kb.e; ++k) {
      const Real z = coords.Xc<3>(k);
      for (int j = jb.s; j <= jb.e; ++j) {
        const Real y = coords.Xc<2>(j);
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real x = coords.Xc<1>(i);
          if (in_sphere && (x*x + y*y + z*z) >= rc_code * rc_code) continue;
          const Real dvx = (vt(0, k, j, i) - vmean[0]) * scale;
          const Real dvy = (vt(1, k, j, i) - vmean[1]) * scale;
          const Real dvz = (vt(2, k, j, i) - vmean[2]) * scale;
          const Real rho     = u(IDN, k, j, i);
          const Real vx_old  = u(IM1, k, j, i) / rho;
          const Real vy_old  = u(IM2, k, j, i) / rho;
          const Real vz_old  = u(IM3, k, j, i) / rho;
          const Real vx_new  = vx_old + dvx;
          const Real vy_new  = vy_old + dvy;
          const Real vz_new  = vz_old + dvz;
          u(IM1, k, j, i) = rho * vx_new;
          u(IM2, k, j, i) = rho * vy_new;
          u(IM3, k, j, i) = rho * vz_new;
          const Real ke_old = 0.5 * rho * (vx_old*vx_old + vy_old*vy_old + vz_old*vz_old);
          const Real ke_new = 0.5 * rho * (vx_new*vx_new + vy_new*vy_new + vz_new*vz_new);
          u(IEN, k, j, i) += (ke_new - ke_old);
        }
      }
    }
  }
  u_dev.DeepCopy(u);

  // ============================================================================
  // Seed the M1 radiation field in LTE with the initial gas (only when the
  // radiation package is active). Er = arad * T^4, isotropic (Fr = 0). Without
  // this, Er defaults to 0 and the matter-coupling solve would drain a ~arad*T^4
  // chunk of gas thermal energy on the first step (a spurious t=0 cooling
  // transient). Mirrors the ir = T^4 seed used in the Athena++ collapse_rad pgen.
  // ============================================================================
  if (pmb->packages.AllPackages().count("radiation") > 0) {
    auto rad_pkg = pmb->packages.Get("radiation");
    const Real arad = rad_pkg->Param<Real>("arad");
    auto Er = data->Get("rad.Er").data.GetHostMirrorAndCopy();
    auto Fx = data->Get("rad.Fr1").data.GetHostMirrorAndCopy();
    auto Fy = data->Get("rad.Fr2").data.GetHostMirrorAndCopy();
    auto Fz = data->Get("rad.Fr3").data.GetHostMirrorAndCopy();
    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real rho = u(IDN, k, j, i);
          const Real ke  = 0.5 / rho * (u(IM1, k, j, i) * u(IM1, k, j, i)
                                      + u(IM2, k, j, i) * u(IM2, k, j, i)
                                      + u(IM3, k, j, i) * u(IM3, k, j, i));
          const Real me  = mhd ? 0.5 * (u(IB1, k, j, i) * u(IB1, k, j, i)
                                      + u(IB2, k, j, i) * u(IB2, k, j, i)
                                      + u(IB3, k, j, i) * u(IB3, k, j, i))
                               : 0.0;
          const Real p   = gm1 * (u(IEN, k, j, i) - ke - me);
          // T_code (= 1 at the 10 K IC). For the tabulated EOS the ideal relation p/rho
          // does not hold; the BE-sphere IC is isothermal at 10 K everywhere => T_code=1.
          const Real T   = use_h2diss ? 1.0 : p / rho;
          Er(0, k, j, i) = arad * T * T * T * T;
          Fx(0, k, j, i) = 0.0;
          Fy(0, k, j, i) = 0.0;
          Fz(0, k, j, i) = 0.0;
        }
      }
    }
    data->Get("rad.Er").data.DeepCopy(Er);
    data->Get("rad.Fr1").data.DeepCopy(Fx);
    data->Get("rad.Fr2").data.DeepCopy(Fy);
    data->Get("rad.Fr3").data.DeepCopy(Fz);
  }
}

// ---------------------------------------------------------------------------
// Barotropic cooling source term (Masunaga-Inutsuka / Tomida form):
//   e_th = rho/(gamma-1) * sqrt(1 + (rho/rhocrit)^(2(gamma-1)))
// For rho << rhocrit: e_th ~ rho/(gamma-1) (isothermal, c_s = 1)
// For rho >> rhocrit: e_th ~ (rho/(gamma-1)) * (rho/rhocrit)^(gamma-1) (adiabatic)
// This OVERWRITES the thermal energy every step — it's a barotropic EOS
// enforcement, not a differential cooling equation.
//
// Also zeroes momentum outside the sphere (fixed boundary condition on velocity).
// ---------------------------------------------------------------------------
TaskStatus ApplyBarotropicCooling(MeshData<Real> *md, const parthenon::SimTime &tm,
                                   const Real dt) {
  auto pm = md->GetParentPointer();
  auto hydro_pkg = pm->packages.Get("Hydro");
  const Real rhocrit = hydro_pkg->Param<Real>("collapse_be_rhocrit");
  const Real rc = hydro_pkg->Param<Real>("collapse_be_rc");
  const Real gam = hydro_pkg->Param<Real>("collapse_be_gamma");
  const bool mhd = hydro_pkg->Param<bool>("collapse_be_mhd");
  const Real gm1 = gam - 1.0;
  const Real igm1 = 1.0 / gm1;
  const Real rcsq = rc * rc;

  // When the M1 radiation package is active it OWNS the gas thermal energy: skip the
  // barotropic e_th overwrite (Radiation::MatterCoupling sets T from absorption/emission
  // instead). The outside-sphere momentum BC below still applies in both cases.
  const bool rad_owns_energy = pm->packages.AllPackages().count("radiation") > 0;

  const auto &cons_pack = md->PackVariables(std::vector<std::string>{"cons"});

  parthenon::IndexRange ib = md->GetBoundsI(parthenon::IndexDomain::interior);
  parthenon::IndexRange jb = md->GetBoundsJ(parthenon::IndexDomain::interior);
  parthenon::IndexRange kb = md->GetBoundsK(parthenon::IndexDomain::interior);
  const int nblocks = md->NumBlocks();

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "BEcollapse::BaroCool", parthenon::DevExecSpace(),
      0, nblocks - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        auto &cons = cons_pack(b);
        const auto &coords = cons_pack.GetCoords(b);
        const Real x = coords.template Xc<1>(i);
        const Real y = coords.template Xc<2>(j);
        const Real z = coords.template Xc<3>(k);
        const Real r2 = x*x + y*y + z*z;

        // Outside sphere: zero velocity (fixed BC on momentum)
        if (r2 > rcsq) {
          // When radiation owns the gas energy the barotropic overwrite below is skipped,
          // so the kinetic energy of the killed momentum must be removed from IEN too --
          // otherwise it is silently converted to heat (a spurious secular heating of the
          // ambient gas every step). In the non-RT branch IEN is fully overwritten below,
          // so this subtraction is a no-op there by construction.
          if (rad_owns_energy) {
            const Real rho_bc = cons(IDN, k, j, i);
            cons(IEN, k, j, i) -=
                0.5 / rho_bc *
                (cons(IM1, k, j, i) * cons(IM1, k, j, i) +
                 cons(IM2, k, j, i) * cons(IM2, k, j, i) +
                 cons(IM3, k, j, i) * cons(IM3, k, j, i));
          }
          cons(IM1, k, j, i) = 0.0;
          cons(IM2, k, j, i) = 0.0;
          cons(IM3, k, j, i) = 0.0;
        }

        // Barotropic EOS enforcement (everywhere) -- SKIPPED when radiation owns energy.
        if (!rad_owns_energy) {
          const Real rho = cons(IDN, k, j, i);
          const Real ke = 0.5 / rho * (cons(IM1, k, j, i) * cons(IM1, k, j, i)
                                     + cons(IM2, k, j, i) * cons(IM2, k, j, i)
                                     + cons(IM3, k, j, i) * cons(IM3, k, j, i));
          const Real me = mhd ? 0.5 * (cons(IB1, k, j, i) * cons(IB1, k, j, i)
                                     + cons(IB2, k, j, i) * cons(IB2, k, j, i)
                                     + cons(IB3, k, j, i) * cons(IB3, k, j, i))
                               : 0.0;
          const Real te = igm1 * rho * std::sqrt(1.0 + std::pow(rho / rhocrit, 2.0 * gm1));
          cons(IEN, k, j, i) = te + ke + me;
        }
      });
  return TaskStatus::complete;
}

} // namespace collapse_be
