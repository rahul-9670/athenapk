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
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include <parthenon/package.hpp>

#include "../main.hpp"
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
  const Real gam = pin->GetReal("hydro", "gamma");
  const Real gm1 = gam - 1.0;
  const Real igm1 = 1.0 / gm1;

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
    const Real code_bfield_cgs   = std::sqrt(rho0) * (l0 / t0);

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

        // Pressure: c_s = 1 everywhere initially (isothermal T=10K-equivalent)
        // p = rho * c_s^2 / gamma (adiabatic relation — cooling fn will enforce barotropic)
        const Real p = rho / gam;

        u(IDN, k, j, i) = rho;
        u(IM1, k, j, i) = rho * v1;
        u(IM2, k, j, i) = rho * v2;
        u(IM3, k, j, i) = rho * v3;
        u(IEN, k, j, i) = p * igm1 + 0.5 * rho * (v1*v1 + v2*v2 + v3*v3);
        if (mhd) {
          u(IB1, k, j, i) = 0.0;
          u(IB2, k, j, i) = 0.0;
          u(IB3, k, j, i) = B0z;
          u(IPS, k, j, i) = 0.0;
          u(IEN, k, j, i) += 0.5 * B0z * B0z;
        }
      }
    }
  }
  u_dev.DeepCopy(u);
}

// ---------------------------------------------------------------------------
// Barotropic cooling source term (Masunaga-Inutsuka / Tomida form):
//   e_th = rho/(gamma-1) * sqrt(1 + (rho/rhocrit)^(2(gamma-1)))
// For rho << rhocrit: e_th ~ rho/(gamma-1) (isothermal, c_s = 1)
// For rho >> rhocrit: e_th ~ (rho/(gamma-1)) * (rho/rhocrit)^(gamma-1) (adiabatic)
// This OVERWRITES the thermal energy every step — it's a barotropic EOS
// enforcement, not a differential cooling equation. Same as Athena++ setup.
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
          cons(IM1, k, j, i) = 0.0;
          cons(IM2, k, j, i) = 0.0;
          cons(IM3, k, j, i) = 0.0;
        }

        // Barotropic EOS enforcement (everywhere)
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
      });
  return TaskStatus::complete;
}

} // namespace collapse_be
