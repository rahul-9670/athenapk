// AthenaPK - a performance portable block structured AMR MHD code
// Copyright (c) 2020-2024, Athena Parthenon Collaboration. All rights reserved.
// Licensed under the 3-Clause License (the "LICENSE")

#ifndef MAIN_HPP_
#define MAIN_HPP_

#include <limits> // numeric limits

// Parthenon headers
#include "basic_types.hpp" // Real
#include <parthenon/package.hpp>

// TODO(pgrete) There's a compiler bug in nvcc < 11.2 that precludes the use
// of C++17 with relaxed-constexpr in Kokkos,
// see https://github.com/kokkos/kokkos/issues/3496
// This also precludes our downstream use of constexpr int here.
// Update once nvcc/cuda >= 11.2 is more widely available on machine.
enum {
  IDN = 0,
  IM1 = 1,
  IM2 = 2,
  IM3 = 3,
  IEN = 4,
  NHYDRO = 5,
  IB1 = 5,
  IB2 = 6,
  IB3 = 7,
  IPS = 8
};

// array indices for 1D primitives: velocity, transverse components of field
enum { IV1 = 1, IV2 = 2, IV3 = 3, IPR = 4 };

// 2026-08-10, flagship reduction: `hllc` and `Fluid::euler` are gone. They were the
// pure-hydro path -- the flagship is MHD (glmmhd + CT) and no production deck ever set
// `hydro/fluid = euler`. Removing euler also removed the three hydro-only Riemann solvers
// (hydro_hlle/hydro_hllc/hydro_dc_llf) and ~14 template instantiations of the flux kernel.
// `Cooling` went with the tabular optically-thin cooling, which was cluster/ISM physics;
// the flagship cools barotropically and through M1 radiation.
// All of it is restorable from git tag `validation-complete-2026-08-10`.
enum class RiemannSolver { undefined, none, hlle, llf, hlld };
enum class Reconstruction { undefined, dc, plm, ppm, wenoz, weno3, limo3 };
enum class Integrator { undefined, rk1, rk2, vl2, rk3 };
enum class Fluid { undefined, glmmhd };
enum class Conduction { none, isotropic, anisotropic };
enum class ConductionCoeff { none, fixed, spitzer };
enum class Viscosity { none, isotropic };
enum class ViscosityCoeff { none, fixed };
enum class Resistivity { none, ohmic };
// ionization = equilibrium self-consistent x_e(rho,T); ionization_chem = use the
// time-dependent x_e evolved by the chemistry package (gow17_reduced) in the Wardle tensor.
enum class ResistivityCoeff { none, fixed, spitzer, ionization, ionization_chem };
enum class Ambipolar { none, ambipolar };
// ionization = equilibrium self-consistent x_e(rho,T); ionization_chem = use the
// time-dependent x_e evolved by the chemistry package (gow17_reduced) for rho_i.
enum class AmbipolarCoeff { none, fixed, ionization, ionization_chem };
enum class Hall { none, hall };
// ionization = equilibrium x_e(rho,T); ionization_chem = chemistry-evolved x_e in the tensor.
enum class HallCoeff { none, fixed, ionization, ionization_chem };
enum class DiffInt { none, unsplit, rkl2 };

enum class Hst { idx, ekin, emag, divb, divb_max };

enum class CartesianDir { x, y, z };

constexpr parthenon::Real float_min{std::numeric_limits<float>::min()};

using InitPackageDataFun_t =
    std::function<void(parthenon::ParameterInput *pin, parthenon::StateDescriptor *pkg)>;

#endif // MAIN_HPP_
