//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2024, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-Clause License (the "LICENSE").
//========================================================================================
//! \file cshock.cpp
//! \brief 1D oblique MHD shock tube for validating ambipolar diffusion (C-shock test).
//!
//! Two uniform left/right states are separated by a discontinuity at x = xshock. The
//! magnetic field is oblique to the shock normal (Bx is uniform, as required by div(B)=0
//! in 1D; By/Bz may jump). In the ideal limit this relaxes to the usual MHD shock
//! structure; with ambipolar diffusion enabled (diffusion/ambipolar=ambipolar,
//! diffusion/integrator=unsplit) the magnetic precursor smooths the transition into a
//! continuous C-type shock whose width scales with the ambipolar coefficient. Run to a
//! steady state and compare the relaxed profile to the (semi-)analytic C-shock solution.
//!
//! All states are read from the <problem/cshock> block; defaults form a strong fast shock.
//! For a near-isothermal C-shock, run with a small adiabatic index (e.g. hydro/gamma~1.01).

// C++ headers
#include <cmath>

// Parthenon headers
#include "mesh/mesh.hpp"
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

// AthenaPK headers
#include "../main.hpp"
#include "utils/error_checking.hpp"

namespace cshock {
using namespace parthenon::driver::prelude;

void ProblemGenerator(MeshBlock *pmb, ParameterInput *pin) {
  auto hydro_pkg = pmb->packages.Get("Hydro");
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  auto &mbd = pmb->meshblock_data.Get();
  auto &u = mbd->Get("cons").data;

  const auto gamma = pin->GetReal("hydro", "gamma");
  const bool mhd_enabled = hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd;
  PARTHENON_REQUIRE_THROWS(mhd_enabled,
                           "The cshock problem requires MHD (hydro/fluid=glmmhd).");

  // Location of the initial discontinuity.
  const auto xshock = pin->GetOrAddReal("problem/cshock", "xshock", 0.0);
  // Uniform normal field component (constant for div(B)=0 in 1D).
  const auto bx = pin->GetOrAddReal("problem/cshock", "bx", 1.0);

  // Left (upstream) primitive state. Defaults: supersonic inflow toward the shock.
  const auto dl = pin->GetOrAddReal("problem/cshock", "dl", 1.0);
  const auto vxl = pin->GetOrAddReal("problem/cshock", "vxl", 5.0);
  const auto vyl = pin->GetOrAddReal("problem/cshock", "vyl", 0.0);
  const auto vzl = pin->GetOrAddReal("problem/cshock", "vzl", 0.0);
  const auto byl = pin->GetOrAddReal("problem/cshock", "byl", 1.0);
  const auto bzl = pin->GetOrAddReal("problem/cshock", "bzl", 0.0);
  const auto pl = pin->GetOrAddReal("problem/cshock", "pl", 1.0);

  // Right (downstream) primitive state. Defaults: compressed, decelerated post-shock gas.
  const auto dr = pin->GetOrAddReal("problem/cshock", "dr", 3.0);
  const auto vxr = pin->GetOrAddReal("problem/cshock", "vxr", 5.0 / 3.0);
  const auto vyr = pin->GetOrAddReal("problem/cshock", "vyr", 0.0);
  const auto vzr = pin->GetOrAddReal("problem/cshock", "vzr", 0.0);
  const auto byr = pin->GetOrAddReal("problem/cshock", "byr", 3.0);
  const auto bzr = pin->GetOrAddReal("problem/cshock", "bzr", 0.0);
  const auto pr = pin->GetOrAddReal("problem/cshock", "pr", 16.0);

  const Real gm1 = gamma - 1.0;
  auto &coords = pmb->coords;

  pmb->par_for(
      "ProblemGenerator: cshock", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        const bool left = coords.Xc<1>(i) <= xshock;

        const Real d = left ? dl : dr;
        const Real vx = left ? vxl : vxr;
        const Real vy = left ? vyl : vyr;
        const Real vz = left ? vzl : vzr;
        const Real by = left ? byl : byr;
        const Real bz = left ? bzl : bzr;
        const Real p = left ? pl : pr;

        u(IDN, k, j, i) = d;
        u(IM1, k, j, i) = d * vx;
        u(IM2, k, j, i) = d * vy;
        u(IM3, k, j, i) = d * vz;

        u(IB1, k, j, i) = bx;
        u(IB2, k, j, i) = by;
        u(IB3, k, j, i) = bz;

        u(IEN, k, j, i) = p / gm1 + 0.5 * d * (SQR(vx) + SQR(vy) + SQR(vz)) +
                          0.5 * (SQR(bx) + SQR(by) + SQR(bz));
      });
}
} // namespace cshock
