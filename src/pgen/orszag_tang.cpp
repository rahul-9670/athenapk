
//========================================================================================
// AthenaPK - a performance portable block structured AMR MHD code
// Copyright (c) 2021-2023, Athena Parthenon Collaboration. All rights reserved.
// Licensed under the 3-Clause License (the "LICENSE")
//========================================================================================
//! \file orszag_tang.cpp
//! \brief Problem generator for the Orszag Tang vortex.
//!
//! REFERENCE: Orszag & Tang (J. Fluid Mech., 90, 129, 1998) and
//! https://www.astro.princeton.edu/~jstone/Athena/tests/orszag-tang/pagesource.html
//========================================================================================

// Parthenon headers
#include "mesh/mesh.hpp"
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

// AthenaPK headers
#include "../hydro/ct/ct.hpp"
#include "../main.hpp"

namespace orszag_tang {
using namespace parthenon::driver::prelude;

void ProblemGenerator(MeshBlock *pmb, ParameterInput *pin) {
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  auto &mbd = pmb->meshblock_data.Get();
  auto &u = mbd->Get("cons").data;

  Real gm1 = pin->GetReal("hydro", "gamma") - 1.0;
  Real B0 = 1.0 / std::sqrt(4.0 * M_PI);
  Real d0 = 25.0 / (36.0 * M_PI);
  Real v0 = 1.0;
  Real p0 = 5.0 / (12.0 * M_PI);

  auto &coords = pmb->coords;

  pmb->par_for(
      "ProblemGenerator: Orszag-Tang", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        u(IDN, k, j, i) = d0;
        // Note the different signs in this pgen compared to the the eqn mentioned in the
        // original paper (and other codes).
        // They are related to our domain going from -0.5 to 0.5 (for symmetry reason)
        // rather than 0  to 2pi (i.e., the sign for single wave sinus is flipped).
        u(IM1, k, j, i) = d0 * v0 * std::sin(2.0 * M_PI * coords.Xc<2>(j));
        u(IM2, k, j, i) = -d0 * v0 * std::sin(2.0 * M_PI * coords.Xc<1>(i));
        u(IM3, k, j, i) = 0.0;

        u(IB1, k, j, i) = B0 * std::sin(2.0 * M_PI * coords.Xc<2>(j));
        u(IB2, k, j, i) = B0 * std::sin(4.0 * M_PI * coords.Xc<1>(i));
        u(IB3, k, j, i) = 0.0;

        u(IEN, k, j, i) =
            p0 / gm1 +
            0.5 * (SQR(u(IB1, k, j, i)) + SQR(u(IB2, k, j, i)) + SQR(u(IB3, k, j, i)) +
                   (SQR(u(IM1, k, j, i)) + SQR(u(IM2, k, j, i)) + SQR(u(IM3, k, j, i))) /
                       u(IDN, k, j, i));
      });

  // ---- Constrained Transport: face field Bf from the analytic vector potential ----
  // A_z(x,y) = B0/(4 pi) cos(4 pi x) - B0/(2 pi) cos(2 pi y); B = curl(A_z e_z) recovers
  // the cell-centered field above. Discrete curl of nodal A_z => div(B)_face = 0.
  auto hydro_pkg = pmb->packages.Get("Hydro");
  if (hydro_pkg->Param<bool>("use_ct")) {
    using TE = parthenon::TopologicalElement;
    const Real B0_l = B0;
    auto desc = parthenon::MakePackDescriptor<Hydro::CT::Bf>(mbd.get());
    auto pack = desc.GetPack(mbd.get());
    const int bidx = 0;
    auto Az = KOKKOS_LAMBDA(const Real x, const Real y)->Real {
      return B0_l / (4.0 * M_PI) * std::cos(4.0 * M_PI * x) -
             B0_l / (2.0 * M_PI) * std::cos(2.0 * M_PI * y);
    };
    IndexRange fib = pmb->cellbounds.GetBoundsI(IndexDomain::interior, TE::F1);
    IndexRange fjb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior, TE::F1);
    IndexRange fkb = pmb->cellbounds.GetBoundsK(IndexDomain::interior, TE::F1);
    pmb->par_for(
        "OT CT F1", fkb.s, fkb.e, fjb.s, fjb.e, fib.s, fib.e,
        KOKKOS_LAMBDA(const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(bidx);
          const Real x = c.X<X1DIR, TE::F1>(k, j, i);
          const Real ylo = c.X<X2DIR, TE::F2>(k, j, i);
          const Real yhi = c.X<X2DIR, TE::F2>(k, j + 1, i);
          pack(bidx, TE::F1, Hydro::CT::Bf(), k, j, i) =
              (Az(x, yhi) - Az(x, ylo)) / (yhi - ylo);
        });
    fib = pmb->cellbounds.GetBoundsI(IndexDomain::interior, TE::F2);
    fjb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior, TE::F2);
    fkb = pmb->cellbounds.GetBoundsK(IndexDomain::interior, TE::F2);
    pmb->par_for(
        "OT CT F2", fkb.s, fkb.e, fjb.s, fjb.e, fib.s, fib.e,
        KOKKOS_LAMBDA(const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(bidx);
          const Real y = c.X<X2DIR, TE::F2>(k, j, i);
          const Real xlo = c.X<X1DIR, TE::F1>(k, j, i);
          const Real xhi = c.X<X1DIR, TE::F1>(k, j, i + 1);
          pack(bidx, TE::F2, Hydro::CT::Bf(), k, j, i) =
              -(Az(xhi, y) - Az(xlo, y)) / (xhi - xlo);
        });
    fib = pmb->cellbounds.GetBoundsI(IndexDomain::interior, TE::F3);
    fjb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior, TE::F3);
    fkb = pmb->cellbounds.GetBoundsK(IndexDomain::interior, TE::F3);
    pmb->par_for(
        "OT CT F3", fkb.s, fkb.e, fjb.s, fjb.e, fib.s, fib.e,
        KOKKOS_LAMBDA(const int k, const int j, const int i) {
          pack(bidx, TE::F3, Hydro::CT::Bf(), k, j, i) = 0.0;
        });
  }
}
} // namespace orszag_tang
