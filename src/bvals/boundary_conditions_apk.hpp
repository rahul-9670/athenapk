//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2025, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file boundary_conditions_apk.chpp
//  \brief AthenaPK specific boundary conditions
//

#ifndef BVALS_BOUNDARY_CONDITIONS_APK_HPP_
#define BVALS_BOUNDARY_CONDITIONS_APK_HPP_

#include <memory>
#include <string>
#include <vector>

// Parthenon headers
#include <parthenon/package.hpp>

#include "basic_types.hpp"
#include "bvals/boundary_conditions_generic.hpp"
#include "mesh/domain.hpp"
#include "mesh/mesh.hpp"
#include "mesh/meshblock.hpp"
#include "utils/error_checking.hpp"

#include "../main.hpp"

namespace Hydro {
namespace BoundaryFunction {

using namespace parthenon::package::prelude;
using parthenon::CoordinateDirection;
// using parthenon::MeshBlockData;
// using parthenon::Real;
using parthenon::BoundaryFunction::BCSide;

template <CoordinateDirection DIR, BCSide SIDE>
void ReflectBC(std::shared_ptr<MeshBlockData<Real>> &mbd, bool coarse) {
  // make sure DIR is X[123]DIR so we don't have to check again
  static_assert(DIR == X1DIR || DIR == X2DIR || DIR == X3DIR, "DIR must be X[123]DIR");

  MeshBlock *pmb = mbd->GetBlockPointer();

  auto hydro_pkg = pmb->packages.Get("Hydro");
  auto fluid = hydro_pkg->Param<Fluid>("fluid");
  PARTHENON_REQUIRE_THROWS(
      fluid == Fluid::euler,
      "Reflecting boundary conditions for MHD need special treatment.");

  // convenient shorthands
  constexpr bool X1 = (DIR == X1DIR);
  constexpr bool X2 = (DIR == X2DIR);
  constexpr bool X3 = (DIR == X3DIR);
  constexpr bool INNER = (SIDE == BCSide::Inner);

  const auto &bounds = coarse ? pmb->c_cellbounds : pmb->cellbounds;

  const auto &range = X1 ? bounds.GetBoundsI(IndexDomain::interior)
                         : (X2 ? bounds.GetBoundsJ(IndexDomain::interior)
                               : bounds.GetBoundsK(IndexDomain::interior));
  const int ref = INNER ? range.s : range.e;

  constexpr IndexDomain domain =
      INNER ? (X1 ? IndexDomain::inner_x1
                  : (X2 ? IndexDomain::inner_x2 : IndexDomain::inner_x3))
            : (X1 ? IndexDomain::outer_x1
                  : (X2 ? IndexDomain::outer_x2 : IndexDomain::outer_x3));

  // used for reflections
  const int offset = (2 * ref) + (INNER ? -1 : 1);

  auto cons = mbd->PackVariables(std::vector<std::string>{"cons"}, coarse);
  const bool fine = false; // no usage of fine fields in AthenaPK for now

  const auto nv = IndexRange{0, cons.GetDim(4) - 1};
  pmb->par_for_bndry(
      "ReflectBC", nv, domain, parthenon::TopologicalElement::CC, coarse, fine,
      KOKKOS_LAMBDA(const int &v, const int &k, const int &j, const int &i) {
        const bool reflect = v == DIR;
        cons(v, k, j, i) =
            (reflect ? -1.0 : 1.0) *
            cons(v, X3 ? offset - k : k, X2 ? offset - j : j, X1 ? offset - i : i);
      });
}

//----------------------------------------------------------------------------------------
//! \brief "Diode" outflow: zero-gradient copy, but the normal momentum is clamped so that
//!        material may leave the domain and may NOT enter it.
//!
//! VALIDATION B1 (2026-08-02). Parthenon's `outflow` is a plain zero-gradient ghost copy with
//! no inflow suppression, so a ghost cell holding an inward-pointing velocity drives material
//! INTO the box. Measured on the production deck via WP-6's `cons-Mout` (verified against an
//! independent numpy surface integral to 1 part in 1e6): the column is NEGATIVE, about -178 on
//! each of the six faces, i.e. net INFLOW through nominally-outflow boundaries -- which is the
//! source of the +1.25e-3 mass rise that WP-6 flagged as "anomaly 2". The floors were excluded
//! as the source by direct measurement (`cons-nfloor` = 0 on every row).
//!
//! This BC is the standard fix (Athena++ calls the same thing a "diode"): copy as for outflow,
//! then set the normal momentum to zero if it points inward. Tangential momenta, density and
//! energy are untouched, so it is strictly less permissive than `outflow` and reduces to it
//! exactly whenever the flow is already outward.
//!
//! OPT-IN. Registered under the name "diode"; `outflow` keeps Parthenon's semantics unchanged.
//! A deck that does not say `diode` is bit-identical. This is deliberate: switching the
//! production BC is RESULT-CHANGING and would invalidate the WP-7/8/18 re-baseline that was
//! just established, so the switch is the user's call, made once, explicitly.
//!
//! NOTE this operates on `cons`, so the clamped quantity is momentum (rho*v); zeroing it is
//! equivalent to zeroing the velocity at fixed density.
template <CoordinateDirection DIR, BCSide SIDE>
void DiodeBC(std::shared_ptr<MeshBlockData<Real>> &mbd, bool coarse) {
  static_assert(DIR == X1DIR || DIR == X2DIR || DIR == X3DIR, "DIR must be X[123]DIR");

  MeshBlock *pmb = mbd->GetBlockPointer();

  constexpr bool X1 = (DIR == X1DIR);
  constexpr bool X2 = (DIR == X2DIR);
  constexpr bool X3 = (DIR == X3DIR);
  constexpr bool INNER = (SIDE == BCSide::Inner);

  const auto &bounds = coarse ? pmb->c_cellbounds : pmb->cellbounds;
  const auto &range = X1 ? bounds.GetBoundsI(IndexDomain::interior)
                         : (X2 ? bounds.GetBoundsJ(IndexDomain::interior)
                               : bounds.GetBoundsK(IndexDomain::interior));
  const int ref = INNER ? range.s : range.e;

  constexpr IndexDomain domain =
      INNER ? (X1 ? IndexDomain::inner_x1
                  : (X2 ? IndexDomain::inner_x2 : IndexDomain::inner_x3))
            : (X1 ? IndexDomain::outer_x1
                  : (X2 ? IndexDomain::outer_x2 : IndexDomain::outer_x3));

  auto cons = mbd->PackVariables(std::vector<std::string>{"cons"}, coarse);
  const bool fine = false;

  const auto nv = IndexRange{0, cons.GetDim(4) - 1};
  pmb->par_for_bndry(
      "DiodeBC", nv, domain, parthenon::TopologicalElement::CC, coarse, fine,
      KOKKOS_LAMBDA(const int &v, const int &k, const int &j, const int &i) {
        // Zero-gradient copy from the nearest interior cell (identical to outflow).
        const Real val = cons(v, X3 ? ref : k, X2 ? ref : j, X1 ? ref : i);
        // The normal momentum component is IM1/IM2/IM3 == DIR (see main.hpp: IDN=0,
        // IM1=1, IM2=2, IM3=3, so component index v == DIR selects the normal momentum).
        // "Inward" is +normal on an inner face and -normal on an outer face.
        const bool normal_mom = (v == DIR);
        const bool inward = INNER ? (val > 0.0) : (val < 0.0);
        cons(v, k, j, i) = (normal_mom && inward) ? 0.0 : val;
      });
}

} // namespace BoundaryFunction
} // namespace Hydro

#endif // BVALS_BOUNDARY_CONDITIONS_APK_HPP_
