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

// ReflectBC was removed 2026-08-10 with the pure-hydro path. It opened with
// PARTHENON_REQUIRE_THROWS(fluid == Fluid::euler, "Reflecting boundary conditions for
// MHD need special treatment."), so in this MHD-only build it could only ever abort.
// Restore from git tag validation-complete-2026-08-10 if a hydro problem ever needs it.

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

  // WP-13b (2026-08-08). FIRST do the plain zero-gradient copy for EVERY variable, exactly as
  // Parthenon's stock `outflow` does (same GenericBC, same variable_names::any pack); THEN apply
  // the diode momentum clamp to `cons` below.
  //
  // WHY. This function used to pack ONLY "cons", so with `ix1_bc = diode` the domain-boundary
  // ghost zones of every OTHER FillGhost field -- in particular all of the M1 radiation moments
  // rad.Er / rad.Fr1..3 and their per-group copies -- were never written by any boundary
  // condition at all. A FRESH run got away with it because the problem generator initialises the
  // whole array, ghosts included, so ghost ~ interior and the face flux stays ~1e-4. A RESTART
  // does not: restart files store INTERIOR cells only, so the ghosts came up zero-initialised and
  // the first CalculateRadFluxes after the restart saw a full-amplitude jump across the domain
  // face. Measured on the flagship deck (job 2491776, ensemble point000, 4xH100, 31 cycles):
  // restart-vs-fresh differed by 7.41e+01 in rad.Fr1/2/3 against a non-determinism floor of
  // 8.6e-09 (1e10x), confined EXACTLY to the outermost meshblock layer, and predicted by
  // chat*Er/2 = 1578*0.104/2 ~ 79. With the radiation package off the same test passes at the
  // floor, and nsub, block ordering and the restart file were all separately exonerated.
  // This reaches the flagship and all 24 ensemble members: every one uses `diode` + radiation.
  //
  // The cons path is UNCHANGED -- the same copy from the same reference cell, then the same
  // clamp -- so a radiation-free run is bit-identical to before this fix.
  parthenon::BoundaryFunction::GenericBC<DIR, SIDE, parthenon::BoundaryFunction::BCType::Outflow,
                                         parthenon::variable_names::any>(mbd, coarse);

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
