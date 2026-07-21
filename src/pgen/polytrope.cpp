//========================================================================================
// AthenaPK - n=1 polytrope (Lane-Emden) self-gravity benchmark problem generator
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
//! \file polytrope.cpp
//  \brief Static n=1 polytrope sphere for an APPLES-TO-APPLES GPU throughput
//         comparison against artemis' `polytrope` pgen (artemis/src/pgen/polytrope.hpp).
//
//  Replicates artemis iprob=1 EXACTLY so only the codebase differs:
//      alpha   = sqrt(1/2)
//      ar      = alpha * |x - x0|         (x0 = origin)
//      theta   = sin(ar)/ar               (Lane-Emden n=1 profile)
//      cutoff  = 0.75*pi
//      inside  = (ar < cutoff)
//      rho     = inside ? theta : rho_amb         (rho_amb = 1e-3)
//      sie     = inside ? theta : sie_amb         (sie_amb = 1e2)   [specific internal energy]
//      v       = 0
//  artemis stores (rho, sie); AthenaPK stores conserved (rho, mom, Etot). With v=0 the
//  internal energy density is rho*sie, so for an ideal EOS P = (gamma-1)*rho*sie and
//      IEN = P/(gamma-1) + KE = rho*sie.
//  For the matched benchmark use gamma=2.0 + four_pi_G=1 (the n=1 polytrope is the
//  P = K rho^2 hydrostatic sphere); this mirrors artemis' poly.in.

#include <cmath>
#include <string>

#include <parthenon/package.hpp>

#include "../main.hpp"
#include "pgen.hpp"

namespace polytrope {
using namespace parthenon::driver::prelude;

void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
  auto hydro_pkg = pmb->packages.Get("Hydro");
  const bool mhd = (hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd);

  // Polytrope params (artemis-matched defaults)
  const Real rho_amb = pin->GetOrAddReal("problem/polytrope", "rho_amb", 1.0e-3);
  const Real sie_amb = pin->GetOrAddReal("problem/polytrope", "sie_amb", 1.0e2);
  const Real x0 = pin->GetOrAddReal("problem/polytrope", "x10a", 0.0);
  const Real y0 = pin->GetOrAddReal("problem/polytrope", "x20a", 0.0);
  const Real z0 = pin->GetOrAddReal("problem/polytrope", "x30a", 0.0);

  // n=1 Lane-Emden hardcoded params (artemis polytrope.hpp:72-73)
  const Real alpha = std::sqrt(0.5);
  const Real cutoff = 0.75 * M_PI;

  auto &data = pmb->meshblock_data.Get();
  auto &u_dev = data->Get("cons").data;
  auto u = u_dev.GetHostMirrorAndCopy();

  parthenon::IndexRange ib = pmb->cellbounds.GetBoundsI(parthenon::IndexDomain::entire);
  parthenon::IndexRange jb = pmb->cellbounds.GetBoundsJ(parthenon::IndexDomain::entire);
  parthenon::IndexRange kb = pmb->cellbounds.GetBoundsK(parthenon::IndexDomain::entire);

  auto &coords = pmb->coords;
  for (int k = kb.s; k <= kb.e; ++k) {
    for (int j = jb.s; j <= jb.e; ++j) {
      for (int i = ib.s; i <= ib.e; ++i) {
        const Real dx = coords.Xc<1>(i) - x0;
        const Real dy = coords.Xc<2>(j) - y0;
        const Real dz = coords.Xc<3>(k) - z0;
        const Real ar = alpha * std::sqrt(dx * dx + dy * dy + dz * dz);
        const Real theta = (ar > 0.0) ? std::sin(ar) / ar : 1.0;  // sin(ar)/ar -> 1
        const bool inside = (ar < cutoff);

        const Real rho = inside ? theta : rho_amb;
        const Real sie = inside ? theta : sie_amb;
        // internal energy density u_int = rho * sie (v=0); IEN = u_int + KE = rho*sie
        const Real eint = rho * sie;

        u(IDN, k, j, i) = rho;
        u(IM1, k, j, i) = 0.0;
        u(IM2, k, j, i) = 0.0;
        u(IM3, k, j, i) = 0.0;
        u(IEN, k, j, i) = eint;
        if (mhd) {
          u(IB1, k, j, i) = 0.0;
          u(IB2, k, j, i) = 0.0;
          u(IB3, k, j, i) = 0.0;
          u(IPS, k, j, i) = 0.0;
        }
      }
    }
  }
  u_dev.DeepCopy(u);
}

} // namespace polytrope
