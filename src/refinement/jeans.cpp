//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD
// code. Copyright (c) 2025, Athena-Parthenon Collaboration. All rights
// reserved. Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include <algorithm>
#include <cmath>
#include <limits>

#include "../hydro/hydro.hpp"
#include "../main.hpp"
#include "refinement.hpp"

namespace refinement {
namespace jeans {

using parthenon::IndexDomain;
using parthenon::IndexRange;

// Jeans-length refinement criterion (Truelove 1997).
// Refines when the Jeans length is resolved by fewer than `njeans` cells.
//
// lambda_J = c_s * sqrt(pi / (G rho)) = 2*pi * c_s / sqrt(four_pi_G * rho), i.e.
// in code units with 4*pi*G = 1:
//   lambda_J = 2*pi * c_s / sqrt(rho)    (hydro)
//   lambda_J = 2*pi * (c_s + v_A) / sqrt(rho)   (MHD, Athena++ convention)
// where v_A = sqrt(B^2 / rho) is the Alfven speed (Heaviside-Lorentz, no 4 pi).
// four_pi_G is taken from the self_gravity package when it is active, so the
// criterion stays correct for a normalization other than 4*pi*G = 1.
//
// nj = lambda_J / dx is "Jeans length in cells"
//   nj < njeans         -> refine
//   nj > 2.5 * njeans   -> derefine
parthenon::AmrTag Jeans(MeshBlockData<Real> *rc) {
  auto pmb = rc->GetBlockPointer();
  auto w = rc->Get("prim").data;

  auto hydro_pkg = pmb->packages.Get("Hydro");
  const Real njeans = hydro_pkg->Param<Real>("refinement/njeans");
  const Real gam = pmb->packages.Get("Hydro")->Param<Real>("AdiabaticIndex");
  const bool mhd = (hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd);

  // Use the LARGEST cell spacing so the criterion is conservative on non-cubic
  // cells (for the usual cubic cells this is just dx).
  const Real dx =
      std::max({pmb->coords.Dxc<1>(0), pmb->coords.Dxc<2>(0), pmb->coords.Dxc<3>(0)});

  // 4 pi G in code units; 1 when self-gravity is not active (external potential).
  const auto &pkgs = pmb->packages.AllPackages();
  const Real four_pi_G = pkgs.count("self_gravity") > 0
                             ? pkgs.at("self_gravity")->Param<Real>("four_pi_G")
                             : 1.0;

  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  const Real fac = 2.0 * M_PI / (dx * std::sqrt(four_pi_G));

  Real njmin = std::numeric_limits<Real>::max();
  pmb->par_reduce(
      "jeans refinement", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int k, const int j, const int i, Real &lnjmin) {
        const Real rho = w(IDN, k, j, i);
        const Real p = w(IPR, k, j, i);
        const Real cs = std::sqrt(gam * p / rho);
        Real v = cs;
        if (mhd) {
          const Real bsq = w(IB1, k, j, i) * w(IB1, k, j, i) +
                           w(IB2, k, j, i) * w(IB2, k, j, i) +
                           w(IB3, k, j, i) * w(IB3, k, j, i);
          const Real va = std::sqrt(bsq / rho);
          v += va;
        }
        const Real nj = v / std::sqrt(rho);
        lnjmin = std::min(lnjmin, nj);
      },
      Kokkos::Min<Real>(njmin));

  njmin *= fac;

  if (njmin < njeans) return parthenon::AmrTag::refine;
  if (njmin > 2.5 * njeans) return parthenon::AmrTag::derefine;
  return parthenon::AmrTag::same;
}

} // namespace jeans
} // namespace refinement
