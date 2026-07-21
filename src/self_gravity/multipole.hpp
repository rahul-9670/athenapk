//========================================================================================
// AthenaPK - Self-gravity package
// Multipole (Cartesian, through quadrupole) exterior boundary condition for phi.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
#ifndef SELF_GRAVITY_MULTIPOLE_HPP_
#define SELF_GRAVITY_MULTIPOLE_HPP_

// WS-5a of PHYSICS_COMPLETION_PLAN.md. Replaces the phi=0 Dirichlet walls with the
// exterior Cartesian multipole expansion of the interior mass (monopole + traceless
// quadrupole about the center of mass), so an isolated collapsing cloud can live in a
// box only a few times its own size without the BC clipping the potential well.
//
// Units: code units, four_pi_G = 1 in production. The potential prefactor is
//   G = four_pi_G / (4 pi).
// The dipole term is identically zero about the center of mass and is omitted.
//
// Ghost fill uses FixedFace semantics (matches the zero-Dirichlet path in
// poisson_equation.hpp): the domain-boundary FACE value is pinned to the expansion,
//   ghost = 2 * Phi_mp(face) - interior_mirror.

#include <cmath>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <bvals/boundary_conditions_generic.hpp>
#include <coordinates/coordinates.hpp>
#include <parthenon/package.hpp>

#include "self_gravity.hpp"

namespace SelfGravity {

// POD carried in a mutable package param, recomputed once per step in FillPoissonRHS
// and captured by value into the (device) boundary kernels. Q is the traceless
// symmetric quadrupole about the center of mass (cx,cy,cz), stored xx,yy,zz,xy,xz,yz.
struct MultipoleMoments {
  Real M;
  Real cx, cy, cz;
  Real Qxx, Qyy, Qzz, Qxy, Qxz, Qyz;
  Real four_pi_G;
};

// Exterior potential at physical point (x,y,z). Valid outside the mass distribution.
KOKKOS_INLINE_FUNCTION Real MultipolePhi(const MultipoleMoments &mm, Real x, Real y,
                                         Real z) {
  constexpr Real kInv4Pi = 0.25 / M_PI;
  const Real dx = x - mm.cx;
  const Real dy = y - mm.cy;
  const Real dz = z - mm.cz;
  const Real r2 = dx * dx + dy * dy + dz * dz;
  const Real r = std::sqrt(r2);
  if (r <= 0.0) return 0.0; // degenerate; a domain boundary is never at the COM
  const Real invr = 1.0 / r;
  const Real invr5 = invr / (r2 * r2);
  const Real dQd = mm.Qxx * dx * dx + mm.Qyy * dy * dy + mm.Qzz * dz * dz +
                   2.0 * (mm.Qxy * dx * dy + mm.Qxz * dx * dz + mm.Qyz * dy * dz);
  const Real G = mm.four_pi_G * kInv4Pi;
  return -G * (mm.M * invr + 0.5 * dQd * invr5);
}

// Per-block physical-boundary function for phi enrolled as a UserBoundaryFunction (the
// packed solver-internal fast path lives in PoissonEquation::SetBoundary). Mirrors
// parthenon's GenericBC<...,FixedFace,...> loop but with a spatially varying face value
// = MultipolePhi(boundary face point). Cell-centered (phi) only.
template <parthenon::CoordinateDirection DIR, parthenon::BoundaryFunction::BCSide SIDE>
void MultipoleBC(std::shared_ptr<parthenon::MeshBlockData<Real>> &rc, bool coarse) {
  using namespace parthenon;
  using TE = TopologicalElement;
  constexpr bool X1 = (DIR == X1DIR);
  constexpr bool X2 = (DIR == X2DIR);
  constexpr bool X3 = (DIR == X3DIR);
  constexpr bool INNER = (SIDE == BoundaryFunction::BCSide::Inner);
  constexpr auto el = TE::CC;

  MeshBlock *pmb = rc->GetBlockPointer();
  // The coarse BC is only invoked on blocks with coarser neighbors, so pmr is valid
  // there (boundary_conditions.cpp:46). The pack always carries FINE coords, so we
  // fetch the matching coordinate object explicitly.
  const Coordinates_t coords =
      (coarse && pmb->pmr) ? pmb->pmr->GetCoarseCoords() : pmb->coords;

  auto pkg = pmb->packages.Get("self_gravity");
  const MultipoleMoments mm = pkg->Param<MultipoleMoments>("grav_multipole_moments");

  std::set<PDOpt> opts = coarse ? std::set<PDOpt>{PDOpt::Coarse} : std::set<PDOpt>{};
  std::vector<MetadataFlag> flags{Metadata::FillGhost, Metadata::Cell};
  const auto desc = MakePackDescriptor<grav::phi>(rc.get(), flags, opts);
  auto q = desc.GetPack(rc.get());
  const int b = 0;
  const int lstart = q.GetLowerBoundHost(b);
  const int lend = q.GetUpperBoundHost(b);
  if (lend < lstart) return; // phi not allocated on this block

  const auto &bounds = coarse ? pmb->c_cellbounds : pmb->cellbounds;
  const auto range = X1 ? bounds.GetBoundsI(IndexDomain::interior, el)
                        : (X2 ? bounds.GetBoundsJ(IndexDomain::interior, el)
                              : bounds.GetBoundsK(IndexDomain::interior, el));
  const int ref = INNER ? range.s : range.e;
  const int offset = 2 * ref + (INNER ? -1 : 1);
  // Lower-face index (along DIR) of the domain boundary face.
  const int faceidx = INNER ? ref : ref + 1;

  constexpr IndexDomain domain =
      INNER ? (X1 ? IndexDomain::inner_x1
                  : (X2 ? IndexDomain::inner_x2 : IndexDomain::inner_x3))
            : (X1 ? IndexDomain::outer_x1
                  : (X2 ? IndexDomain::outer_x2 : IndexDomain::outer_x3));

  auto nb = IndexRange{lstart, lend};
  pmb->par_for_bndry(
      "SG::MultipoleBC", nb, domain, el, coarse, /*fine=*/false,
      KOKKOS_LAMBDA(const int &l, const int &k, const int &j, const int &i) {
        // Boundary-face point aligned with this ghost column: normal coordinate at the
        // domain face, transverse coordinates at the ghost cell center.
        const Real px = X1 ? coords.Xf<X1DIR>(k, j, faceidx) : coords.Xc<X1DIR>(k, j, i);
        const Real py = X2 ? coords.Xf<X2DIR>(k, faceidx, i) : coords.Xc<X2DIR>(k, j, i);
        const Real pz = X3 ? coords.Xf<X3DIR>(faceidx, j, i) : coords.Xc<X3DIR>(k, j, i);
        const Real val = MultipolePhi(mm, px, py, pz);
        q(b, el, l, k, j, i) =
            2.0 * val - q(b, el, l, X3 ? offset - k : k, X2 ? offset - j : j,
                          X1 ? offset - i : i);
      });
}

} // namespace SelfGravity

#endif // SELF_GRAVITY_MULTIPOLE_HPP_
