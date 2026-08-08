//========================================================================================
// AthenaPK - magnetic-transport (fossil-field) history diagnostics. See mag_diag.hpp.
//========================================================================================

#include <cmath>
#include <string>
#include <vector>

#include <parthenon/package.hpp>

#include "../hydro/diffusion/diffusion.hpp"
#include "../main.hpp"
#include "mag_diag.hpp"

namespace Diagnostics {

Real MagDiagReduce(MeshData<Real> *md, MagDiag which) {
  const auto &prim = md->PackVariables(std::vector<std::string>{"prim"});

  const bool need_eta = (which == MagDiag::dissO || which == MagDiag::dissA ||
                         which == MagDiag::dissOcap || which == MagDiag::dissAcap ||
                         which == MagDiag::dissOhi || which == MagDiag::dissOlo ||
                         which == MagDiag::dissAhi || which == MagDiag::dissAlo ||
                         which == MagDiag::dissOsq || which == MagDiag::dissAsq);
  const bool need_cap = (which == MagDiag::dissOcap || which == MagDiag::dissAcap);
  // WP-8 density split. Read here rather than passed in so every call site stays unchanged.
  // 0 (the default) means the split columns are never registered, so this is inert then.
  const Real rho_split =
      md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro")->Param<Real>(
          "mag_diag_rho_split");
  // WP-8 current-sheet split (2026-08-08). Dimensionless threshold on s = |J|*dx/|B|; see the
  // header. 0 (the default) means these columns are never registered, so this is inert then.
  const Real sheet_thresh =
      md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro")->Param<Real>(
          "mag_diag_sheet_thresh");
  // Dummy-alias to `prim` when unused: PackVariables on a field that is not registered
  // would throw, and the lambda must still capture *something*.
  const auto &eta =
      need_eta ? md->PackVariables(std::vector<std::string>{"nonideal_eta"}) : prim;
  const auto &capd =
      need_cap ? md->PackVariables(std::vector<std::string>{"diff.capdiag"}) : prim;

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real sum = 0.0;
  Kokkos::parallel_reduce(
      "MagDiagReduce",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExecSpace(), {0, kb.s, jb.s, ib.s},
          {prim.GetDim(5), kb.e + 1, jb.e + 1, ib.e + 1}, {1, 1, 1, ib.e + 1 - ib.s}),
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &lsum) {
        const auto &w = prim(b);
        const auto &coords = prim.GetCoords(b);
        const Real dV = coords.CellVolume(k, j, i);
        const Real bx = w(IB1, k, j, i), by = w(IB2, k, j, i), bz = w(IB3, k, j, i);

        if (which == MagDiag::MEtor || which == MagDiag::MEpol) {
          // Toroidal direction about the z axis: phihat = (-y, x, 0)/R. The z axis is both
          // the rotation axis (omegatff) and the initial field direction (B0z), so this is
          // the physically meaningful split for the FHC problem. On the axis (R -> 0) the
          // toroidal component is undefined; there B is essentially poloidal, so assign
          // the whole cell to poloidal.
          const Real x = coords.Xc<1>(i), y = coords.Xc<2>(j);
          const Real R2 = x * x + y * y;
          Real b_tor2 = 0.0;
          if (R2 > 0.0) {
            const Real b_tor = (-y * bx + x * by) / std::sqrt(R2);
            b_tor2 = b_tor * b_tor;
          }
          const Real b2 = bx * bx + by * by + bz * bz;
          lsum += 0.5 * ((which == MagDiag::MEtor) ? b_tor2 : (b2 - b_tor2)) * dV;
          return;
        }

        // J = curl B, central differences on cell-centered primitive B. Same stencil as the
        // jeans_nonideal current-sheet criterion (prim is FillGhost, so ghosts are valid).
        const Real dx = coords.Dxc<1>(k, j, i);
        const Real dy = coords.Dxc<2>(k, j, i);
        const Real dz = coords.Dxc<3>(k, j, i);
        const Real Jx = (w(IB3, k, j + 1, i) - w(IB3, k, j - 1, i)) / (2.0 * dy) -
                        (w(IB2, k + 1, j, i) - w(IB2, k - 1, j, i)) / (2.0 * dz);
        const Real Jy = (w(IB1, k + 1, j, i) - w(IB1, k - 1, j, i)) / (2.0 * dz) -
                        (w(IB3, k, j, i + 1) - w(IB3, k, j, i - 1)) / (2.0 * dx);
        const Real Jz = (w(IB2, k, j, i + 1) - w(IB2, k, j, i - 1)) / (2.0 * dx) -
                        (w(IB1, k, j + 1, i) - w(IB1, k, j - 1, i)) / (2.0 * dy);
        const Real Jsq = Jx * Jx + Jy * Jy + Jz * Jz;

        // --- WP-8 density-split + carrying-volume variants -------------------------------
        // Each is an integral over a region defined by PHYSICS (a density threshold) rather
        // than by the grid, so the core and envelope budgets can be tracked separately
        // instead of being summed into one ill-conditioned scalar. See mag_diag.hpp.
        if (which == MagDiag::Vhi) {
          if (w(IDN, k, j, i) > rho_split) lsum += dV;
          return;
        }
        if (which == MagDiag::dissOhi || which == MagDiag::dissOlo ||
            which == MagDiag::dissOsq) {
          const Real q = eta(b, NonidealEtaIdx::O, k, j, i) * Jsq;
          if (which == MagDiag::dissOsq) {
            lsum += q * q * dV;
          } else {
            const bool hi = w(IDN, k, j, i) > rho_split;
            if (hi == (which == MagDiag::dissOhi)) lsum += q * dV;
          }
          return;
        }
        if (which == MagDiag::dissAhi || which == MagDiag::dissAlo ||
            which == MagDiag::dissAsq) {
          const Real b2s = bx * bx + by * by + bz * bz;
          Real Jp2 = Jsq;
          if (b2s > 0.0) {
            const Real Jpar = (Jx * bx + Jy * by + Jz * bz) / std::sqrt(b2s);
            Jp2 = Jsq - Jpar * Jpar;
            if (Jp2 < 0.0) Jp2 = 0.0;
          }
          const Real q = eta(b, NonidealEtaIdx::A, k, j, i) * Jp2;
          if (which == MagDiag::dissAsq) {
            lsum += q * q * dV;
          } else {
            const bool hi = w(IDN, k, j, i) > rho_split;
            if (hi == (which == MagDiag::dissAhi)) lsum += q * dV;
          }
          return;
        }

        // --- WP-8 current-sheet split ----------------------------------------------------
        // s = |J| * dx_min / |B| = the fraction of the local field that reverses across one
        // cell. The MINIMUM of the three spacings is used so that on an anisotropic cell the
        // indicator reports the most grid-limited direction rather than an average that could
        // hide it. Unlike the density split, this separates the pathological grid-scale part
        // of Jsq by construction -- that part lives in the diffuse ENVELOPE, which is exactly
        // where a density threshold puts everything in one bin.
        if (which == MagDiag::Jsqsheet || which == MagDiag::Jsqsmth ||
            which == MagDiag::Vsheet) {
          const Real bmag = std::sqrt(bx * bx + by * by + bz * bz);
          // Ternary rather than std::fmin: this lambda is compiled for the CUDA device by
          // nvcc_wrapper, and the CPU build cannot catch a host-only std:: overload. The
          // surrounding code sticks to std::sqrt/std::pow, which nvcc does map; fmin is not
          // worth the risk for a two-way minimum.
          const Real dmin = (dx < dy ? (dx < dz ? dx : dz) : (dy < dz ? dy : dz));
          // |B| == 0 with J != 0 is a pure grid artefact, so it counts as sheet (s = inf).
          const bool sheet =
              (bmag > 0.0) ? (std::sqrt(Jsq) * dmin / bmag > sheet_thresh) : (Jsq > 0.0);
          if (which == MagDiag::Vsheet) {
            if (sheet) lsum += dV;
          } else if (sheet == (which == MagDiag::Jsqsheet)) {
            lsum += Jsq * dV;
          }
          return;
        }
        if (which == MagDiag::Jsqsq) {
          lsum += Jsq * Jsq * dV; // int |J|^4 dV -- gives f_eff(Jsq); see the header
          return;
        }

        if (which == MagDiag::Jsq) {
          lsum += Jsq * dV;
        } else if (which == MagDiag::Hc) {
          lsum += (bx * Jx + by * Jy + bz * Jz) * dV;
        } else if (which == MagDiag::dissO || which == MagDiag::dissOcap) {
          if (which == MagDiag::dissOcap &&
              capd(b, CapDiagIdx::flagO, k, j, i) == 0.0) {
            return;
          }
          lsum += eta(b, NonidealEtaIdx::O, k, j, i) * Jsq * dV;
        } else if (which == MagDiag::dissA || which == MagDiag::dissAcap) {
          if (which == MagDiag::dissAcap &&
              capd(b, CapDiagIdx::flagA, k, j, i) == 0.0) {
            return;
          }
          // Ambipolar heats on the PERPENDICULAR current only: J_perp = J - (J.bhat)bhat.
          const Real b2 = bx * bx + by * by + bz * bz;
          Real Jperp2 = Jsq;
          if (b2 > 0.0) {
            const Real Jpar = (Jx * bx + Jy * by + Jz * bz) / std::sqrt(b2);
            Jperp2 = Jsq - Jpar * Jpar;
            if (Jperp2 < 0.0) Jperp2 = 0.0; // round-off guard
          }
          lsum += eta(b, NonidealEtaIdx::A, k, j, i) * Jperp2 * dV;
        }
      },
      sum);
  return sum;
}

} // namespace Diagnostics
