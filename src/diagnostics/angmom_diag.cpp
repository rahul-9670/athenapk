//========================================================================================
// AthenaPK - angular-momentum history diagnostics. See angmom_diag.hpp.
//========================================================================================

#include <cmath>
#include <string>
#include <vector>

#include <parthenon/package.hpp>

#include "../main.hpp"
#include "angmom_diag.hpp"

namespace Diagnostics {

Real AngMomReduce(MeshData<Real> *md, AngMom which) {
  const auto &prim = md->PackVariables(std::vector<std::string>{"prim"});

  // Moment origin = centre of the box (see the header). Taken from mesh_size so it is
  // correct for a non-symmetric box too, and captured by value for the device lambda.
  auto *pm = md->GetBlockData(0)->GetBlockPointer()->pmy_mesh;
  const Real x0 = 0.5 * (pm->mesh_size.xmin(X1DIR) + pm->mesh_size.xmax(X1DIR));
  const Real y0 = 0.5 * (pm->mesh_size.xmin(X2DIR) + pm->mesh_size.xmax(X2DIR));
  const Real z0 = 0.5 * (pm->mesh_size.xmin(X3DIR) + pm->mesh_size.xmax(X3DIR));
  const Real xlo = pm->mesh_size.xmin(X1DIR), xhi = pm->mesh_size.xmax(X1DIR);
  const Real ylo = pm->mesh_size.xmin(X2DIR), yhi = pm->mesh_size.xmax(X2DIR);
  const Real zlo = pm->mesh_size.xmin(X3DIR), zhi = pm->mesh_size.xmax(X3DIR);

  const bool is_flux = (which == AngMom::FLx || which == AngMom::FLy || which == AngMom::FLz);
  const bool is_tflux = (which == AngMom::FTx || which == AngMom::FTy || which == AngMom::FTz);
  // B6: the stage-consistent variant reads the SOLVER's momentum fluxes instead of rebuilding
  // the stress from end-of-step primitives.
  const bool is_sflux = (which == AngMom::FTsolverX || which == AngMom::FTsolverY ||
                         which == AngMom::FTsolverZ);
  // Packed BY NAME, not by the {Independent} metadata flag: grav::phi is also
  // Independent+WithFluxes, so a flag-based pack could put phi's flux at the flat IM1 index.
  // Same trap cons_diag and SelfGravity::ApplyGravitySource both document.
  auto cons_flx = md->PackVariablesAndFluxes(std::vector<std::string>{"cons"},
                                             std::vector<std::string>{"cons"});
  const bool is_torque =
      (which == AngMom::TmagX || which == AngMom::TmagY || which == AngMom::TmagZ);
  const bool is_grav =
      (which == AngMom::TgravX || which == AngMom::TgravY || which == AngMom::TgravZ);

  const bool is_split =
      (which == AngMom::Lxhi || which == AngMom::Lxlo || which == AngMom::Lyhi ||
       which == AngMom::Lylo || which == AngMom::Lzhi || which == AngMom::Lzlo);

  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  const bool mhd = hydro_pkg->Param<bool>("angmom_diag_mhd");
  // Captured by value for the device lambda. 0 => the split columns are never registered
  // (see the header), so this value is only ever read when it is > 0.
  const Real rho_split = hydro_pkg->Param<Real>("angmom_diag_rho_split");
  // Dummy-alias to `prim` when unused: PackVariables on an unregistered field throws and
  // the device lambda must still capture something. Same idiom as mag_diag.cpp.
  const auto &phi =
      is_grav ? md->PackVariables(std::vector<std::string>{"grav.phi"}) : prim;

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real sum = 0.0;
  Kokkos::parallel_reduce(
      "AngMomReduce",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExecSpace(), {0, kb.s, jb.s, ib.s},
          {prim.GetDim(5), kb.e + 1, jb.e + 1, ib.e + 1}, {1, 1, 1, ib.e + 1 - ib.s}),
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &lsum) {
        const auto &w = prim(b);
        const auto &coords = prim.GetCoords(b);
        const Real dV = coords.CellVolume(k, j, i);
        const Real dx = coords.Dxc<1>(k, j, i);
        const Real dy = coords.Dxc<2>(k, j, i);
        const Real dz = coords.Dxc<3>(k, j, i);

        // Position relative to the moment origin.
        const Real x = coords.Xc<1>(i) - x0;
        const Real y = coords.Xc<2>(j) - y0;
        const Real z = coords.Xc<3>(k) - z0;

        if (is_torque) {
          // J = curl B, 2*dx cell-centered central difference -- the same stencil as
          // mag_diag and the jeans_nonideal criterion (prim is FillGhost, ghosts valid).
          const Real Jx = (w(IB3, k, j + 1, i) - w(IB3, k, j - 1, i)) / (2.0 * dy) -
                          (w(IB2, k + 1, j, i) - w(IB2, k - 1, j, i)) / (2.0 * dz);
          const Real Jy = (w(IB1, k + 1, j, i) - w(IB1, k - 1, j, i)) / (2.0 * dz) -
                          (w(IB3, k, j, i + 1) - w(IB3, k, j, i - 1)) / (2.0 * dx);
          const Real Jz = (w(IB2, k, j, i + 1) - w(IB2, k, j, i - 1)) / (2.0 * dx) -
                          (w(IB1, k, j + 1, i) - w(IB1, k, j - 1, i)) / (2.0 * dy);
          const Real bx = w(IB1, k, j, i), by = w(IB2, k, j, i), bz = w(IB3, k, j, i);
          // Lorentz force density, Heaviside-Lorentz: F = J x B (no 4pi/c).
          const Real Fx = Jy * bz - Jz * by;
          const Real Fy = Jz * bx - Jx * bz;
          const Real Fz = Jx * by - Jy * bx;
          Real t = 0.0;
          if (which == AngMom::TmagX) t = y * Fz - z * Fy;
          else if (which == AngMom::TmagY) t = z * Fx - x * Fz;
          else t = x * Fy - y * Fx;
          lsum += t * dV;
          return;
        }

        const Real rho = w(IDN, k, j, i);

        if (is_grav) {
          // Gravitational torque density rho * (r x (-grad Phi)). grav.phi is FillGhost so
          // the 2*dx central difference is valid on interior cells. In the continuum the
          // TOTAL self-gravitational torque vanishes (central pair forces); what this
          // measures is how far the discrete solve+gradient depart from that.
          const Real gx = -(phi(b, 0, k, j, i + 1) - phi(b, 0, k, j, i - 1)) / (2.0 * dx);
          const Real gy = -(phi(b, 0, k, j + 1, i) - phi(b, 0, k, j - 1, i)) / (2.0 * dy);
          const Real gz = -(phi(b, 0, k + 1, j, i) - phi(b, 0, k - 1, j, i)) / (2.0 * dz);
          Real t = 0.0;
          if (which == AngMom::TgravX) t = y * gz - z * gy;
          else if (which == AngMom::TgravY) t = z * gx - x * gz;
          else t = x * gy - y * gx;
          lsum += rho * t * dV;
          return;
        }

        if (which == AngMom::Mhi) {
          if (rho > rho_split) lsum += rho * dV;
          return;
        }

        const Real vx = w(IV1, k, j, i), vy = w(IV2, k, j, i), vz = w(IV3, k, j, i);
        // Specific angular momentum r x v.
        const Real lx = y * vz - z * vy;
        const Real ly = z * vx - x * vz;
        const Real lz = x * vy - y * vx;

        if (is_split) {
          // Same integrand as Lx/Ly/Lz, restricted to one side of the density split, so
          // (Lxhi + Lxlo) reconstructs Lx. Each cell contributes to exactly one side and the
          // per-cell FP ops are identical to the unsplit branch; the two partial sums are
          // still separate Kokkos reductions, so bit-exactness is not guaranteed under a
          // non-deterministic reduction order. MEASURED (runs/wp4_split_gate/on, 32^3 glmmhd,
          // rho_split=5.0): |hi+lo-L| is at most 0.28x the .hst print-resolution floor on all
          // three axes, i.e. consistent with exact summation -- but the .hst writes 6 sig
          // figs, so it cannot resolve any finer, and exactness is NOT claimed here.
          const bool hi = rho > rho_split;
          Real l;
          if (which == AngMom::Lxhi || which == AngMom::Lxlo) l = lx;
          else if (which == AngMom::Lyhi || which == AngMom::Lylo) l = ly;
          else l = lz;
          const bool want_hi =
              (which == AngMom::Lxhi || which == AngMom::Lyhi || which == AngMom::Lzhi);
          if (hi == want_hi) lsum += rho * l * dV;
          return;
        }

        if (is_sflux) {
          // B6 -- STAGE-CONSISTENT angular-momentum flux. cons.flux(dir, IM1+c) is the momentum
          // flux the Riemann solver actually applied on that face; it already carries the
          // pressure and Maxwell stress, so nothing is reconstructed here and the two cannot
          // drift apart. Sign convention matches FT*: outward-positive, hence the minus on the
          // lower faces (where the outward normal is -e_dir).
          const auto &cfb = cons_flx(b);
          auto face = [&](const int dir, const int kk, const int jj, const int ii,
                          const Real sgn, const Real dA) {
            const Real tx = cfb.flux(dir, IM1, kk, jj, ii);
            const Real ty = cfb.flux(dir, IM2, kk, jj, ii);
            const Real tz = cfb.flux(dir, IM3, kk, jj, ii);
            Real t;
            if (which == AngMom::FTsolverX) t = y * tz - z * ty;
            else if (which == AngMom::FTsolverY) t = z * tx - x * tz;
            else t = x * ty - y * tx;
            return sgn * t * dA;
          };
          Real f = 0.0;
          if (std::fabs(coords.Xf<1>(i) - xlo) < 1.0e-9 * dx)
            f += face(X1DIR, k, j, i, -1.0, dV / dx);
          if (std::fabs(coords.Xf<1>(i + 1) - xhi) < 1.0e-9 * dx)
            f += face(X1DIR, k, j, i + 1, 1.0, dV / dx);
          if (std::fabs(coords.Xf<2>(j) - ylo) < 1.0e-9 * dy)
            f += face(X2DIR, k, j, i, -1.0, dV / dy);
          if (std::fabs(coords.Xf<2>(j + 1) - yhi) < 1.0e-9 * dy)
            f += face(X2DIR, k, j + 1, i, 1.0, dV / dy);
          if (std::fabs(coords.Xf<3>(k) - zlo) < 1.0e-9 * dz)
            f += face(X3DIR, k, j, i, -1.0, dV / dz);
          if (std::fabs(coords.Xf<3>(k + 1) - zhi) < 1.0e-9 * dz)
            f += face(X3DIR, k + 1, j, i, 1.0, dV / dz);
          lsum += f;
          return;
        }

        if (!is_flux && !is_tflux) {
          Real l = (which == AngMom::Lx) ? lx : ((which == AngMom::Ly) ? ly : lz);
          lsum += rho * l * dV;
          return;
        }

        if (is_tflux) {
          // TOTAL angular-momentum flux r x (T.nhat) with the full stress
          //     T_ij = rho v_i v_j + (P + B^2/2) delta_ij - B_i B_j
          // This is the budget term: d/dt L = -surf r x (T.nhat) dA + Tgrav. The pressure
          // part is the one the advective-only column omits, and it dominates under
          // outflow BCs -- see the header.
          const Real pres = w(IPR, k, j, i);
          Real bx = 0.0, by = 0.0, bz = 0.0;
          if (mhd) {
            bx = w(IB1, k, j, i);
            by = w(IB2, k, j, i);
            bz = w(IB3, k, j, i);
          }
          const Real ptot = pres + 0.5 * (bx * bx + by * by + bz * bz);
          auto face = [&](const Real nx, const Real ny, const Real nz, const Real dA) {
            const Real vn = vx * nx + vy * ny + vz * nz;
            const Real bn = bx * nx + by * ny + bz * nz;
            // (T.nhat)_i for i = x,y,z
            const Real tx = rho * vx * vn + ptot * nx - bx * bn;
            const Real ty = rho * vy * vn + ptot * ny - by * bn;
            const Real tz = rho * vz * vn + ptot * nz - bz * bn;
            Real t;
            if (which == AngMom::FTx) t = y * tz - z * ty;
            else if (which == AngMom::FTy) t = z * tx - x * tz;
            else t = x * ty - y * tx;
            return t * dA;
          };
          Real f = 0.0;
          if (std::fabs(coords.Xf<1>(i) - xlo) < 1.0e-9 * dx) f += face(-1, 0, 0, dV / dx);
          if (std::fabs(coords.Xf<1>(i + 1) - xhi) < 1.0e-9 * dx) f += face(1, 0, 0, dV / dx);
          if (std::fabs(coords.Xf<2>(j) - ylo) < 1.0e-9 * dy) f += face(0, -1, 0, dV / dy);
          if (std::fabs(coords.Xf<2>(j + 1) - yhi) < 1.0e-9 * dy) f += face(0, 1, 0, dV / dy);
          if (std::fabs(coords.Xf<3>(k) - zlo) < 1.0e-9 * dz) f += face(0, 0, -1, dV / dz);
          if (std::fabs(coords.Xf<3>(k + 1) - zhi) < 1.0e-9 * dz) f += face(0, 0, 1, dV / dz);
          lsum += f;
          return;
        }

        // --- advective flux through the PHYSICAL domain boundary, outward-positive -------
        // A cell is on a domain face when its own face coordinate coincides with the mesh
        // extent. Corner cells sit on two or three faces and contribute to each.
        const Real l = (which == AngMom::FLx) ? lx : ((which == AngMom::FLy) ? ly : lz);
        const Real rl = rho * l;
        Real f = 0.0;
        if (std::fabs(coords.Xf<1>(i) - xlo) < 1.0e-9 * dx) f += rl * (-vx) * (dV / dx);
        if (std::fabs(coords.Xf<1>(i + 1) - xhi) < 1.0e-9 * dx) f += rl * (vx) * (dV / dx);
        if (std::fabs(coords.Xf<2>(j) - ylo) < 1.0e-9 * dy) f += rl * (-vy) * (dV / dy);
        if (std::fabs(coords.Xf<2>(j + 1) - yhi) < 1.0e-9 * dy) f += rl * (vy) * (dV / dy);
        if (std::fabs(coords.Xf<3>(k) - zlo) < 1.0e-9 * dz) f += rl * (-vz) * (dV / dz);
        if (std::fabs(coords.Xf<3>(k + 1) - zhi) < 1.0e-9 * dz) f += rl * (vz) * (dV / dz);
        lsum += f;
      },
      sum);

  return sum;
}

} // namespace Diagnostics
