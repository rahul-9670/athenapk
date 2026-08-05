//========================================================================================
// AthenaPK - Poisson / self-gravity boundary-condition test problem
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
//! \file poisson_test.cpp
//  \brief Static uniform-density sphere for validating the self-gravity Poisson solve
//         and its boundary conditions (WS-5a multipole BC gate (i)). The gas is set at
//         rest with a uniform pressure; only the FIRST Poisson solve (from the analytic
//         IC density) is of interest. Compare grav.phi to the analytic uniform-sphere
//         potential:  phi_in(r) = -(2/3)pi G rho (3R^2 - r^2),  phi_out(r) = -G M / r,
//         with G = four_pi_G/(4 pi), M = (4/3) pi R^3 rho.

#include <cmath>

#include <parthenon/package.hpp>

#include "../main.hpp"
#include "pgen.hpp"

namespace poisson_test {
using namespace parthenon::driver::prelude;

void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
  const Real gam = pin->GetReal("hydro", "gamma");
  const Real gm1 = gam - 1.0;

  const Real R       = pin->GetOrAddReal("problem/poisson_test", "radius", 4.0);
  const Real rho_in  = pin->GetOrAddReal("problem/poisson_test", "rho_in", 1.0);
  const Real rho_out = pin->GetOrAddReal("problem/poisson_test", "rho_out", 1.0e-8);
  const Real pres    = pin->GetOrAddReal("problem/poisson_test", "pressure", 1.0);
  const Real cx      = pin->GetOrAddReal("problem/poisson_test", "cx", 0.0);
  const Real cy      = pin->GetOrAddReal("problem/poisson_test", "cy", 0.0);
  const Real cz      = pin->GetOrAddReal("problem/poisson_test", "cz", 0.0);
  // Optional uniform bulk velocity so a density bump advects (used to force dynamic AMR
  // regridding in the sink swarm test); default 0 preserves the gas-at-rest gate (i) setup.
  const Real bvx     = pin->GetOrAddReal("problem/poisson_test", "bulk_vx", 0.0);
  const Real bvy     = pin->GetOrAddReal("problem/poisson_test", "bulk_vy", 0.0);
  const Real bvz     = pin->GetOrAddReal("problem/poisson_test", "bulk_vz", 0.0);
  // Optional Gaussian density peak + radial inflow (converging flow) for the sink-creation
  // test: rho = rho_out + (rho_in-rho_out) exp(-r^2/2 sigma^2), v = -v_radial * r_hat.
  const bool gaussian = pin->GetOrAddString("problem/poisson_test", "profile",
                                            "tophat") == "gaussian";
  const Real sigma   = pin->GetOrAddReal("problem/poisson_test", "sigma", 1.0);
  const Real v_radial = pin->GetOrAddReal("problem/poisson_test", "v_radial", 0.0);
  // AUDIT 2026-08-05: optional uniform B along z, so this test problem can be run under
  // fluid=glmmhd with a real magnetic energy density (Heaviside-Lorentz, ME = B^2/2). Needed
  // by the A4 sink-accretion regression, which has to reach the ME > e_th regime. Default 0
  // reproduces the previous behaviour exactly, so every existing poisson_test deck is
  // bit-identical.
  const Real b0z = pin->GetOrAddReal("problem/poisson_test", "b0z", 0.0);

  auto hydro_pkg = pmb->packages.Get("Hydro");
  const bool mhd = (hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd);

  auto &data = pmb->meshblock_data.Get();
  auto &u_dev = data->Get("cons").data;
  auto u = u_dev.GetHostMirrorAndCopy();
  auto &coords = pmb->coords;

  parthenon::IndexRange ib = pmb->cellbounds.GetBoundsI(parthenon::IndexDomain::entire);
  parthenon::IndexRange jb = pmb->cellbounds.GetBoundsJ(parthenon::IndexDomain::entire);
  parthenon::IndexRange kb = pmb->cellbounds.GetBoundsK(parthenon::IndexDomain::entire);
  const Real R2 = R * R;

  for (int k = kb.s; k <= kb.e; ++k) {
    const Real z = coords.Xc<3>(k) - cz;
    for (int j = jb.s; j <= jb.e; ++j) {
      const Real y = coords.Xc<2>(j) - cy;
      for (int i = ib.s; i <= ib.e; ++i) {
        const Real x = coords.Xc<1>(i) - cx;
        const Real r2 = x * x + y * y + z * z;
        Real rho, vx, vy, vz;
        if (gaussian) {
          rho = rho_out + (rho_in - rho_out) * std::exp(-r2 / (2.0 * sigma * sigma));
          const Real rr = std::sqrt(r2);
          const Real vr = (rr > 1e-12) ? -v_radial / rr : 0.0; // radial inflow -> div v < 0
          vx = vr * x + bvx; vy = vr * y + bvy; vz = vr * z + bvz;
        } else {
          rho = (r2 < R2) ? rho_in : rho_out;
          vx = bvx; vy = bvy; vz = bvz;
        }
        u(IDN, k, j, i) = rho;
        u(IM1, k, j, i) = rho * vx;
        u(IM2, k, j, i) = rho * vy;
        u(IM3, k, j, i) = rho * vz;
        u(IEN, k, j, i) = pres / gm1 + 0.5 * rho * (vx * vx + vy * vy + vz * vz);
        if (mhd) {
          u(IB1, k, j, i) = 0.0;
          u(IB2, k, j, i) = 0.0;
          u(IB3, k, j, i) = b0z;
          u(IPS, k, j, i) = 0.0;
          // total energy carries the magnetic term (0 when b0z = 0 => bit-identical)
          u(IEN, k, j, i) += 0.5 * b0z * b0z;
        }
      }
    }
  }
  u_dev.DeepCopy(u);
}

} // namespace poisson_test
