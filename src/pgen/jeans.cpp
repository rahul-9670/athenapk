//========================================================================================
// AthenaPK - Jeans dispersion problem generator for self-gravity validation
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
//! \file jeans.cpp
//  \brief Uniform periodic box with a sinusoidal density perturbation of wavelength
//         equal to the box size in x1, for validating the self-gravity module
//         against the Jeans dispersion relation:
//             omega^2 = k^2 c_s^2 - 4 pi G rho0

#include <cmath>
#include <string>

#include <parthenon/package.hpp>

#include "../main.hpp"
#include "pgen.hpp"

namespace jeans {
using namespace parthenon::driver::prelude;

void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
  auto hydro_pkg = pmb->packages.Get("Hydro");
  const bool mhd = (hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd);
  const Real gam = pin->GetReal("hydro", "gamma");
  const Real gm1 = gam - 1.0;

  // Problem parameters
  const Real rho0 = pin->GetOrAddReal("problem/jeans", "rho0", 1.0);
  const Real cs   = pin->GetOrAddReal("problem/jeans", "cs",   1.0);
  const Real amp  = pin->GetOrAddReal("problem/jeans", "amp",  1.0e-4);
  // Wavenumber in units of 2*pi / Lx (so n_wave=1 means wavelength = Lx)
  const int  nwave = pin->GetOrAddInteger("problem/jeans", "nwave", 1);
  const Real B0z = mhd ? pin->GetOrAddReal("problem/jeans", "B0z", 0.0) : 0.0;

  // Uniform pressure consistent with given rho0 and cs (adiabatic):
  //   c_s^2 = gamma * p / rho  ->  p = rho * c_s^2 / gamma
  const Real p0 = rho0 * cs * cs / gam;

  // Box length in x1 (from mesh) — used to set k
  const Real x1min = pin->GetReal("parthenon/mesh", "x1min");
  const Real x1max = pin->GetReal("parthenon/mesh", "x1max");
  const Real Lx = x1max - x1min;
  const Real k_wave = 2.0 * M_PI * static_cast<Real>(nwave) / Lx;

  auto &data = pmb->meshblock_data.Get();
  auto &u_dev = data->Get("cons").data;
  // Mirror to host for init, then copy back
  auto u = u_dev.GetHostMirrorAndCopy();

  parthenon::IndexRange ib = pmb->cellbounds.GetBoundsI(parthenon::IndexDomain::entire);
  parthenon::IndexRange jb = pmb->cellbounds.GetBoundsJ(parthenon::IndexDomain::entire);
  parthenon::IndexRange kb = pmb->cellbounds.GetBoundsK(parthenon::IndexDomain::entire);

  auto &coords = pmb->coords;
  for (int k = kb.s; k <= kb.e; ++k) {
    for (int j = jb.s; j <= jb.e; ++j) {
      for (int i = ib.s; i <= ib.e; ++i) {
        const Real x1 = coords.Xc<1>(i);
        const Real rho = rho0 * (1.0 + amp * std::sin(k_wave * x1));
        const Real v1 = 0.0;
        const Real v2 = 0.0;
        const Real v3 = 0.0;
        const Real p = p0 * (1.0 + amp * std::sin(k_wave * x1));  // isothermal-like init

        u(IDN, k, j, i) = rho;
        u(IM1, k, j, i) = rho * v1;
        u(IM2, k, j, i) = rho * v2;
        u(IM3, k, j, i) = rho * v3;
        u(IEN, k, j, i) = p / gm1 + 0.5 * rho * (v1*v1 + v2*v2 + v3*v3);
        if (mhd) {
          u(IB1, k, j, i) = 0.0;
          u(IB2, k, j, i) = 0.0;
          u(IB3, k, j, i) = B0z;
          u(IPS, k, j, i) = 0.0;
          u(IEN, k, j, i) += 0.5 * B0z * B0z;  // add B^2/2 to total energy
        }
      }
    }
  }
  u_dev.DeepCopy(u);
}

} // namespace jeans
