//========================================================================================
// AthenaPK - radiation transport validation problem (M1, increment 2b).
// A Gaussian radiation-energy pulse on a static, optically-thin uniform gas, with the
// radiation flux set to a beam of reduced flux f (default f=1 => free-streaming).
// For f=1 the M1 signal speed is exactly chat, so the pulse must advect at chat
// without changing shape. Set reduced_flux=0 instead to test the diffusion-limit
// spreading speed chat/sqrt(3). Requires <physics> radiation = true.
//========================================================================================

#include <cmath>
#include <string>

#include <parthenon/package.hpp>

#include "../main.hpp"
#include "pgen.hpp"

namespace rad_pulse {
using namespace parthenon::driver::prelude;

void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
  auto hydro_pkg = pmb->packages.Get("Hydro");
  const bool mhd = (hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd);
  const Real gam = pin->GetReal("hydro", "gamma");
  const Real gm1 = gam - 1.0;

  // Static uniform background gas
  const Real rho0 = pin->GetOrAddReal("problem/rad_pulse", "rho0", 1.0);
  const Real p0 = pin->GetOrAddReal("problem/rad_pulse", "p0", 1.0);

  // Radiation pulse: Er = E0 + Eamp * exp(-(x-x0)^2 / 2 width^2); beam reduced flux fred
  const Real E0 = pin->GetOrAddReal("problem/rad_pulse", "E0", 1.0e-10);
  const Real Eamp = pin->GetOrAddReal("problem/rad_pulse", "Eamp", 1.0);
  const Real x0 = pin->GetOrAddReal("problem/rad_pulse", "x0", 0.0);
  const Real width = pin->GetOrAddReal("problem/rad_pulse", "width", 0.05);
  const Real fred = pin->GetOrAddReal("problem/rad_pulse", "reduced_flux", 1.0);

  auto rad_pkg = pmb->packages.Get("radiation");
  const Real chat = rad_pkg->Param<Real>("chat");

  auto &data = pmb->meshblock_data.Get();
  auto u = data->Get("cons").data.GetHostMirrorAndCopy();
  auto Er = data->Get("rad.Er").data.GetHostMirrorAndCopy();
  auto Fx = data->Get("rad.Fr1").data.GetHostMirrorAndCopy();
  auto Fy = data->Get("rad.Fr2").data.GetHostMirrorAndCopy();
  auto Fz = data->Get("rad.Fr3").data.GetHostMirrorAndCopy();

  parthenon::IndexRange ib = pmb->cellbounds.GetBoundsI(parthenon::IndexDomain::entire);
  parthenon::IndexRange jb = pmb->cellbounds.GetBoundsJ(parthenon::IndexDomain::entire);
  parthenon::IndexRange kb = pmb->cellbounds.GetBoundsK(parthenon::IndexDomain::entire);
  auto &coords = pmb->coords;

  for (int k = kb.s; k <= kb.e; ++k) {
    for (int j = jb.s; j <= jb.e; ++j) {
      for (int i = ib.s; i <= ib.e; ++i) {
        const Real x1 = coords.Xc<1>(i);
        // static uniform gas, v = 0
        u(IDN, k, j, i) = rho0;
        u(IM1, k, j, i) = 0.0;
        u(IM2, k, j, i) = 0.0;
        u(IM3, k, j, i) = 0.0;
        u(IEN, k, j, i) = p0 / gm1;
        if (mhd) {
          u(IB1, k, j, i) = 0.0;
          u(IB2, k, j, i) = 0.0;
          u(IB3, k, j, i) = 0.0;
          u(IPS, k, j, i) = 0.0;
        }
        // radiation pulse + beam
        const Real arg = (x1 - x0) / width;
        const Real E = E0 + Eamp * std::exp(-0.5 * arg * arg);
        Er(0, k, j, i) = E;
        Fx(0, k, j, i) = fred * chat * E; // reduced flux fred in +x1
        Fy(0, k, j, i) = 0.0;
        Fz(0, k, j, i) = 0.0;
      }
    }
  }
  data->Get("cons").data.DeepCopy(u);
  data->Get("rad.Er").data.DeepCopy(Er);
  data->Get("rad.Fr1").data.DeepCopy(Fx);
  data->Get("rad.Fr2").data.DeepCopy(Fy);
  data->Get("rad.Fr3").data.DeepCopy(Fz);
}

} // namespace rad_pulse
