//========================================================================================
// AthenaPK - M1 radiation transport validation: 1D radiative shock (colliding flows).
// Two equal, opposed supersonic streams (v = +v0 for x<0, -v0 for x>0) collide at x=0,
// driving a symmetric pair of shocks that propagate outward. The hot post-shock slab
// radiates; with matter coupling + finite absorption opacity a RADIATIVE PRECURSOR forms
// (upstream gas is preheated by downstream radiation AHEAD of the density jump), and behind
// each shock the gas and radiation temperatures relax toward equilibrium (T_gas -> T_rad,
// Er -> arad*T^4). Both x-boundaries are outflow; the symmetry keeps the shock structure at
// the center. Requires <physics> radiation = true. Set <radiation> arad_code=1 to put the
// radiation and gas energies on comparable footing (radiation-coupled regime).
//========================================================================================

#include <cmath>
#include <string>

#include <parthenon/package.hpp>

#include "../main.hpp"
#include "pgen.hpp"

namespace rad_shock {
using namespace parthenon::driver::prelude;

void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
  const Real gam = pin->GetReal("hydro", "gamma");
  const Real gm1 = gam - 1.0;

  const Real rho0 = pin->GetOrAddReal("problem/rad_shock", "rho0", 1.0);
  const Real T0 = pin->GetOrAddReal("problem/rad_shock", "T0", 0.2);
  const Real v0 = pin->GetOrAddReal("problem/rad_shock", "v0", 2.5);
  const Real smooth = pin->GetOrAddReal("problem/rad_shock", "smooth", 0.008);
  const Real p0 = rho0 * T0; // ideal: T = p/rho (code units)

  auto rad_pkg = pmb->packages.Get("radiation");
  const Real arad = rad_pkg->Param<Real>("arad");
  const Real Er0 = arad * T0 * T0 * T0 * T0; // equilibrium with the cold upstream gas

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
        // v = -v0*tanh(x/smooth): +v0 far left, -v0 far right (colliding at x=0)
        const Real v = -v0 * std::tanh(x1 / smooth);
        u(IDN, k, j, i) = rho0;
        u(IM1, k, j, i) = rho0 * v;
        u(IM2, k, j, i) = 0.0;
        u(IM3, k, j, i) = 0.0;
        u(IEN, k, j, i) = p0 / gm1 + 0.5 * rho0 * v * v;
        Er(0, k, j, i) = Er0;
        Fx(0, k, j, i) = 0.0;
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

} // namespace rad_shock
