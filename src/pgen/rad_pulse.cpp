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
#include <vector>

#include <parthenon/package.hpp>

#include "../main.hpp"
#include "../radiation/radiation_groups.hpp" // GroupFieldNames (multigroup seeding)
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
  // Multigroup seeding: split the pulse EQUALLY across the n_group groups (each group gets
  // 1/n_group of E and F, i.e. the SAME reduced flux). Then sum_g E_g,0 = E (the gray IC),
  // and because every group has identical reduced flux the M1 closure is identical per group,
  // so under free-streaming each group evolves as (1/n_group)*gray and the group SUM tracks
  // the gray run cell-for-cell -- the multigroup==gray equivalence gate. n_group=1 => share=1
  // => group 0 carries the full pulse (bit-identical to the pre-multigroup pgen).
  const int n_group = rad_pkg->Param<int>("n_group");
  const Real share = 1.0 / static_cast<Real>(n_group);

  auto &data = pmb->meshblock_data.Get();
  auto u = data->Get("cons").data.GetHostMirrorAndCopy();
  // Per-group host mirrors of (Er, Fr1, Fr2, Fr3).
  std::vector<decltype(data->Get("rad.Er").data.GetHostMirrorAndCopy())> gEr, gFx, gFy, gFz;
  for (int g = 0; g < n_group; ++g) {
    const auto nm = Radiation::GroupFieldNames(g);
    gEr.push_back(data->Get(nm[0]).data.GetHostMirrorAndCopy());
    gFx.push_back(data->Get(nm[1]).data.GetHostMirrorAndCopy());
    gFy.push_back(data->Get(nm[2]).data.GetHostMirrorAndCopy());
    gFz.push_back(data->Get(nm[3]).data.GetHostMirrorAndCopy());
  }

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
        // radiation pulse + beam, split equally across groups
        const Real arg = (x1 - x0) / width;
        const Real E = share * (E0 + Eamp * std::exp(-0.5 * arg * arg));
        for (int g = 0; g < n_group; ++g) {
          gEr[g](0, k, j, i) = E;
          gFx[g](0, k, j, i) = fred * chat * E; // reduced flux fred in +x1
          gFy[g](0, k, j, i) = 0.0;
          gFz[g](0, k, j, i) = 0.0;
        }
      }
    }
  }
  data->Get("cons").data.DeepCopy(u);
  for (int g = 0; g < n_group; ++g) {
    const auto nm = Radiation::GroupFieldNames(g);
    data->Get(nm[0]).data.DeepCopy(gEr[g]);
    data->Get(nm[1]).data.DeepCopy(gFx[g]);
    data->Get(nm[2]).data.DeepCopy(gFy[g]);
    data->Get(nm[3]).data.DeepCopy(gFz[g]);
  }
}

} // namespace rad_pulse
