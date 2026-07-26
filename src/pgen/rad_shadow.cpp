//========================================================================================
// AthenaPK - M1 radiation transport validation: 2D shadow test.
// A steady radiation beam is injected from the -x boundary (free-streaming, Fr = chat*Er in
// +x) into a static gas containing a dense, optically-thick circular clump. The clump (high
// rho => high absorption optical depth tau = kappa_a*rho*2r >> 1) blocks the beam and must
// cast a SHADOW behind it (Er_behind << Er_beside). This is the canonical M1 test: a single
// source + single obstacle produces a proper shadow (unlike flux-limited diffusion, which
// washes it out). The known M1 failure mode (two crossing beams merging) is NOT exercised
// here. The clump is given a huge heat capacity (high rho) so absorption re-emission stays
// negligible and the shadow stays dark. Requires <physics> radiation = true.
//========================================================================================

#include <algorithm>
#include <cmath>
#include <string>

#include <parthenon/package.hpp>

#include "../main.hpp"
#include "pgen.hpp"

namespace rad_shadow {
using namespace parthenon::driver::prelude;

// Beam state (code units), set in InitUserMeshData, read by the inflow BC.
Real E_beam = 1.0;
Real rho_bg = 1.0;
Real eint_bg = 1.0; // p_bg/(gamma-1), precomputed so the BC needs no package lookup

void InitUserMeshData(Mesh *mesh, parthenon::ParameterInput *pin) {
  E_beam = pin->GetOrAddReal("problem/rad_shadow", "E_beam", 1.0);
  rho_bg = pin->GetOrAddReal("problem/rad_shadow", "rho_bg", 1.0);
  const Real p_bg = pin->GetOrAddReal("problem/rad_shadow", "p_bg", 1.0);
  const Real gm1 = pin->GetReal("hydro", "gamma") - 1.0;
  eint_bg = p_bg / gm1;
}

void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
  const Real gam = pin->GetReal("hydro", "gamma");
  const Real gm1 = gam - 1.0;

  const Real rho0 = pin->GetOrAddReal("problem/rad_shadow", "rho_bg", 1.0);
  const Real p0 = pin->GetOrAddReal("problem/rad_shadow", "p_bg", 1.0);
  const Real rho_cl = pin->GetOrAddReal("problem/rad_shadow", "rho_clump", 1.0e3);
  const Real xc = pin->GetOrAddReal("problem/rad_shadow", "xc", -0.1);
  const Real yc = pin->GetOrAddReal("problem/rad_shadow", "yc", 0.0);
  const Real rcl = pin->GetOrAddReal("problem/rad_shadow", "r_clump", 0.1);
  const Real steep = pin->GetOrAddReal("problem/rad_shadow", "steepness", 40.0);
  const Real E0 = pin->GetOrAddReal("problem/rad_shadow", "E_floor_init", 1.0e-6);

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
        const Real x2 = coords.Xc<2>(j);
        // Smooth dense clump: rho = rho0 + (rho_cl-rho0) * 1/(1+exp(steep*(r-rcl))).
        const Real rr = std::sqrt((x1 - xc) * (x1 - xc) + (x2 - yc) * (x2 - yc));
        const Real w = 1.0 / (1.0 + std::exp(steep * (rr - rcl)));
        const Real rho = rho0 + (rho_cl - rho0) * w;
        u(IDN, k, j, i) = rho;
        u(IM1, k, j, i) = 0.0;
        u(IM2, k, j, i) = 0.0;
        u(IM3, k, j, i) = 0.0;
        u(IEN, k, j, i) = p0 / gm1; // static, uniform pressure
        // Dark initial radiation field; the beam enters through the -x boundary.
        Er(0, k, j, i) = E0;
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

//! Steady free-streaming radiation beam injected at the inner-x1 face. Fills the hydro
//! ghosts with the static background and the radiation ghosts with Er=E_beam,
//! Fr1=chat*E_beam (reduced flux f=1 => pure +x beam).
void InflowBeamX1(std::shared_ptr<MeshBlockData<Real>> &mbd, bool coarse) {
  auto pmb = mbd->GetBlockPointer();
  const auto nb = IndexRange{0, 0};
  const bool fine = false;

  // The hydro ("cons") and the operator-split radiation ("rad.*") variables are exchanged on
  // SEPARATE containers, so this BC is invoked once per container. Fill only what is present.
  if (mbd->Contains(std::string("cons"))) {
    auto cons = mbd->PackVariables(std::vector<std::string>{"cons"}, coarse);
    const Real rhob = rho_bg;
    const Real eintb = eint_bg;
    pmb->par_for_bndry(
        "rad_shadow::InflowBeamX1_hydro", nb, parthenon::IndexDomain::inner_x1,
        parthenon::TopologicalElement::CC, coarse, fine,
        KOKKOS_LAMBDA(const int &, const int &k, const int &j, const int &i) {
          cons(IDN, k, j, i) = rhob;
          cons(IM1, k, j, i) = 0.0;
          cons(IM2, k, j, i) = 0.0;
          cons(IM3, k, j, i) = 0.0;
          cons(IEN, k, j, i) = eintb;
        });
  }
  if (mbd->Contains(std::string("rad.Er"))) {
    const Real chat = pmb->packages.Get("radiation")->Param<Real>("chat");
    auto rad = mbd->PackVariables(
        std::vector<std::string>{"rad.Er", "rad.Fr1", "rad.Fr2", "rad.Fr3"}, coarse);
    const Real Eb = E_beam;
    pmb->par_for_bndry(
        "rad_shadow::InflowBeamX1_rad", nb, parthenon::IndexDomain::inner_x1,
        parthenon::TopologicalElement::CC, coarse, fine,
        KOKKOS_LAMBDA(const int &, const int &k, const int &j, const int &i) {
          rad(0, k, j, i) = Eb;        // rad.Er
          rad(1, k, j, i) = chat * Eb; // rad.Fr1  (f=1 beam in +x)
          rad(2, k, j, i) = 0.0;       // rad.Fr2
          rad(3, k, j, i) = 0.0;       // rad.Fr3
        });
  }
}

} // namespace rad_shadow
