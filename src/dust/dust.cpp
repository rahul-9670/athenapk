//========================================================================================
// AthenaPK WS-4 dust package — implementation. Two passive scalars carry the dust-to-gas
// ratio f_dg and the characteristic grain radius a_c [cm]; this package adds the operator-
// split growth + sublimation source (DustModel::integrate_cell). Gated behind <physics>
// dust=true (default false) => bit-identical to pre-dust runs.
//========================================================================================
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

#include "dust_pkg.hpp"
#include "dust.hpp" // DustModel (after parthenon => DUST_FN = device)
// NSPEC_H2 / NSPEC_GOW for the F3 scalar-index collision guard (after parthenon headers
// so the networks' CHEMG_FN/CHEM_FN resolve to the device-callable form).
#include "../chemistry/network_h2.hpp"
#include "../chemistry/network_gow17_reduced.hpp"
#include "../main.hpp"

using namespace parthenon::package::prelude;

namespace Dust {

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {
  auto pkg = std::make_shared<StateDescriptor>("dust");
  const std::string bn = "dust";

  const bool evolve = pin->GetOrAddBoolean(bn, "evolve", false);
  pkg->AddParam("evolve", evolve);

  DustModel d;
  d.rho_unit = pin->GetOrAddReal(bn, "rho_unit_cgs", 5.467e-19);
  d.t_unit = pin->GetOrAddReal(bn, "t_unit_cgs", 1.476998822551e12);
  d.rho_grain = pin->GetOrAddReal(bn, "rho_grain_cgs", 3.0);
  d.a_ref = pin->GetOrAddReal(bn, "a_ref_cm", 1.0e-5);
  d.f_dg_ref = pin->GetOrAddReal(bn, "f_dg_ref", 0.01);
  d.alpha_turb = pin->GetOrAddReal(bn, "alpha_turb", 1.0e-2);
  d.omega_ref = pin->GetOrAddReal(bn, "omega_ref_cgs", 1.0e-13);
  d.mu_gas = pin->GetOrAddReal(bn, "mu_n", 2.33);
  d.T_subl = pin->GetOrAddReal(bn, "T_subl_K", 1.5e3);
  d.growth_cap = pin->GetOrAddReal(bn, "growth_cap", 0.10);
  // zero-growth switch for the inc-2 bit-identical gate (evolve on, but grains frozen).
  const bool freeze_growth = pin->GetOrAddBoolean(bn, "freeze_growth", false);
  pkg->AddParam("freeze_growth", freeze_growth);
  pkg->AddParam("model", d);
  pkg->AddParam("T_unit_K", pin->GetOrAddReal(bn, "T_unit_K", 10.015));
  pkg->AddParam("gamma", pin->GetReal("hydro", "gamma"));
  // Scalar block index of f_dg (a_c = f_dg index + 1). Placed after the chemistry species.
  pkg->AddParam("scalar_index", pin->GetOrAddInteger(bn, "scalar_index", 0));

  if (evolve) {
    const int nscalars = pin->GetOrAddInteger("hydro", "nscalars", 0);
    const int si = pin->GetOrAddInteger(bn, "scalar_index", 0);
    PARTHENON_REQUIRE(nscalars >= si + 2,
                      "dust: <hydro> nscalars must be >= dust/scalar_index + 2 (f_dg, a_c).");
    // F3 fix (audit 2026-07-18): when chemistry is active its species occupy scalar
    // indices [0, NSPEC); the default scalar_index=0 would alias f_dg/a_c onto x_H2/x_Hp
    // with no error. Require the dust pair to sit past the chemistry block.
    if (pin->GetOrAddBoolean("physics", "chemistry", false)) {
      const std::string net = pin->GetOrAddString("chemistry", "network", "H2");
      const int nspec_chem =
          (net == "gow17_reduced") ? Chemistry::NSPEC_GOW : Chemistry::NSPEC_H2;
      PARTHENON_REQUIRE(si >= nspec_chem,
                        "dust: dust/scalar_index (" + std::to_string(si) +
                            ") collides with the chemistry species block; it must be >= " +
                            std::to_string(nspec_chem) + " for network '" + net + "'.");
    }
  }
  return pkg;
}

TaskStatus ReactDust(MeshData<Real> *md, const Real dt) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  auto pkg = pmb->pmy_mesh->packages.Get("dust");
  if (!pkg->Param<bool>("evolve")) return TaskStatus::complete;

  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  const int nhydro = pmb->pmy_mesh->packages.Get("Hydro")->Param<int>("nhydro");
  const bool mhd = (nhydro > NHYDRO);
  const DustModel d = pkg->Param<DustModel>("model");
  const bool freeze = pkg->Param<bool>("freeze_growth");
  const Real gm1 = pkg->Param<Real>("gamma") - 1.0;
  const Real T_unit = pkg->Param<Real>("T_unit_K");
  const int sidx = pkg->Param<int>("scalar_index");

  static const std::vector<std::string> names{"cons"};
  parthenon::PackIndexMap imap;
  auto pack = md->PackVariables(names, imap);
  const int ic = imap["cons"].first;
  const int is = ic + nhydro + sidx; // f_dg scalar_density; a_c at is+1

  const Real rho_to_cgs = d.rho_unit;
  const Real dt_s = static_cast<Real>(dt) * d.t_unit;

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Dust::React", parthenon::DevExecSpace(), 0, pack.GetDim(5) - 1,
      kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const Real rho = pack(b, ic + IDN, k, j, i);
        const Real inv_rho = 1.0 / rho;
        const Real m1 = pack(b, ic + IM1, k, j, i);
        const Real m2 = pack(b, ic + IM2, k, j, i);
        const Real m3 = pack(b, ic + IM3, k, j, i);
        Real eint = pack(b, ic + IEN, k, j, i) - 0.5 * (m1 * m1 + m2 * m2 + m3 * m3) * inv_rho;
        if (mhd) {
          const Real b1 = pack(b, ic + IB1, k, j, i);
          const Real b2 = pack(b, ic + IB2, k, j, i);
          const Real b3 = pack(b, ic + IB3, k, j, i);
          eint -= 0.5 * (b1 * b1 + b2 * b2 + b3 * b3);
        }
        Real T_code = gm1 * eint * inv_rho;
        if (T_code < 1.0e-4) T_code = 1.0e-4;
        const double T_K = static_cast<double>(T_code * T_unit);
        const double rho_cgs = static_cast<double>(rho * rho_to_cgs);

        double f_dg = static_cast<double>(pack(b, is, k, j, i) * inv_rho);
        double a_c = static_cast<double>(pack(b, is + 1, k, j, i) * inv_rho);
        if (freeze) {
          // bit-identical guard: only apply sublimation switch, no growth.
          f_dg = d.f_dg_ref * d.subl_factor(T_K);
        } else {
          d.integrate_cell(a_c, f_dg, rho_cgs, T_K, static_cast<double>(dt_s));
        }
        pack(b, is, k, j, i) = rho * static_cast<Real>(f_dg);
        pack(b, is + 1, k, j, i) = rho * static_cast<Real>(a_c);
      });
  return TaskStatus::complete;
}

void AddDustTasks(TaskCollection &tc, Mesh *pmesh, const Real dt) {
  using namespace parthenon;
  if (!pmesh->packages.Get("dust")->Param<bool>("evolve")) return;
  auto partitions = pmesh->GetDefaultBlockPartitions();
  const int num_partitions = partitions.size();
  TaskRegion &region = tc.AddRegion(num_partitions);
  for (int i = 0; i < num_partitions; ++i) {
    TaskList &tl = region[i];
    auto &md = pmesh->mesh_data.Add("base", partitions[i]);
    tl.AddTask(TaskID(0), ReactDust, md.get(), dt);
  }
}

} // namespace Dust
