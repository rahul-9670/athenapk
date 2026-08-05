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
#include "../eos/adiabatic_glmmhd.hpp" // EOS-table T under eos=hydrogen (audit N2)
#include "../main.hpp"
#include "../units/physical_units.hpp" // the one authoritative unit system (audit N3)

using namespace parthenon::package::prelude;

namespace Dust {

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {
  auto pkg = std::make_shared<StateDescriptor>("dust");
  const std::string bn = "dust";

  const bool evolve = pin->GetOrAddBoolean(bn, "evolve", false);
  pkg->AddParam("evolve", evolve);

  // AUDIT 2026-08-05 (N3) -- ONE UNIT SYSTEM, NOT A PRIVATE COPY.
  // These used to be read from <dust> with hardcoded defaults (rho 5.467e-19, t
  // 1.476998822551e12, T 10.015) -- the rounded BE values for ONE initial condition
  // (mass=6, temperature=10, f=5). Every other package (chemistry, radiation, the
  // diffusion ionization model) was routed through PhysUnits::BuildPhysicalUnits in the
  // flagship Phase-1 consolidation; dust was missed. That mattered structurally rather
  // than numerically: BuildPhysicalUnits deliberately HARD-FAILS on a <units> base-scale
  // override for collapse_be so the microphysics cannot desync from the IC, and <dust>'s
  // private knobs walked straight past that guard. Drift for the current deck is small
  // (rho +0.0031%, T -0.0014%, t -0.0000%) only because the defaults happen to be this
  // IC's rounded values; change mass/temperature/f -- a mu-ladder, a resolution ladder --
  // and dust would silently keep the 6/10/5 numbers while the dynamics moved.
  // The REQUIRE_THROWS below breaks no existing deck: a scan of all 295 *.in files on
  // 2026-08-05 found ZERO that set any of these three keys inside <dust>.
  const auto U = PhysUnits::BuildPhysicalUnits(pin);
  PARTHENON_REQUIRE_THROWS(
      !(pin->DoesParameterExist(bn, "rho_unit_cgs") ||
        pin->DoesParameterExist(bn, "t_unit_cgs") || pin->DoesParameterExist(bn, "T_unit_K")),
      "<dust> rho_unit_cgs / t_unit_cgs / T_unit_K are no longer read: the dust package now "
      "takes its cgs scales from the single authoritative unit system "
      "(PhysUnits::BuildPhysicalUnits, set by <units> or by the BE IC). Remove them so the "
      "deck cannot claim a normalisation the code does not use.");
  DustModel d;
  d.rho_unit = U.rho_unit;
  d.t_unit = U.time_unit;
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
  pkg->AddParam("T_unit_K", U.temperature_unit);
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

  auto hydro_pkg = pmb->pmy_mesh->packages.Get("Hydro");
  const int nhydro = hydro_pkg->Param<int>("nhydro");
  // audit N7: ask the package what fluid it is, rather than inferring it from the component
  // count. Same value today (glmmhd carries IB1..IB3,IPS past the hydro vars), one less
  // invariant to break silently when a variable is added.
  const bool mhd = (hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd);
  // AUDIT 2026-08-05 (N2) -- THE SUBLIMATION TEMPERATURE MUST COME FROM THE EOS.
  // This derived T as (gamma-1) e/rho with the IDEAL gamma while every production deck sets
  // <hydro> eos = hydrogen, then compared it with T_subl = 1500 K. That is the same defect as
  // audit A1/A6, and it lands in the WORST possible place: 1500 K is inside the H2
  // dissociation zone, exactly where Gamma_1 -> 1.14 and the ideal mapping over-estimates T.
  // Dust therefore sublimated EARLY, removing opacity and grain-surface ionisation with it.
  // Chemistry (chemistry.cpp) and radiation (radiation_moments.cpp) both got this fix; dust
  // was the one module left on the ideal mapping. Ideal-gas path is bit-identical.
  bool use_h2 = false;
  EOSTable::EosTable eos_tab;
  if (mhd) {
    const auto &eos = hydro_pkg->Param<AdiabaticGLMMHDEOS>("eos");
    use_h2 = eos.UseH2Diss();
    if (use_h2) eos_tab = eos.GetEosTable();
  }
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
        // audit N2: tabulated EOS temperature when eos=hydrogen (guarded against eint <= 0,
        // where the table interpolates in log10(eint/rho) and would return NaN); otherwise
        // the ideal mapping, bit-identical to before.
        double T_K;
        if (use_h2 && eint > 0.0) {
          T_K = static_cast<double>(eos_tab.TemperatureK(rho, eint));
          if (T_K < 1.0e-4 * T_unit) T_K = 1.0e-4 * T_unit; // floor (matches ideal path)
        } else {
          Real T_code = gm1 * eint * inv_rho;
          if (T_code < 1.0e-4) T_code = 1.0e-4;
          T_K = static_cast<double>(T_code * T_unit);
        }
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
