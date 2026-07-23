//========================================================================================
// AthenaPK GPU chemistry package — implementation.
//
// The species are AthenaPK passive scalars: scalar_density_i (conserved, advected with
// the hydro) packs rho*x_i, so prim scalar_i = scalar_density_i/rho = abundance x_i.
// Reactions modify x_i in place at fixed rho (operator split). Two networks are
// selectable via <chemistry> network:
//   "H2"            -> minimal H/H2 (NSPEC=2)               [network_h2.hpp]
//   "gow17_reduced" -> reduced H-C-O ionization (NSPEC=5),  [network_gow17_reduced.hpp]
//                      evolves x_e (scalar index 4) which feeds the non-ideal MHD
//                      diffusivities when <diffusion> ...=ionization_chem is set.
//========================================================================================

#include <iostream>
#include <memory>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

#include "chemistry.hpp"
// parthenon (above) defines KOKKOS_INLINE_FUNCTION, so the networks' CHEM_FN/CHEMG_FN
// resolve to the device-callable form. Include AFTER the parthenon headers.
#include "network_h2.hpp"
#include "network_gow17_reduced.hpp"
#include "thermo.hpp" // WS-2 gas thermochemistry (Gamma-Lambda, T_eq)
#include "../eos/adiabatic_glmmhd.hpp" // EOS-table T for reaction rates under eos=hydrogen
#include "../main.hpp" // IDN, IM1..IM3, IEN, IB1..IB3, NHYDRO
#include "../units/physical_units.hpp"       // flagship Phase 1: one unit system
#include "../units/ionization_environment.hpp" // flagship Phase 1: shared CR rate

using namespace parthenon::package::prelude;

namespace Chemistry {

enum class NetworkType { h2, gow17_reduced };

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {
  auto pkg = std::make_shared<StateDescriptor>("chemistry");
  const std::string bn = "chemistry";

  const std::string network = pin->GetOrAddString(bn, "network", "H2");
  PARTHENON_REQUIRE(network == "H2" || network == "gow17_reduced",
                    "chemistry/network must be 'H2' or 'gow17_reduced' (got: " + network + ").");

  // Flagship Phase 1: read the shared code-unit scales from the single authoritative unit
  // system (was: rho_unit/t_unit/T_unit hardcoded here, in ionization.hpp, and derived in
  // radiation.cpp -- audit findings #3/#4). Option A: for collapse_be these are the EXACT BE
  // normalization (rho_unit=5.4668e-19, not the old rounded 5.467e-19), so the network's cgs
  // densities match the dynamics; T_unit is the single T0 (~10.015 K); t_unit = length/v.
  const auto U = PhysUnits::BuildPhysicalUnits(pin);
  const Real rho_unit = U.rho_unit;
  const Real t_unit = U.time_unit;
  const Real T_unit = U.temperature_unit;
  // mu_n (neutral mean molecular weight) is composition, shared with the diffusion
  // ionization model via IonizationEnvironment (which mirrors this same <chemistry> mu_n).
  const Real mu_n = pin->GetOrAddReal(bn, "mu_n", 2.33);
  pkg->AddParam("T_unit_K", T_unit);
  pkg->AddParam("gamma", pin->GetReal("hydro", "gamma"));

  int nspec = 0;
  Real zeta_cr_cgs = 0.0; // resolved CR ionization rate [s^-1] (audit #1: single source of truth)
  if (network == "H2") {
    H2Network net;
    net.rho_unit = rho_unit;
    net.t_unit = t_unit;
    net.mu_n = mu_n;
    net.kgr = pin->GetOrAddReal(bn, "kgr_cgs", 3.0e-17);
    net.xi_cr = pin->GetOrAddReal(bn, "zeta_cr_cgs", 2.0e-16);
    zeta_cr_cgs = net.xi_cr;
    net.kcr = pin->GetOrAddReal(bn, "kcr_factor", 3.0) * net.xi_cr;
    net.x_floor = pin->GetOrAddReal(bn, "x_floor", 1.0e-20);
    net.nsub_max = pin->GetOrAddInteger(bn, "nsub_max", 200);
    net.cfl = pin->GetOrAddReal(bn, "cfl_cool", 0.1);
    pkg->AddParam("network_h2", net);
    pkg->AddParam("network_type", NetworkType::h2);
    nspec = NSPEC_H2;
  } else { // gow17_reduced
    Gow17ReducedNetwork net;
    net.rho_unit = rho_unit;
    net.t_unit = t_unit;
    net.mu_n = mu_n;
    // Cosmic-ray rate from the shared IonizationEnvironment (flagship Phase 1). This is the
    // gow x_e-evolving network that feeds ambipolar_coeff=ionization_chem, so its CR rate
    // must equal the diffusion equilibrium model's -- the environment (and its guard,
    // generalizing audit fix #1) enforces that. Default 1e-16 => bit-identical.
    const auto env = PhysUnits::BuildIonizationEnvironment(pin);
    net.zeta = env.zeta_cr_cgs;
    zeta_cr_cgs = net.zeta;
    net.kgr = pin->GetOrAddReal(bn, "kgr_cgs", 3.0e-17);
    net.x_Ctot = pin->GetOrAddReal(bn, "x_Ctot", 1.6e-4);
    net.xe_floor = pin->GetOrAddReal(bn, "xe_floor", 1.0e-15);
    net.x_floor = pin->GetOrAddReal(bn, "x_floor", 1.0e-20);
    net.nsub_max = pin->GetOrAddInteger(bn, "nsub_max", 400);
    net.cfl = pin->GetOrAddReal(bn, "cfl_cool", 0.1);
    pkg->AddParam("network_gow", net);
    pkg->AddParam("network_type", NetworkType::gow17_reduced);
    nspec = NSPEC_GOW;
  }
  pkg->AddParam("nspec", nspec);

  // Expose the resolved chemistry CR rate for output/back-compat. The chemistry<->diffusion
  // consistency guard (audit fix #1) now lives in BuildIonizationEnvironment, which both the
  // gow network path (above) and the diffusion ionization model call, so a mismatched
  // <diffusion> ion_zeta hard-fails from the shared object rather than a local assertion.
  pkg->AddParam("zeta_cr_cgs", zeta_cr_cgs);

  // Running count of cells whose reaction integration was truncated at nsub_max
  // (under-integrated abundances). Used to warn once on first occurrence.
  pkg->AddParam<>("chem_trunc_total", 0.0, Params::Mutability::Mutable);

  const int nscalars = pin->GetOrAddInteger("hydro", "nscalars", 0);
  PARTHENON_REQUIRE(nscalars >= nspec,
                    "chemistry: <hydro> nscalars (" + std::to_string(nscalars) +
                        ") must be >= number of chemical species (" +
                        std::to_string(nspec) + ") for network '" + network + "'.");

  // --- WS-2 gas thermochemistry (heating/cooling) -----------------------------------
  // thermo=true adds the Gamma-Lambda gas-energy source ON TOP of the abundance network.
  // It needs the gow17 abundances (x_C+, x_CO) and RT (for T_dust = (Er/arad)^{1/4} and
  // because RT owns e_th; the barotropic overwrite would erase any Gamma-Lambda work).
  const bool thermo = pin->GetOrAddBoolean(bn, "thermo", false);
  pkg->AddParam("thermo", thermo);
  if (thermo) {
    PARTHENON_REQUIRE(network == "gow17_reduced",
                      "chemistry/thermo=true requires network=gow17_reduced (needs x_C+, x_CO).");
    PARTHENON_REQUIRE(pin->GetOrAddBoolean("physics", "radiation", false),
                      "chemistry/thermo=true requires <physics> radiation=true (RT owns e_th "
                      "and provides T_dust; the barotropic path would overwrite cooling).");
    // The thermo sub-cycler (AdvanceThermoEnergy) maps e <-> T with the ideal gamma;
    // under eos=hydrogen that mapping is wrong in the dissociation zone (same class as
    // audit F1). Forbid the combination until the EOS table is plumbed through thermo.hpp.
    PARTHENON_REQUIRE(pin->GetString("hydro", "eos") != "hydrogen",
                      "chemistry/thermo=true is not yet compatible with hydro/eos=hydrogen "
                      "(AdvanceThermoEnergy uses the ideal-gas e<->T mapping).");
    ThermoParams tp;
    tp.zeta = pin->GetOrAddReal(bn, "zeta_cr_cgs", 1.0e-16); // consistent with the network
    tp.kgr = pin->GetOrAddReal(bn, "kgr_cgs", 3.0e-17);
    tp.q_cr_erg = pin->GetOrAddReal(bn, "q_cr_ev", 20.0) * thermo_cgs::eV;
    tp.h2_heat_erg = pin->GetOrAddReal(bn, "h2_heat_ev", 0.2) * thermo_cgs::eV;
    tp.pe_eps = pin->GetOrAddReal(bn, "pe_eps", 0.05);
    tp.G0 = pin->GetOrAddReal(bn, "G0", 1.7);
    tp.av0 = pin->GetOrAddReal(bn, "av0", 0.5);
    tp.av_n0 = pin->GetOrAddReal(bn, "av_n0", 1.0e3);
    tp.av_exp = pin->GetOrAddReal(bn, "av_exp", 2.0 / 3.0);
    tp.x_O = pin->GetOrAddReal(bn, "x_O", 3.2e-4);
    tp.co_tau_coeff = pin->GetOrAddReal(bn, "co_tau_coeff", 5.5);
    tp.alpha_gd = pin->GetOrAddReal(bn, "alpha_gd", 3.2e-34);
    pkg->AddParam("thermo_params", tp);
    // Gamma-Lambda [erg cm^-3 s^-1] -> code energy-density rate: * t_unit / e_unit.
    // e_unit (energy-density unit = rho*v^2) from the single authoritative unit system.
    pkg->AddParam("rate_to_code", t_unit / U.energy_density_unit);
    pkg->AddParam("thermo_nsub_max", pin->GetOrAddInteger(bn, "thermo_nsub_max", 200));
    pkg->AddParam("thermo_cfl", pin->GetOrAddReal(bn, "thermo_cfl", 0.1));
    // Internal-energy-density floor for the thermal integrator (code units); chemistry-local
    // so it does not collide with hydro/pfloor's own default.
    pkg->AddParam("thermo_efloor", pin->GetOrAddReal(bn, "thermo_efloor_code", 1.0e-15));
  }

  return pkg;
}

TaskStatus ReactScalars(MeshData<Real> *md, const Real dt) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  auto pkg = pmb->pmy_mesh->packages.Get("chemistry");
  const auto ntype = pkg->Param<NetworkType>("network_type");

  // scalar_density_i are components [nhydro, nhydro+nscalars) of the "cons" field.
  auto hydro_pkg = pmb->pmy_mesh->packages.Get("Hydro");
  const int nhydro = hydro_pkg->Param<int>("nhydro");
  const bool mhd = (nhydro > NHYDRO); // glmmhd carries IB1..IB3(,IPS) past the hydro vars

  // F1 fix (audit 2026-07-18): under eos=hydrogen the reaction temperature must come from
  // the EOS table -- mu varies with H2 dissociation/ionization, so T != gm1*eint/rho there
  // (the ideal identity over-estimates T, hence x_e, in the ~2000 K dissociation zone).
  // Mirrors the PrecomputeNonidealEta idiom. Ideal-gas path is bit-identical.
  bool use_h2 = false;
  EOSTable::EosTable eos_tab;
  if (mhd) {
    const auto &eos = hydro_pkg->Param<AdiabaticGLMMHDEOS>("eos");
    use_h2 = eos.UseH2Diss();
    if (use_h2) eos_tab = eos.GetEosTable();
  }

  // WS-2 thermochemistry needs T_dust from the M1 field -> also pack rad.Er.
  const bool thermo = pkg->Param<bool>("thermo");
  std::vector<std::string> names{"cons"};
  if (thermo) names.push_back("rad.Er");
  parthenon::PackIndexMap imap;
  auto pack = md->PackVariables(names, imap);
  const int ic = imap["cons"].first; // cons base: rho = pack(ic+IDN)
  const int is = ic + nhydro;        // first scalar_density component
  const int iEr = thermo ? imap["rad.Er"].first : -1;

  // Count cells whose sub-stepped integration hit nsub_max before covering dt (their
  // abundances are silently under-integrated); warn on the first occurrence per run.
  int ntrunc = 0;

  if (ntype == NetworkType::h2) {
    const H2Network net = pkg->Param<H2Network>("network_h2");
    parthenon::par_reduce(
        parthenon::loop_pattern_mdrange_tag, "Chem::ReactH2", parthenon::DevExecSpace(),
        0, pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, int &ltrunc) {
          const Real rho = pack(b, ic + IDN, k, j, i);
          const Real inv_rho = 1.0 / rho;
          double y[NSPEC_H2];
          for (int s = 0; s < NSPEC_H2; ++s)
            y[s] = static_cast<double>(pack(b, is + s, k, j, i) * inv_rho);
          ltrunc +=
              net.integrate_cell(y, static_cast<double>(rho), static_cast<double>(dt));
          for (int s = 0; s < NSPEC_H2; ++s)
            pack(b, is + s, k, j, i) = rho * static_cast<Real>(y[s]);
        },
        Kokkos::Sum<int>(ntrunc));
  } else { // gow17_reduced — needs the gas temperature for recombination rates
    const Gow17ReducedNetwork net = pkg->Param<Gow17ReducedNetwork>("network_gow");
    const Real gm1 = pkg->Param<Real>("gamma") - 1.0;
    const Real T_unit = pkg->Param<Real>("T_unit_K");
    // WS-2 thermochemistry (captured by value; all POD -> device-safe). Inert if !thermo.
    ThermoParams tp;
    Real arad = 0.0, rate_to_code = 0.0, thermo_cfl = 0.1, thermo_efloor = 0.0;
    int thermo_nsub = 1;
    if (thermo) {
      tp = pkg->Param<ThermoParams>("thermo_params");
      arad = pmb->pmy_mesh->packages.Get("radiation")->Param<Real>("arad");
      rate_to_code = pkg->Param<Real>("rate_to_code");
      thermo_nsub = pkg->Param<int>("thermo_nsub_max");
      thermo_cfl = pkg->Param<Real>("thermo_cfl");
      thermo_efloor = pkg->Param<Real>("thermo_efloor");
    }
    parthenon::par_reduce(
        parthenon::loop_pattern_mdrange_tag, "Chem::ReactGow", parthenon::DevExecSpace(),
        0, pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, int &ltrunc) {
          const Real rho = pack(b, ic + IDN, k, j, i);
          const Real inv_rho = 1.0 / rho;
          // Gas temperature from the conserved energy (code units: p = rho*T_code).
          const Real m1 = pack(b, ic + IM1, k, j, i);
          const Real m2 = pack(b, ic + IM2, k, j, i);
          const Real m3 = pack(b, ic + IM3, k, j, i);
          Real eint = pack(b, ic + IEN, k, j, i) - 0.5 * (m1 * m1 + m2 * m2 + m3 * m3) * inv_rho;
          if (mhd) {
            const Real b1 = pack(b, ic + IB1, k, j, i);
            const Real b2 = pack(b, ic + IB2, k, j, i);
            const Real b3 = pack(b, ic + IB3, k, j, i);
            eint -= 0.5 * (b1 * b1 + b2 * b2 + b3 * b3); // HL: P_mag = B^2/2
          }
          double T_K;
          if (use_h2 && eint > 0.0) {
            // EOS-table temperature (guarded against eint <= 0: the table interpolates
            // in log10(eint/rho) -- fall through to the floored ideal branch there).
            T_K = static_cast<double>(eos_tab.TemperatureK(rho, eint));
            if (T_K < 1.0e-4 * T_unit) T_K = 1.0e-4 * T_unit; // floor (matches ideal path)
          } else {
            Real T_code = gm1 * eint * inv_rho;
            if (T_code < 1.0e-4) T_code = 1.0e-4; // floor
            T_K = static_cast<double>(T_code * T_unit);
          }

          double y[NSPEC_GOW];
          for (int s = 0; s < NSPEC_GOW; ++s)
            y[s] = static_cast<double>(pack(b, is + s, k, j, i) * inv_rho);
          ltrunc += net.integrate_cell(y, static_cast<double>(rho), T_K,
                                       static_cast<double>(dt));
          for (int s = 0; s < NSPEC_GOW; ++s)
            pack(b, is + s, k, j, i) = rho * static_cast<Real>(y[s]);

          // WS-2: gas-energy source Gamma-Lambda (after abundances; RT coupling ran first).
          if (thermo) {
            const double n_H = net.nH(static_cast<double>(rho));
            double x_H = 1.0 - 2.0 * y[gH2] - y[gHp];
            if (x_H < 0.0) x_H = 0.0;
            const Real Er = pack(b, iEr, k, j, i);
            const double T_dust =
                (Er > 0.0 ? std::pow(static_cast<double>(Er / arad), 0.25) : 0.0) * T_unit;
            int lt = 0;
            const double e_new = AdvanceThermoEnergy(
                tp, static_cast<double>(rho), static_cast<double>(eint), n_H, x_H, y[gCp],
                y[gCO], T_dust, static_cast<double>(gm1), static_cast<double>(T_unit),
                static_cast<double>(rate_to_code), static_cast<double>(dt), thermo_nsub,
                static_cast<double>(thermo_cfl), static_cast<double>(thermo_efloor), &lt);
            pack(b, ic + IEN, k, j, i) += static_cast<Real>(e_new) - eint;
            ltrunc += lt;
          }
        },
        Kokkos::Sum<int>(ntrunc));
  }

  if (ntrunc > 0) {
    const Real prev = pkg->Param<Real>("chem_trunc_total");
    pkg->UpdateParam("chem_trunc_total", prev + static_cast<Real>(ntrunc));
    if (prev == 0.0) {
      std::cout << "### WARNING Chemistry: " << ntrunc
                << " cell(s) hit nsub_max before covering the full dt; abundances in "
                   "those cells are under-integrated. Consider raising "
                   "chemistry/nsub_max or cfl_cool. (Warning printed once per run; "
                   "total truncations accumulate in the chem_trunc_total param.)"
                << std::endl;
    }
  }
  return TaskStatus::complete;
}

void AddChemistryTasks(TaskCollection &tc, Mesh *pmesh, const Real dt) {
  using namespace parthenon;
  TaskID none(0);

  auto partitions = pmesh->GetDefaultBlockPartitions();
  const int num_partitions = partitions.size();
  TaskRegion &region = tc.AddRegion(num_partitions);
  for (int i = 0; i < num_partitions; ++i) {
    TaskList &tl = region[i];
    auto &md = pmesh->mesh_data.Add("base", partitions[i]);
    tl.AddTask(none, ReactScalars, md.get(), dt);
  }
}

} // namespace Chemistry
