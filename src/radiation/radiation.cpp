//========================================================================================
// AthenaPK - M1 moment radiation transport package.
// Ported from Artemis (LANL, BSD) src/radiation/moments. Licensed BSD 3-Clause.
//
// INCREMENT 1 (scaffold): package registration only. Registers the M1 fields and
// constant/opacity params. No transport yet — AddRadiationTasks is a no-op.
//========================================================================================

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

#include "../eos/eos_table.hpp"
#include "../units/physical_units.hpp" // flagship Phase 1: one unit system (audit #4)
#include "radiation.hpp"
#include "radiation_closure.hpp"
#include "radiation_groups.hpp"
#include "radiation_opacity.hpp"

namespace Radiation {

using namespace parthenon::package::prelude;

namespace {
// Force instantiation of the M1 closure templates so they are compile-checked now
// (host instantiation; device kernels in increment 2b exercise __device__ paths).
[[maybe_unused]] void rad_closure_compile_check() {
  const Real f = 0.5;
  const auto chi = ThriceEddingtonFactor<Closure::M1>(f);
  const auto P = EddingtonTensor<Closure::M1>({0.1, 0.2, 0.3});
  const auto lam = WaveSpeed<Closure::M1>(0.5, f);
  const auto nf = NormalizeFlux(0.1, 0.2, 0.3);
  const auto ff = FleckFactor(1.0, 1.0, 1.0);
  (void)chi;
  (void)P;
  (void)lam;
  (void)nf;
  (void)ff;
}
} // namespace

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {
  auto pkg = std::make_shared<StateDescriptor>("radiation");
  const std::string bn = "radiation";

  // --- Code-unit calibration -------------------------------------------------
  // The radiation energy density Er lives in the gas energy-density unit
  //   e_unit = rho_unit * v_unit^2  (= the pressure unit, shared by both codes),
  // and the temperature unit follows from c_s,iso = 1 at the reference temperature:
  //   T_unit = mu * m_H * v_unit^2 / k_B   (=> T_unit ~ 10 K for the FHC IC).
  // From these, the physical radiation constant / speed of light / dust opacity map
  // into code units with NO free knobs:
  //   arad_code   = a_R * T_unit^4 / e_unit            (so E_eq = arad_code * T_code^4)
  //   c_code      = c_light / v_unit                   (~1.58e6 => RSLA is mandatory)
  //   kappa_conv  = rho_unit * l_unit                  (kappa_code = kappa_phys*kappa_conv)
  // Defaults are the FHC shared units; a test can override e.g. arad_code=c_code=1.
  // Flagship Phase 1 (roots audit finding #4): the temperature/energy/opacity/light units
  // come from the SINGLE authoritative unit system, not a radiation-local mu. The RT
  // temperature calibration is mu_thermal = 2.29 (=> T0 = 10.015 K), the SAME T0 chemistry
  // and the non-ideal-MHD ionization model use. Production previously overrode <radiation>
  // mu = 2.33 here, giving T0 ~= 10.19 K and shifting a_R*T0^4/e0 (and every opacity
  // temperature) by ~7% relative to the rest of the microphysics -- that split is now
  // structurally impossible. mu_thermal (T-unit) is deliberately NOT the neutral mu_n = 2.33
  // used for number densities (that lives in IonizationEnvironment).
  const auto U = PhysUnits::BuildPhysicalUnits(pin);
  const Real rho_unit = U.rho_unit;   // g/cm^3
  const Real v_unit = U.v_unit;       // cm/s
  const Real l_unit = U.length_unit;  // cm
  const Real T_unit = U.temperature_unit; // K
  const Real e_unit = U.energy_density_unit; // erg/cm^3
  const Real arad_default = U.arad_code();
  const Real c_default = U.c_code();
  const Real kappa_conv = U.opacity_unit;
  // A legacy <radiation> mu is now IGNORED (the T-unit is single-sourced). Warn loudly so a
  // deck carrying the old mu=2.33 override is not silently misread as still active.
  if (pin->DoesParameterExist(bn, "mu")) {
    const Real mu_legacy = pin->GetReal(bn, "mu");
    if (parthenon::Globals::my_rank == 0 &&
        std::abs(mu_legacy - U.mu_thermal) > 1.0e-6 * U.mu_thermal) {
      std::cout << "### WARNING Radiation: <radiation> mu=" << mu_legacy
                << " is IGNORED. The radiation temperature unit is now single-sourced from "
                   "PhysicalUnits (mu_thermal=" << U.mu_thermal << ", T0=" << T_unit
                << " K). Set <units> mu_thermal to change it." << std::endl;
    }
  }

  //   c      : speed of light  (code units; default = calibrated c_light/v_unit)
  //   arad   : radiation constant a_R (so E_eq = arad * T^4; default = calibrated)
  //   creduc : reduced-speed-of-light factor (RSLA); chat = c / creduc (>=1)
  const Real c_code = pin->GetOrAddReal(bn, "c_code", c_default);
  const Real arad_code = pin->GetOrAddReal(bn, "arad_code", arad_default);
  const Real creduc = pin->GetOrAddReal(bn, "creduc", 1.0);
  PARTHENON_REQUIRE(creduc >= 1.0, "radiation/creduc must be >= 1 (chat = c/creduc <= c).");
  pkg->AddParam("c", c_code);
  pkg->AddParam("arad", arad_code);
  pkg->AddParam("chat", c_code / creduc);
  pkg->AddParam("creduc", creduc);
  pkg->AddParam("T_unit", T_unit);

  // --- Gray opacity model (code units) ---------------------------------------
  // "constant": kappa_a_code / kappa_s_code placeholders (per unit mass).
  // "dust":     Bell & Lin (1994) low-T ice-grain law kappa = k0 * T_phys^2 [cm^2/g],
  //             converted to code units internally (this is what traps the FHC).
  OpacityParams op;
  const std::string opmodel = pin->GetOrAddString(bn, "opacity_model", "constant");
  PARTHENON_REQUIRE(opmodel == "constant" || opmodel == "dust" || opmodel == "belllin" ||
                        opmodel == "tabulated",
                    "radiation/opacity_model must be 'constant', 'dust', 'belllin', or "
                    "'tabulated'.");
  op.model = (opmodel == "tabulated")
                 ? OpacityModel::tabulated
                 : (opmodel == "belllin")
                       ? OpacityModel::belllin
                       : (opmodel == "dust") ? OpacityModel::dust : OpacityModel::constant;
  op.kappa_a0 = pin->GetOrAddReal(bn, "kappa_a_code", 0.0);
  op.kappa_s0 = pin->GetOrAddReal(bn, "kappa_s_code", 0.0);
  op.dust_k0_cgs = pin->GetOrAddReal(bn, "dust_kappa0_cgs", 2.0e-4);
  op.dust_Tmax = pin->GetOrAddReal(bn, "dust_Tmax_K", 1.5e3);
  op.kappa_s_dust = pin->GetOrAddReal(bn, "kappa_s_code", 0.0);
  op.T_unit = T_unit;
  op.rho_unit = rho_unit;
  op.kappa_conv = kappa_conv;
  // Planck/Rosseland split: kappa_P = planck_ross_ratio * kappa_R for the emission term.
  // Default 1.0 => bit-identical to the single-opacity behavior.
  op.planck_ross_ratio = pin->GetOrAddReal(bn, "planck_ross_ratio", 1.0);
  PARTHENON_REQUIRE(op.planck_ross_ratio > 0.0,
                    "radiation/planck_ross_ratio must be > 0.");
  // Regime-skip-robust Bell&Lin walk (default false = bit-identical; fixes the off-track
  // low-rho/high-T kappa discontinuity). Only meaningful for opacity_model=belllin.
  op.bell_lin_fix_regime_skip =
      pin->GetOrAddBoolean(bn, "bell_lin_fix_regime_skip", false);
  // Tabulated model: load the offline frequency-resolved gray magnitude table (true Planck +
  // Rosseland means; self-consistent, replacing the planck_ross_ratio fudge). The per-group
  // multipliers from the SAME file feed GroupOpacityTable below.
  const std::string opac_table_file = pin->GetOrAddString(
      bn, "opacity_table_file",
      "/beegfs/u/bbg6470/athenapk/src/radiation/opacity_table.bin");
  if (op.model == OpacityModel::tabulated) op.table.Load(opac_table_file, rho_unit, T_unit);
  pkg->AddParam("opacity", op);
  // Back-compat scalar params (constant model still reads these directly).
  pkg->AddParam("kappa_a", op.kappa_a0);
  pkg->AddParam("kappa_s", op.kappa_s0);

  // --- Closure model ---------------------------------------------------------
  const std::string closure = pin->GetOrAddString(bn, "closure", std::string("M1"));
  PARTHENON_REQUIRE(closure == "M1" || closure == "P1",
                    "radiation/closure must be M1 or P1 (got: " + closure + ").");
  pkg->AddParam("closure", closure);

  // --- Multigroup frequency structure (Phase 4) ------------------------------
  // n_group=1 (default) => a single group spanning [0,inf) == the gray path, bit-identical.
  // n_group>1 partitions the spectrum into log-spaced groups [nu_min,nu_max] Hz (edges 0 and
  // +inf on the ends so no energy is lost); each group carries its own (Er_g,Fr_g) M1 system,
  // the closure is frequency-independent (reused per group), and the groups couple ONLY
  // through the matter temperature in MatterCoupling. Monochromatic opacity nu-dependence
  // kappa_nu ~ (nu/nu_ref)^opacity_nu_index (index 0 => gray-per-group).
  const int n_group = pin->GetOrAddInteger(bn, "n_group", 1);
  PARTHENON_REQUIRE(n_group >= 1 && n_group <= MAX_GROUP,
                    "radiation/n_group must be in [1, MAX_GROUP=8].");
  const Real nu_min_hz = pin->GetOrAddReal(bn, "nu_min_hz", 1.0e12); // ~40 K blackbody peak
  const Real nu_max_hz = pin->GetOrAddReal(bn, "nu_max_hz", 1.0e15); // ~far-UV
  RadGroups groups = BuildRadGroups(n_group, nu_min_hz, nu_max_hz);
  // Monochromatic dust opacity spectral index beta in kappa_nu = kappa_gray*(nu/nu_ref)^beta.
  // The per-group Planck/Rosseland means are the band-limited averages of kappa_nu with the
  // correct weights (B_nu for emission/absorption, dB_nu/dT for flux), normalized so the
  // full-spectrum means reproduce kappa_gray. beta=0 => every group mean == kappa_gray =>
  // multigroup reduces EXACTLY to n_group copies of gray (equivalence gate). `opacity_beta`
  // is the key; `opacity_nu_index` (the earlier representative-frequency exponent) is accepted
  // as a fallback so existing decks keep working.
  const Real beta_fallback = pin->GetOrAddReal(bn, "opacity_nu_index", 0.0);
  DustOpacityModel dmodel;
  dmodel.beta = pin->GetOrAddReal(bn, "opacity_beta", beta_fallback);
  dmodel.nu_break_hz = pin->GetOrAddReal(bn, "opacity_nu_break_hz", 1.0e14);
  dmodel.T_sub_lo = pin->GetOrAddReal(bn, "opacity_T_sub_lo_K", 1400.0);
  dmodel.T_sub_hi = pin->GetOrAddReal(bn, "opacity_T_sub_hi_K", 1600.0);
  pkg->AddParam("groups", groups);
  pkg->AddParam("n_group", n_group);
  pkg->AddParam("opacity_beta", dmodel.beta);
  pkg->AddParam("dust_opacity_model", dmodel);
  // Precompute the per-group Semenov-class band-mean multipliers on a log(T) grid (device
  // table) so the matter coupling does an O(1) interpolation instead of a per-cell quadrature.
  // Inactive (multiplier == 1) for gray/beta=0, so the equivalence gate is untouched.
  const int op_nT = pin->GetOrAddInteger(bn, "opacity_table_nT", 256);
  const Real op_Tmin = pin->GetOrAddReal(bn, "opacity_table_Tmin_K", 3.0);
  const Real op_Tmax = pin->GetOrAddReal(bn, "opacity_table_Tmax_K", 1.0e6);
  // Tabulated model: per-group multipliers come from the SAME frequency-resolved file (real
  // <psi>_band/<psi>_full ratios) instead of the analytic DustOpacityModel band means. Gray
  // (n_group=1) needs only the gray magnitude table (loaded above) -> inactive multiplier table,
  // so DON'T read the file's per-group block (and don't require the file's ng to match 1).
  if (op.model == OpacityModel::tabulated && n_group > 1)
    pkg->AddParam("optable", BuildGroupOpacityTableFromFile(opac_table_file, n_group));
  else
    pkg->AddParam("optable", BuildGroupOpacityTable(dmodel, groups, op_Tmin, op_Tmax, op_nT));

  // --- Transport controls (increment 2b) -------------------------------------
  pkg->AddParam("cfl", pin->GetOrAddReal(bn, "cfl", 0.4));       // radiation CFL number
  pkg->AddParam("efloor", pin->GetOrAddReal(bn, "efloor", 1.0e-15)); // Er floor

  // Spatial reconstruction of (Er, Fr) at cell faces for the M1 transport flux (WS-3b):
  //   "dc"  = donor cell / piecewise-constant (1st order, default; bit-identical to pre-3b).
  //   "plm" = piecewise-linear + minmod (2nd order; +-2 stencil, needs nghost >= 2).
  const std::string rad_recon = pin->GetOrAddString(bn, "reconstruction", std::string("dc"));
  PARTHENON_REQUIRE(rad_recon == "dc" || rad_recon == "plm",
                    "radiation/reconstruction must be 'dc' or 'plm' (got: " + rad_recon + ").");
  const bool rad_plm = (rad_recon == "plm");
  if (rad_plm) {
    const int nghost = pin->GetInteger("parthenon/mesh", "nghost");
    PARTHENON_REQUIRE(nghost >= 2,
                      "radiation/reconstruction=plm needs parthenon/mesh/nghost >= 2 "
                      "(the +-2 face stencil); increase nghost.");
  }
  pkg->AddParam("rad_recon_plm", rad_plm);

  // --- Matter coupling (increment 3) -----------------------------------------
  // When the radiation package is active it OWNS the gas thermal energy: the
  // absorption/emission source drives gas e_th toward the radiation field
  // (E_eq = arad T^4), replacing collapse_be's barotropic floor. Needs the gas
  // adiabatic index to map internal energy <-> temperature: for the FHC code
  // units (c_s,iso = 1 at 10 K), pressure p = rho*T_code with T_code = 1 <=> 10 K,
  // so e_int = rho*T/(gamma-1) and the heat capacity per volume Cv = rho/(gamma-1).
  pkg->AddParam("matter_coupling", pin->GetOrAddBoolean(bn, "matter_coupling", true));
  pkg->AddParam("gamma", pin->GetReal("hydro", "gamma"));
  // H2-dissociation general EOS: when active, RT owns the gas energy and must use the
  // SAME Saha EOS as the hydro to map internal-energy <-> temperature (else RT and hydro
  // disagree on T and the dissociation buffering is inconsistent). Mirror the hydro-block
  // keys (defaults = FHC-calibrated). T(e) uses T_code = T_phys[K]/T_unit.
  const bool use_h2diss = (pin->GetOrAddString("hydro", "eos", "adiabatic") == "hydrogen");
  pkg->AddParam("use_h2diss", use_h2diss);
  EOSTable::EosTable eos_tab;
  if (use_h2diss) {
    eos_tab.Load(pin->GetOrAddString("hydro", "eos_table_file",
                                     "/beegfs/u/bbg6470/athenapk/src/eos/eos_table.bin"));
  }
  pkg->AddParam("eos_tab", eos_tab);
  pkg->AddParam("tfloor", pin->GetOrAddReal(bn, "tfloor_code", 1.0e-3)); // T floor (code)
  pkg->AddParam("inner_iteration_max", pin->GetOrAddInteger(bn, "inner_iteration_max", 100));
  pkg->AddParam("inner_iteration_tol", pin->GetOrAddReal(bn, "inner_iteration_tol", 1.0e-8));

  // --- Fields: M1 two-moment radiation (Er, Fr1, Fr2, Fr3) -------------------
  // Independent (evolved state: communicated, refined, checkpointed) + WithFluxes (our
  // own HLL flux arrays) + FillGhost. The OperatorSplit user-flag is the key: AthenaPK's
  // hydro integrator packs {WithFluxes,Cell} and would otherwise sweep these into its
  // low-storage update (and zero them via the un-seeded u1 register). The vendored
  // parthenon FluxDivergence / UpdateWithFluxDivergence now Exclude OperatorSplit, so
  // hydro leaves radiation alone and AddRadiationTasks advances it operator-split.
  const parthenon::MetadataFlag OperatorSplit =
      parthenon::Metadata::FlagNameExists("OperatorSplit")
          ? parthenon::Metadata::GetUserFlag("OperatorSplit")
          : parthenon::Metadata::AddUserFlag("OperatorSplit");
  std::vector<parthenon::MetadataFlag> flags{Metadata::Cell, Metadata::Independent,
                                             Metadata::FillGhost, Metadata::WithFluxes,
                                             OperatorSplit};
  Metadata m(flags);
  m.RegisterRefinementOps<parthenon::refinement_ops::ProlongatePiecewiseConstant,
                          parthenon::refinement_ops::RestrictAverage>();
  pkg->AddField<rad::Er>(m);
  pkg->AddField<rad::Fr1>(m);
  pkg->AddField<rad::Fr2>(m);
  pkg->AddField<rad::Fr3>(m);
  // Multigroup: register the extra groups' moments by string name (same Metadata as the
  // gray fields). Group 0 is the four fields above (original names) so n_group=1 adds
  // nothing here -> registration is byte-for-byte the gray layout.
  for (int g = 1; g < n_group; ++g)
    for (const auto &nm : GroupFieldNames(g)) pkg->AddField(nm, m);

  return pkg;
}

} // namespace Radiation
