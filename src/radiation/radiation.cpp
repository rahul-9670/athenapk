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

  // AUDIT 2026-08-05 (N4): <radiation> rho_unit_cgs / v_unit_cgs / length_unit_cgs are DEAD
  // KEYS and have been since the Phase-1 unit consolidation -- the base scales come from
  // PhysUnits::BuildPhysicalUnits, which reads block <units> (or derives them from the BE IC).
  // 165 of the 295 decks in the tree still carry them, where they read as authoritative
  // (counted 2026-08-05; all 165 occurrences are in <radiation>). Worse, the
  // hard-fail that protects the BE normalisation only inspects <units>, so writing the
  // override in <radiation> gets SILENCE rather than the intended error. Warn instead of
  // failing: the values in every existing deck happen to agree with the derived ones to
  // rounding, and hard-failing would break decks belonging to queued jobs.
  if (parthenon::Globals::my_rank == 0) {
    const char *keys[3] = {"rho_unit_cgs", "v_unit_cgs", "length_unit_cgs"};
    const Real derived[3] = {rho_unit, v_unit, l_unit};
    for (int q = 0; q < 3; ++q) {
      if (!pin->DoesParameterExist(bn, keys[q])) continue;
      const Real given = pin->GetReal(bn, keys[q]);
      std::cout << "### WARNING Radiation: <radiation> " << keys[q] << " = " << given
                << " is IGNORED (dead key). The base scales are single-sourced from "
                   "PhysicalUnits; this run uses "
                << derived[q] << " (relative difference "
                << (derived[q] != 0.0 ? (given - derived[q]) / derived[q] : 0.0)
                << "). Remove it from the deck, or set it in <units> -- where, for a "
                   "collapse_be problem, it is correctly REJECTED rather than ignored."
                << std::endl;
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
  // Default flipped false -> true 2026-08-08; see OpacityParams in radiation_opacity.hpp for the
  // measurement. Bit-identical below ~100 K (the whole ensemble epoch); up to 8.59 decades of
  // kappa error avoided once the run is hot at low density.
  op.bell_lin_fix_regime_skip =
      pin->GetOrAddBoolean(bn, "bell_lin_fix_regime_skip", true);
  // Tabulated model: load the offline frequency-resolved gray magnitude table (true Planck +
  // Rosseland means; self-consistent, replacing the planck_ross_ratio fudge). The per-group
  // multipliers from the SAME file feed GroupOpacityTable below.
  // AUDIT N13, DEFAULT CORRECTED 2026-08-08. This default used to be `opacity_table.bin`, which
  // is a LEGACY v1 file: kappa stored in the generator's own rounded code units. Measured on the
  // two files directly (all three arrays kP/kR/ks, every one of 60x120 = 7200 nodes, min == max):
  // v1 is UNIFORMLY +0.13538 % larger than what this run's opacity_unit gives, because the
  // generator hardcoded rho0 = 5.467e-19 and l0 = 2.81e16 (product 0.01536227) where the exact
  // BE normalization gives 0.0153415. The two files are otherwise the SAME PHYSICS: v1 divided by
  // (v2_cgs x the generator's unit) is 1.0000000000 to ten digits, so this switch changes that one
  // factor and nothing else.
  //
  // The production decks (root_ladder/fhc_rootladder.in:203 and all 24 ensemble members) already
  // name opacity_table_v2.bin explicitly and were never affected. The hazard this fixes is the
  // SILENT one: a new deck that sets `opacity_model = tabulated` without naming a file inherited
  // the biased table with nothing in the deck to show it. 28 older decks (wp1_sts, wp2_creduc,
  // wp18_seed_ensemble, prod_flagship_test, convergence_ladder) still name the v1 file EXPLICITLY
  // and are deliberately left alone -- they are completed harnesses whose results must stay
  // reproducible, and their explicit path still resolves to the same file it always did.
  const std::string opac_table_file = pin->GetOrAddString(
      bn, "opacity_table_file",
      "/beegfs/u/bbg6470/athenapk/src/radiation/opacity_table_v2.bin");
  if (op.model == OpacityModel::tabulated) {
    // AUDIT N13 (2026-08-05), FIXED. A v2 table stores kappa in cgs and Load() converts it
    // with THIS run's opacity_unit, so the file is independent of any particular IC. A legacy
    // v1 table stores kappa already in the generator's own rounded code units and is used
    // verbatim -- bit-identical to the historical behaviour, and carrying the historical bias.
    // Warn in that case, quantifying it, so an old table can never be mistaken for a fixed one.
    op.table.Load(opac_table_file, rho_unit, T_unit, kappa_conv);
    if (parthenon::Globals::my_rank == 0 && op.table.used_legacy_code_units) {
      constexpr Real kGenRhoUnit = 5.467e-19; // gen_opacity_table.py RHO_UNIT, format v1
      constexpr Real kGenLenUnit = 2.81e16;   // gen_opacity_table.py LEN_UNIT, format v1
      const Real unit_ratio = kGenRhoUnit * kGenLenUnit / kappa_conv;
      std::cout << "### WARNING Radiation (audit N13): " << opac_table_file
                << " is a LEGACY (v1) opacity table: kappa is stored in the generator's own "
                   "code units and is used verbatim. This run's opacity_unit = "
                << kappa_conv << " vs the generator's " << kGenRhoUnit * kGenLenUnit
                << " => kappa is " << 100.0 * (unit_ratio - 1.0)
                << "% off. Regenerate with gen_opacity_table.py (which now emits cgs) to "
                   "remove the bias; the v2 file is correct for any normalization."
                << std::endl;
    }
  }
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
  PARTHENON_REQUIRE(!(n_group > 1) || (nu_min_hz > 0.0 && nu_max_hz > nu_min_hz),
                    "radiation: multigroup needs 0 < nu_min_hz < nu_max_hz.");
  // audit N11 (2026-08-05): n_group=2 has one interior edge, so nu_max_hz cannot be honoured.
  if (n_group == 2 && parthenon::Globals::my_rank == 0)
    std::cout << "### WARNING Radiation: n_group = 2 has a single interior group edge, which is "
                 "placed at nu_min_hz = "
              << nu_min_hz << " Hz. radiation/nu_max_hz = " << nu_max_hz
              << " is IGNORED. Use n_group >= 3 to span [nu_min_hz, nu_max_hz]." << std::endl;
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

  // Diagnostics, promoted from environment variables to deck keys on 2026-08-10 (they were
  // RAD_DISABLE_TRANSPORT / RAD_PRINT_NSUB). Both default false, so behaviour is unchanged;
  // the point is that setting one now leaves a record in the deck and in the restart file,
  // instead of living only in the submitting shell's environment.
  // disable_transport makes M1 transport a NO-OP -- it is a diagnostic, not a physics switch.
  pkg->AddParam("disable_transport", pin->GetOrAddBoolean(bn, "disable_transport", false));
  pkg->AddParam("print_nsub", pin->GetOrAddBoolean(bn, "print_nsub", false));
  pkg->AddParam("gamma", pin->GetReal("hydro", "gamma"));
  // H2-dissociation general EOS: when active, RT owns the gas energy and must use the
  // SAME Saha EOS as the hydro to map internal-energy <-> temperature (else RT and hydro
  // disagree on T and the dissociation buffering is inconsistent). Mirror the hydro-block
  // keys (defaults = FHC-calibrated). T(e) uses T_code = T_phys[K]/T_unit.
  const bool use_h2diss = (pin->GetOrAddString("hydro", "eos", "adiabatic") == "hydrogen");
  pkg->AddParam("use_h2diss", use_h2diss);
  EOSTable::EosTable eos_tab;
  if (use_h2diss) {
    // Same units as the hydro package's load (audit N14): a v2 cgs table is converted with
    // this run's scales, a legacy table is consumed verbatim. hydro.cpp does the warning.
    // DEFAULT MUST MATCH hydro.cpp:996-998 EXACTLY. Parthenon records the default passed to
    // GetOrAddString and hard-aborts ("has at least two inconsistent default values") when the
    // same (block, key) is registered twice with different defaults -- and it does so even when
    // the deck sets the key explicitly, so an explicit production deck is NOT immune.
    //
    // 2026-08-10: this was live on HEAD. Commit a97b6673 (2026-08-08, audit N14) flipped the
    // hydro side to the v2 cgs table and left this one at the legacy v1 path, which aborted
    // every `hydro/eos = hydrogen` + `physics/radiation = true` run at startup -- the flagship
    // configuration. hydro.cpp:994 had already written down the rule this broke ("Paired with
    // the opacity default in radiation.cpp; changing one without the other would leave half the
    // hazard open"); the opacity pair was updated, this EOS pair was not. Neither GPU binary
    // carried the defect (both predate a97b6673 -- verified: neither contains the string
    // "src/eos/eos_table_v2.bin"), so no production result is affected; only a rebuild from
    // HEAD would have been, and it would have died in the first seconds.
    eos_tab.Load(pin->GetOrAddString("hydro", "eos_table_file",
                                     "/beegfs/u/bbg6470/athenapk/src/eos/eos_table_v2.bin"),
                 U.rho_unit, U.v_unit);
  }
  pkg->AddParam("eos_tab", eos_tab);

  // --- A3 (2026-08-10): the opacity-table domain must not be a SILENT ceiling ---------
  // OpacityTable::bilin (radiation_opacity.hpp:62-77) clamps BOTH the index and the
  // interpolation weight, so a lookup past an edge returns the edge opacity: no warning,
  // no error, no NaN, nothing in the output to show it happened. That made a real ceiling
  // invisible -- opacity_table_v2.bin stops at rho = 1e-2 g/cm^3 and T = 1e5 K, while
  // eos_table_hires_v2.bin covers rho <= 1 g/cm^3 and T <= 3.0e5 K. The EOS is therefore
  // valid through second collapse and the opacity is not, and nothing said so.
  //
  // Measured cost of the clamp at corners v2 could not represent (bilin ported to python,
  // compared against gen_opacity_table.group_opacities directly):
  //   rho = 1e0  g/cm^3, T = 2e3 K : kappa_R 3.74e0 vs 8.00e1 exact  (-95.3 %)
  //   rho = 1e-6 g/cm^3, T = 2e5 K : kappa_R 4.84e1 vs 8.39e0 exact  (+477 %)
  // The T ceiling is the broader of the two: it binds at EVERY density, not just the
  // second-core densities, so it is not only a deep-collapse concern.
  //
  // Fix has two halves. (1) opacity_table_v3.bin extends the grid to the full EOS domain
  // (see src/eos/TABLE_PROVENANCE.md). (2) this report, so the domain is in every run log
  // and a short table announces itself at startup instead of at analysis time. WARN, never
  // abort: a hard failure here would be restart-hostile for every completed deck that
  // legitimately names a narrower table.
  if (op.model == OpacityModel::tabulated && parthenon::Globals::my_rank == 0) {
    const auto &ot = op.table;
    const Real op_rlo = std::pow(10.0, ot.lr0_);
    const Real op_rhi = std::pow(10.0, ot.lr0_ + (ot.nr_ - 1) * ot.dlr_);
    const Real op_tlo = std::pow(10.0, ot.lT0_);
    const Real op_thi = std::pow(10.0, ot.lT0_ + (ot.nT_ - 1) * ot.dlT_);
    std::cout << "### Radiation opacity table: " << opac_table_file << "\n"
              << "###   domain rho[g/cm^3] " << op_rlo << " .. " << op_rhi << "  (nr = "
              << ot.nr_ << ", dlogrho = " << ot.dlr_ << ")\n"
              << "###   domain T[K]        " << op_tlo << " .. " << op_thi << "  (nT = "
              << ot.nT_ << ", dlogT = " << ot.dlT_ << ")\n"
              << "###   lookups outside this box are CLAMPED to the edge value."
              << std::endl;
    if (eos_tab.loaded_) {
      // EosTable stores its rho axis in CODE units (eos_table.hpp:151 shifts a cgs file by
      // log10(rho_unit)); a legacy file carries the generator's own code units, which differ
      // by 0.14 % -- irrelevant against a domain check spanning decades. Its T axis is log10
      // K in both cases. Put both tables in cgs and compare.
      const Real eos_rlo = std::pow(10.0, eos_tab.lr0_) * U.rho_unit;
      const Real eos_rhi =
          std::pow(10.0, eos_tab.lr0_ + (eos_tab.nr_ - 1) * eos_tab.dlr_) * U.rho_unit;
      const Real eos_tlo = std::pow(10.0, eos_tab.lT0_);
      const Real eos_thi =
          std::pow(10.0, eos_tab.lT0_ + (eos_tab.nT_ - 1) * eos_tab.dlT_);
      const bool short_rho = op_rhi < eos_rhi * 0.999;
      const bool short_T = op_thi < eos_thi * 0.999;
      if (short_rho || short_T) {
        std::cout << "### WARNING Radiation (A3): the opacity table does NOT cover the "
                     "thermodynamic range the EOS table claims. The EOS stays valid where "
                     "the opacity is silently clamped to its edge value.\n";
        if (short_rho)
          std::cout << "###   rho: opacity ends at " << op_rhi << " g/cm^3, EOS reaches "
                    << eos_rhi << " g/cm^3 (short by "
                    << std::log10(eos_rhi / op_rhi) << " decades)\n";
        if (short_T)
          std::cout << "###   T:   opacity ends at " << op_thi << " K, EOS reaches "
                    << eos_thi << " K (short by " << std::log10(eos_thi / op_thi)
                    << " decades)\n";
        std::cout << "###   If this run can reach those conditions, regenerate the table:\n"
                     "###     gen_opacity_table.py --out <file> --nr 101 --nT 201 "
                     "--rho_min 1e-20 --rho_max 1e0 --T_min 3.0 --T_max 3.0e5\n"
                     "###   and record it in src/eos/TABLE_PROVENANCE.md."
                  << std::endl;
      } else {
        std::cout << "###   covers the EOS table domain (rho <= " << eos_rhi
                  << " g/cm^3, T <= " << eos_thi << " K): OK." << std::endl;
      }
    }
  }

  pkg->AddParam("tfloor", pin->GetOrAddReal(bn, "tfloor_code", 1.0e-3)); // T floor (code)
  pkg->AddParam("inner_iteration_max", pin->GetOrAddInteger(bn, "inner_iteration_max", 100));
  pkg->AddParam("inner_iteration_tol", pin->GetOrAddReal(bn, "inner_iteration_tol", 1.0e-8));
  // Upper T bracket for the multigroup coupling rtsafe solve (code T). Default = the opacity/EOS
  // table Tmax (stay in range); overridable. Ideal-EOS runs (no table) use a large fixed cap.
  pkg->AddParam("coupling_tmax",
                pin->GetOrAddReal(bn, "coupling_tmax_K", op_Tmax) / T_unit);

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
