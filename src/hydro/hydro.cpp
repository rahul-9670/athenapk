//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2020-2021, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Parthenon headers
#include <parthenon/package.hpp>

// AthenaPK headers
#include "../eos/adiabatic_glmmhd.hpp"
#include "../eos/adiabatic_hydro.hpp"
#include "../main.hpp"
#include "../pgen/pgen.hpp"
#include "../recon/dc_simple.hpp"
#include "../recon/limo3_simple.hpp"
#include "../recon/plm_simple.hpp"
#include "../recon/ppm_simple.hpp"
#include "../recon/weno3_simple.hpp"
#include "../recon/wenoz_simple.hpp"
#include "../refinement/refinement.hpp"
#include "../tracers/tracers.hpp"
#include "../sinks/sinks.hpp"
#include "../units.hpp"
#include "../units/physical_units.hpp"       // flagship Phase 1: one unit system
#include "../units/ionization_environment.hpp" // flagship Phase 1: shared CR rate
#include "defs.hpp"
#include "ct/ct.hpp"
#include "../diagnostics/angmom_diag.hpp"
#include "../diagnostics/cons_diag.hpp"
#include "../diagnostics/mag_diag.hpp"
#include "diffusion/diffusion.hpp"
#include "glmmhd/glmmhd.hpp"
#include "hydro.hpp"
#include "interface/params.hpp"
#include "outputs/outputs.hpp"
#include "prolongation/custom_ops.hpp"
#include "rsolvers/rsolvers.hpp"
#include "srcterms/tabular_cooling.hpp"
#include "utils/error_checking.hpp"
#include "../self_gravity/self_gravity.hpp"
#include "../radiation/radiation.hpp"
#include "../chemistry/chemistry.hpp"
#include "../dust/dust_pkg.hpp"

using namespace parthenon::package::prelude;

// *************************************************//
// define the "physics" package Hydro, which  *//
// includes defining various functions that control*//
// how parthenon functions and any tasks needed to *//
// implement the "physics"                         *//
// *************************************************//

namespace Hydro {

using cooling::TabularCooling;
using parthenon::HistoryOutputVar;

parthenon::Packages_t ProcessPackages(std::unique_ptr<ParameterInput> &pin) {
  parthenon::Packages_t packages;
  packages.Add(Hydro::Initialize(pin.get()));
  packages.Add(Tracers::Initialize(pin.get()));
  packages.Add(Sinks::Initialize(pin.get()));
  if (pin->GetOrAddBoolean("physics", "self_gravity", false)) {
    packages.Add(SelfGravity::Initialize(pin.get()));
  }
  if (pin->GetOrAddBoolean("physics", "radiation", false)) {
    packages.Add(Radiation::Initialize(pin.get()));
  }
  if (pin->GetOrAddBoolean("physics", "chemistry", false)) {
    packages.Add(Chemistry::Initialize(pin.get()));
  }
  if (pin->GetOrAddBoolean("physics", "dust", false)) {
    packages.Add(Dust::Initialize(pin.get()));
  }
  return packages;
}

// Calculate mininum dx, which is used in calculating the divergence cleaning speed c_h
// TODO(PG) eventually move to calculating the timestep once the timestep calc
// has been moved to be done before Step()
Real CalculateGlobalMinDx(MeshData<Real> *md) {
  auto *pmb = md->GetBlockData(0)->GetBlockPointer();
  auto hydro_pkg = pmb->packages.Get("Hydro");

  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real mindx = std::numeric_limits<Real>::max();

  bool nx2 = prim_pack.GetDim(2) > 1;
  bool nx3 = prim_pack.GetDim(3) > 1;
  pmb->par_reduce(
      "CalculateGlobalMinDx", 0, prim_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s,
      ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &lmindx) {
        const auto &coords = prim_pack.GetCoords(b);
        lmindx = fmin(lmindx, coords.Dxc<1>(k, j, i));
        if (nx2) {
          lmindx = fmin(lmindx, coords.Dxc<2>(k, j, i));
        }
        if (nx3) {
          lmindx = fmin(lmindx, coords.Dxc<3>(k, j, i));
        }
      },
      Kokkos::Min<Real>(mindx));

  return mindx;
}

// Using this per cycle function to populate various variables in
// Params that require global reduction *and* need to be set/known when
// the task list is constructed (versus when the task list is being executed).
// TODO(next person touching this function): If more/separate feature are required
// please separate concerns.
void PreStepMeshUserWorkInLoop(Mesh *pmesh, ParameterInput *pin, SimTime &tm) {
  auto hydro_pkg = pmesh->packages.Get("Hydro");

  // Calculate hyperbolic divergence cleaning speed
  // TODO(pgrete) Calculating mindx is only required after remeshing. Need to
  // find a clean solution for this one-off global reduction.
  if (hydro_pkg->Param<bool>("calc_c_h") ||
      hydro_pkg->Param<DiffInt>("diffint") != DiffInt::none) {

    Real mindx = std::numeric_limits<Real>::max();
    // Going over default partitions. Not using a (new) single partition containing
    // all blocks here as this (default) split is also used main Step() function and
    // thus does not create an overhead (such as creating a new MeshBlockPack that is just
    // used here). All partitions are executed sequentially. Given that a par_reduce to a
    // host var is blocking it's save to dirctly use the return value.
    const int num_partitions = pmesh->DefaultNumPartitions();
    for (int i = 0; i < num_partitions; i++) {
      auto &mu0 = pmesh->mesh_data.GetOrAdd("base", i);
      mindx = std::min(mindx, CalculateGlobalMinDx(mu0.get()));
    }
#ifdef MPI_PARALLEL
    Real mins[3];
    mins[0] = mindx;
    mins[1] = hydro_pkg->Param<Real>("dt_hyp");
    mins[2] = hydro_pkg->Param<Real>("dt_diff");
    PARTHENON_MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, mins, 3, MPI_PARTHENON_REAL, MPI_MIN,
                                      MPI_COMM_WORLD));

    hydro_pkg->UpdateParam("mindx", mins[0]);
    hydro_pkg->UpdateParam("dt_hyp", mins[1]);
    hydro_pkg->UpdateParam("dt_diff", mins[2]);
#else
    hydro_pkg->UpdateParam("mindx", mindx);
    // dt_hyp and dt_diff are already set directly in Params when they're calculated
#endif
    // Finally update c_h
    const auto &cfl_hyp = hydro_pkg->Param<Real>("cfl");
    const auto &dt_hyp = hydro_pkg->Param<Real>("dt_hyp");
    mindx = hydro_pkg->Param<Real>("mindx");
    hydro_pkg->UpdateParam("c_h", cfl_hyp * mindx / dt_hyp);
  }
}

template <Hst hst, int idx = -1>
Real HydroHst(MeshData<Real> *md) {
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");

  const auto &cons_pack = md->PackVariables(std::vector<std::string>{"cons"});
  const bool three_d = cons_pack.GetNdim() == 3;

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real sum = 0.0;

  // Sanity checks
  if ((hst == Hst::idx) && (idx < 0)) {
    PARTHENON_FAIL("Idx based hst output needs index >= 0");
  }
  Kokkos::parallel_reduce(
      "HydroHst",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExecSpace(), {0, kb.s, jb.s, ib.s},
          {cons_pack.GetDim(5), kb.e + 1, jb.e + 1, ib.e + 1},
          {1, 1, 1, ib.e + 1 - ib.s}),
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &lsum) {
        const auto &cons = cons_pack(b);
        const auto &coords = cons_pack.GetCoords(b);

        if (hst == Hst::idx) {
          lsum += cons(idx, k, j, i) * coords.CellVolume(k, j, i);
        } else if (hst == Hst::ekin) {
          lsum += 0.5 / cons(IDN, k, j, i) *
                  (SQR(cons(IM1, k, j, i)) + SQR(cons(IM2, k, j, i)) +
                   SQR(cons(IM3, k, j, i))) *
                  coords.CellVolume(k, j, i);
        } else if (hst == Hst::emag) {
          lsum += 0.5 *
                  (SQR(cons(IB1, k, j, i)) + SQR(cons(IB2, k, j, i)) +
                   SQR(cons(IB3, k, j, i))) *
                  coords.CellVolume(k, j, i);
          // relative divergence of B error, i.e., L * |div(B)| / |B|
        } else if (hst == Hst::divb) {
          Real divb =
              (cons(IB1, k, j, i + 1) - cons(IB1, k, j, i - 1)) / coords.Dxc<1>(k, j, i) +
              (cons(IB2, k, j + 1, i) - cons(IB2, k, j - 1, i)) / coords.Dxc<2>(k, j, i);
          if (three_d) {
            divb += (cons(IB3, k + 1, j, i) - cons(IB3, k - 1, j, i)) /
                    coords.Dxc<3>(k, j, i);
          }

          Real abs_b = std::sqrt(SQR(cons(IB1, k, j, i)) + SQR(cons(IB2, k, j, i)) +
                                 SQR(cons(IB3, k, j, i)));

          lsum += (abs_b != 0) ? 0.5 *
                                     (std::sqrt(SQR(coords.Dxc<1>(k, j, i)) +
                                                SQR(coords.Dxc<2>(k, j, i)) +
                                                SQR(coords.Dxc<3>(k, j, i)))) *
                                     std::abs(divb) / abs_b * coords.CellVolume(k, j, i)
                               : 0; // Add zero when abs_b ==0
        }
      },
      sum);

  return sum;
}

// WS-5c: peak (max over all cells) of the relative divergence error L*|div(B)|/|B|, for the
// div B fidelity audit. Registered with UserHistoryOperation::max so Parthenon takes the
// global max across blocks/ranks. Reports the WORST-cell cleaning error (the "relDivB" sum
// above is a volume integral -> a mean; near first-core formation the peak is what matters).
Real HydroHstMaxDivB(MeshData<Real> *md) {
  const auto &cons_pack = md->PackVariables(std::vector<std::string>{"cons"});
  const bool three_d = cons_pack.GetNdim() == 3;
  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real vmax = 0.0;
  Kokkos::parallel_reduce(
      "HydroHstMaxDivB",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExecSpace(), {0, kb.s, jb.s, ib.s},
          {cons_pack.GetDim(5), kb.e + 1, jb.e + 1, ib.e + 1},
          {1, 1, 1, ib.e + 1 - ib.s}),
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &lmax) {
        const auto &cons = cons_pack(b);
        const auto &coords = cons_pack.GetCoords(b);
        Real divb =
            (cons(IB1, k, j, i + 1) - cons(IB1, k, j, i - 1)) / coords.Dxc<1>(k, j, i) +
            (cons(IB2, k, j + 1, i) - cons(IB2, k, j - 1, i)) / coords.Dxc<2>(k, j, i);
        if (three_d) {
          divb += (cons(IB3, k + 1, j, i) - cons(IB3, k - 1, j, i)) / coords.Dxc<3>(k, j, i);
        }
        const Real abs_b = std::sqrt(SQR(cons(IB1, k, j, i)) + SQR(cons(IB2, k, j, i)) +
                                     SQR(cons(IB3, k, j, i)));
        const Real dx = 0.5 * std::sqrt(SQR(coords.Dxc<1>(k, j, i)) +
                                        SQR(coords.Dxc<2>(k, j, i)) +
                                        SQR(coords.Dxc<3>(k, j, i)));
        const Real rel = (abs_b != 0.0) ? dx * std::abs(divb) / abs_b : 0.0;
        if (rel > lmax) lmax = rel;
      },
      Kokkos::Max<Real>(vmax));
  return vmax;
}

// TOOD(pgrete) check is we can enlist this with FillDerived directly
// this is the package registered function to fill derived, here, convert the
// conserved variables to primitives
template <class T>
void ConsToPrim(MeshData<Real> *md) {
  const auto &eos =
      md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro")->Param<T>("eos");
  eos.ConservedToPrimitive(md);
}

// Add unsplit sources, i.e., source that are integrated in all stages of the
// explicit integration scheme.
// Note 1: Given that the sources are integrated in an unsplit manner, ensure
// that potential timestep constrains are also properly enforced when the
// respective source in active.
// Note 2: Directly update the "cons" variables based on the "prim" variables
// as the "cons" variables have already been updated when this function is called.
TaskStatus AddUnsplitSources(MeshData<Real> *md, const SimTime &tm, const Real beta_dt) {
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");

  // GLM/Dedner source terms. Under constrained transport these are not merely redundant,
  // they are ACTIVELY HARMFUL and are skipped (ct_glm_inert):
  //   * the extended (Powell) momentum source -(div B) * B uses a cell-centered 2-delta
  //     divergence of the *projected* face field. The CT invariant is the FACE divergence
  //     (measured 4e-12 here); the cell-centered one is pure truncation noise, O(1) at a
  //     shock. Its acceleration scales as 1/rho, so it runs away exactly in evacuated
  //     magnetized cells: measured dv = 0.22-0.39 per step at rho~1e-3 on the CT flagship,
  //     versus 4e-7 max anywhere on the matched GLM run.
  //   * the extended energy source -1/2 B.grad(psi) is driven by a psi that CT has left
  //     unbounded (see AdiabaticGLMMHDEOS::ConsToPrim): measured ~1e6 x the cell internal
  //     energy per step in those same cells.
  // Together they evacuate a pocket at the accretion-shock edge down to the density/pressure
  // floors, which is what produced the "CT over-magnetizes the core edge" (max ME/E -> 1)
  // artifact and the dt wall. |B| itself was never anomalous. See DEV_LOG 2026-07-29.
  if (hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd &&
      !hydro_pkg->Param<bool>("ct_glm_inert")) {
    hydro_pkg->Param<GLMMHD::SourceFun_t>("glmmhd_source")(md, beta_dt);
  }
  const auto &enable_cooling = hydro_pkg->Param<Cooling>("enable_cooling");

  if (enable_cooling == Cooling::tabular) {
    const TabularCooling &tabular_cooling =
        hydro_pkg->Param<TabularCooling>("tabular_cooling");

    tabular_cooling.SrcTerm(md, beta_dt);
  }
  if (ProblemSourceUnsplit != nullptr) {
    ProblemSourceUnsplit(md, tm, beta_dt);
  }

  return TaskStatus::complete;
}

TaskStatus AddSplitSourcesFirstOrder(MeshData<Real> *md, const SimTime &tm) {
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");

  if (ProblemSourceFirstOrder != nullptr) {
    ProblemSourceFirstOrder(md, tm, tm.dt);
  }
  return TaskStatus::complete;
}

TaskStatus AddSplitSourcesStrang(MeshData<Real> *md, const SimTime &tm) {
  if (ProblemSourceStrangSplit != nullptr) {
    ProblemSourceStrangSplit(md, tm, tm.dt);
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Build the reduced CR+thermal+grain ionization model (NICIL-class, MRN charged grains)
//! shared by the ambipolar / Ohmic / Hall non-ideal terms when their coefficient is
//! "ionization". Reads <diffusion> ion_* keys; all default to the shared FHC code-unit
//! calibration. The MRN grain bins are filled host-side here.
namespace {
Ionization::IonizationModel BuildIonizationModel(ParameterInput *pin) {
  Ionization::IonizationModel ion;
  // Flagship Phase 1: unit conversions come from the single authoritative unit system
  // (was: ion_rho_unit_cgs/ion_T_unit_K/ion_B_unit_G/ion_eta_unit hardcoded defaults here,
  // duplicating chemistry.cpp/radiation.cpp). rho_unit/T_unit bit-identical with chemistry;
  // B_unit now the exact sqrt(4*pi*rho)*v (audit #3); eta_unit = time/length^2.
  const auto U = PhysUnits::BuildPhysicalUnits(pin);
  ion.rho_unit = U.rho_unit;
  ion.T_unit = U.temperature_unit;
  ion.B_unit = U.magnetic_unit;
  ion.eta_unit = U.diffusivity_unit;
  // Cosmic-ray rate + neutral composition from the shared IonizationEnvironment: the same
  // object chemistry consumes, so Ohm/Hall/AD and the chemistry x_e describe one charge
  // population (audit fix #1, now a shared object rather than an init-time assertion).
  const auto env = PhysUnits::BuildIonizationEnvironment(pin);
  ion.zeta = env.zeta_cr_cgs;
  ion.mu_n = env.mu_n;
  ion.m_ion = pin->GetOrAddReal("diffusion", "ion_m_ion", ion.m_ion);
  ion.alpha0 = pin->GetOrAddReal("diffusion", "ion_alpha0", ion.alpha0);
  ion.sigma_en = pin->GetOrAddReal("diffusion", "ion_sigma_en", ion.sigma_en);
  ion.sigv_in = pin->GetOrAddReal("diffusion", "ion_sigv_in", ion.sigv_in);
  ion.gamma_AD = pin->GetOrAddReal("diffusion", "ion_gamma_AD", ion.gamma_AD);
  // grains (MRN distribution)
  ion.f_dg = pin->GetOrAddReal("diffusion", "ion_f_dg", ion.f_dg);
  // Phase-3 fix (2026-07-25): the pre-fix charge-state solvers were both defective -- the
  // relaxed fixed point in SolveCharges does not converge once grains dominate the charge
  // budget (rho >~ 1e-12), and SahaThermal's fixed-step absolute bisection injects spurious
  // electrons in cold gas. Default is now the robust pair. RESULT-CHANGING where grains
  // matter -- set true only to bit-reproduce pre-fix runs.
  ion.legacy_charge_solver =
      pin->GetOrAddBoolean("diffusion", "ion_legacy_charge_solver", false);
  ion.rho_grain = pin->GetOrAddReal("diffusion", "ion_rho_grain", ion.rho_grain);
  ion.T_subl = pin->GetOrAddReal("diffusion", "ion_T_subl", ion.T_subl);
  const Real a_min = pin->GetOrAddReal("diffusion", "ion_a_min_cm", 5.0e-7);
  const Real a_max = pin->GetOrAddReal("diffusion", "ion_a_max_cm", 2.5e-5);
  const Real mrn_p = pin->GetOrAddReal("diffusion", "ion_mrn_p", 3.5);
  Ionization::SetupGrainBins(ion, a_min, a_max, mrn_p);
  // thermal Saha
  ion.x_K = pin->GetOrAddReal("diffusion", "ion_x_K", ion.x_K);
  ion.chi_K = pin->GetOrAddReal("diffusion", "ion_chi_K", ion.chi_K);
  // ambipolar closure: single_fluid (Athena++-matched, default) or tensor (grain-modified)
  const auto cl = pin->GetOrAddString("diffusion", "ion_ad_closure", "single_fluid");
  if (cl == "tensor") {
    ion.ad_closure = Ionization::ADClosure::tensor;
  } else if (cl == "single_fluid") {
    ion.ad_closure = Ionization::ADClosure::single_fluid;
  } else {
    PARTHENON_FAIL("diffusion/ion_ad_closure must be 'single_fluid' or 'tensor'.");
  }
  return ion;
}
} // namespace

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {
  auto pkg = std::make_shared<StateDescriptor>("Hydro");

  Real cfl = pin->GetOrAddReal("parthenon/time", "cfl", 0.3);
  pkg->AddParam<>("cfl", cfl);

  bool pack_in_one = pin->GetOrAddBoolean("parthenon/mesh", "pack_in_one", true);
  pkg->AddParam<>("pack_in_one", pack_in_one);

  // Flagship Phase 1 provenance: stamp the authoritative unit system + shared ionization
  // environment into the Hydro package as Restart-flagged params so EVERY restart/output
  // records the exact code->cgs mapping and cosmic-ray/composition inputs the run used --
  // no more reconstructing units from the deck. Numerically inert (metadata only; no kernel
  // reads these). Reading an older restart that predates these attrs is safe: Parthenon's
  // ReadFromHDF5AllParamsOfType catches the missing-attribute error and keeps the value
  // recomputed here (params.cpp). These are the SAME builders the physics packages consume.
  {
    const auto Uprov = PhysUnits::BuildPhysicalUnits(pin);
    const auto envprov = PhysUnits::BuildIonizationEnvironment(pin);
    using RM = Params::Mutability;
    pkg->AddParam<Real>("phys_units/rho_unit_cgs", Uprov.rho_unit, RM::Restart);
    pkg->AddParam<Real>("phys_units/v_unit_cgs", Uprov.v_unit, RM::Restart);
    pkg->AddParam<Real>("phys_units/length_unit_cgs", Uprov.length_unit, RM::Restart);
    pkg->AddParam<Real>("phys_units/time_unit_cgs", Uprov.time_unit, RM::Restart);
    pkg->AddParam<Real>("phys_units/mass_unit_cgs", Uprov.mass_unit, RM::Restart);
    pkg->AddParam<Real>("phys_units/energy_density_unit_cgs", Uprov.energy_density_unit,
                        RM::Restart);
    pkg->AddParam<Real>("phys_units/magnetic_unit_G", Uprov.magnetic_unit, RM::Restart);
    pkg->AddParam<Real>("phys_units/temperature_unit_K", Uprov.temperature_unit, RM::Restart);
    pkg->AddParam<Real>("phys_units/diffusivity_unit", Uprov.diffusivity_unit, RM::Restart);
    pkg->AddParam<Real>("phys_units/opacity_unit", Uprov.opacity_unit, RM::Restart);
    pkg->AddParam<Real>("phys_units/mu_thermal", Uprov.mu_thermal, RM::Restart);
    pkg->AddParam<Real>("phys_units/zeta_cr_cgs", envprov.zeta_cr_cgs, RM::Restart);
    pkg->AddParam<Real>("phys_units/mu_n", envprov.mu_n, RM::Restart);
  }

  const auto fluid_str = pin->GetOrAddString("hydro", "fluid", "euler");
  auto fluid = Fluid::undefined;
  bool calc_c_h = false; // calculate hyperbolic divergence cleaning speed
  int nhydro = -1;

  if (fluid_str == "euler") {
    fluid = Fluid::euler;
    nhydro = GetNVars<Fluid::euler>();
  } else if (fluid_str == "glmmhd") {
    fluid = Fluid::glmmhd;
    nhydro = GetNVars<Fluid::glmmhd>();
    // TODO(pgrete) reeval default value based on testing
    auto glmmhd_source_str =
        pin->GetOrAddString("hydro", "glmmhd_source", "dedner_plain");
    if (glmmhd_source_str == "dedner_plain") {
      pkg->AddParam<GLMMHD::SourceFun_t>("glmmhd_source", GLMMHD::DednerSource<false>);
    } else if (glmmhd_source_str == "dedner_extended") {
      pkg->AddParam<GLMMHD::SourceFun_t>("glmmhd_source", GLMMHD::DednerSource<true>);
    } else {
      PARTHENON_FAIL("AthenaPK hydro: Unknown glmmhd_source");
    }
    // ratio of diffusive to advective timescale of the divergence cleaning
    auto glmmhd_alpha = pin->GetOrAddReal("hydro", "glmmhd_alpha", 0.1);
    pkg->AddParam<Real>("glmmhd_alpha", glmmhd_alpha);
    calc_c_h = true;
  } else {
    PARTHENON_FAIL("AthenaPK hydro: Unknown fluid method.");
  }
  pkg->AddParam<>("fluid", fluid);
  pkg->AddParam<>("nhydro", nhydro);
  pkg->AddParam<>("calc_c_h", calc_c_h);

  // Divergence control (Phase 2). "glm" = current GLM/Dedner path (default, bit-identical
  // to the previous behavior). "ct" = staggered Constrained Transport (increment 1: ideal
  // MHD, single level). On the CT path a face-centered Bf field is the primary magnetic
  // variable and IB1..IB3 become its projection; the GLM machinery (psi flux + Dedner
  // damping) is left running but decoupled -- the induction update comes from the CT curl
  // and the projection overwrites the cell-centered B each substage. See docs/CT_DESIGN.md.
  const auto divergence_control_str =
      pin->GetOrAddString("hydro", "divergence_control", "glm");
  bool use_ct = false;
  if (divergence_control_str == "glm") {
    use_ct = false;
  } else if (divergence_control_str == "ct") {
    PARTHENON_REQUIRE(fluid == Fluid::glmmhd,
                      "hydro/divergence_control=ct requires fluid=glmmhd.");
    use_ct = true;
  } else {
    PARTHENON_FAIL("AthenaPK hydro: Unknown divergence_control (use 'glm' or 'ct').");
  }
  pkg->AddParam<>("use_ct", use_ct);
  // CT edge-EMF averaging scheme (increment 3). "gs05" = Gardiner & Stone (2005) upwind
  // CT (default; matches Athena++, adds Riemann dissipation via the transverse-B face
  // fluxes); "arithmetic" = Balsara & Spicer (1999) plain average (increment-1 fallback).
  // Both are divergence-free by construction; this only sets the EMF *value*.
  const auto ct_emf_str = pin->GetOrAddString("hydro", "ct_emf", "gs05");
  bool use_ct_gs05 = true;
  if (ct_emf_str == "gs05") {
    use_ct_gs05 = true;
  } else if (ct_emf_str == "arithmetic") {
    use_ct_gs05 = false;
  } else {
    PARTHENON_FAIL("AthenaPK hydro: Unknown ct_emf (use 'gs05' or 'arithmetic').");
  }
  pkg->AddParam<>("use_ct_gs05", use_ct_gs05);

  // Under CT the whole GLM/Dedner apparatus (psi transport + Powell/extended source terms)
  // must be INERT: CT already guarantees div B = 0 in the face sense, and because the
  // cell-centered B is overwritten by the face projection each substage the psi -> B half of
  // the cleaning loop is discarded, leaving psi driven but never relieved. Setting
  // ct_legacy_glm_source=true restores the old (broken) behaviour for A/B bisects only.
  // One-sided internal-energy guard on the CT magnetic-energy replacement, see
  // CT_ProjectBfToCC. Fraction of the pre-projection internal energy below which the
  // projection is not allowed to push the gas. 0 disables (bit-identical, legacy).
  //
  // Calibrated on runs/ct_tests/orszag_tang_ad_{ct,glm}.in (CT+AD+RKL2, 128^2, tlim=0.14):
  //   guard   outcome                          tot-E vs matched GLM
  //   0.0     ABORT, negative pressure @cyc~100        --
  //   0.5     ABORT, guard never fires                 --
  //   0.9     completes                             +0.29%   <- default
  //   0.99    completes                             +2.76%
  // At 0.9 the completed run matches the GLM reference to 0.1% in KE and 1.0% in ME, and
  // ct_maxAbsDivB stays at 3.5e-14. Note the guard is INERT outside the RKL2 path: ideal-CT
  // and CT+AD-unsplit are bit-identical at 0.0 and 0.9. It is a LIMITER, not a cure -- the
  // underlying defect is that cons(IEN) and the CT field are advanced by two independent
  // integrators (see CT_ProjectBfToCC); a consistent-heating / dual-energy formulation is
  // the real fix and is not implemented here.
  const Real ct_eint_guard_frac =
      use_ct ? pin->GetOrAddReal("hydro", "ct_eint_guard_frac", 0.9) : 0.0;
  PARTHENON_REQUIRE(ct_eint_guard_frac >= 0.0 && ct_eint_guard_frac < 1.0,
                    "hydro/ct_eint_guard_frac must be in [0,1).");
  pkg->AddParam<>("ct_eint_guard_frac", ct_eint_guard_frac);

  // CT + RKL2: build the diffusive Poynting energy flux from the same edge EMF that drives
  // the CT induction (CT_AddDiffusivePoynting) instead of the face-based deposits in
  // resistivity.cpp/ambipolar.cpp. Those use a different stencil, and the mismatch lands in
  // the recovered internal energy e = E - KE - ME. Isolated experimentally: with the eint
  // guard OFF, CT+RKL2 Orszag-Tang is clean with Ohmic alone (matched stencils) but aborts on
  // negative pressure with ambipolar (direct-edge EMF vs face-based Poynting). Unsplit CT is
  // unaffected and stays on the face deposits (Bf.flux there carries ideal+diffusive).
  // Set false to restore the old behaviour for A/B bisects.
  // RKL2/STS ONLY: there CT_ZeroEMF -> Ohmic -> ambipolar leaves Bf.flux holding the
  // DIFFUSIVE-ONLY edge EMF, which is what the replacement flux needs. In the unsplit path
  // Bf.flux carries ideal+diffusive and the ideal Poynting is already in the HLLD energy
  // flux, so it must keep the face-based deposits (and it is not affected by the defect --
  // CT+AD-unsplit runs clean). Same default as the real parse site below, so this is a
  // read-only peek at diffusion/integrator.
  const bool ct_diffint_rkl2 =
      pin->GetOrAddString("diffusion", "integrator", "none") == "rkl2";
  const bool ct_edge_poynting =
      use_ct && ct_diffint_rkl2 && pin->GetOrAddBoolean("hydro", "ct_edge_poynting", true);
  pkg->AddParam<>("ct_edge_poynting", ct_edge_poynting);

  // ENERGY-NEUTRAL PROJECTION -- **DEFAULT OFF; EXPERIMENTAL, DO NOT ENABLE FOR SCIENCE.**
  // Carries the projection's magnetic-energy change into the total energy so that
  // e = E - KE - ME is held exactly fixed across the face->cell projection. Motivation: the
  // projection systematically HEATS magnetized low-density cells (+1.5e-3 relative per
  // application at ME/IE>1 in the flagship), and the damage is amplified by ME/IE (~150
  // there), which is the mechanism behind the evacuated-hole defect (DEV_LOG 2026-07-29).
  //
  // WHY IT IS OFF BY DEFAULT: it trades an amplified LOCAL error for a GLOBAL conservation
  // error, and the global cost is unacceptable. Measured on orszag_tang_ad_ct to t=0.20
  // (identical ICs, tot-E = 3.492570e-01 at t=0): ON drifts to 3.318490e-01, **-4.98% total
  // energy**, while OFF and the GLM reference both conserve to +0.0000%. dt also falls 2.7x
  // (4.97e-4 vs GLM 1.32e-3). An a-priori estimate from the flagship dump suggested ~1% per
  // 6.6e6 cycles; that was WRONG because it assumed 2 projections per cycle (VL2). Under RKL2
  // the projection runs at EVERY super-time-step substage -- 10-26 per cycle here -- so the
  // real rate is ~1e4 times higher. Lesson: count the actual call sites, not the stage count.
  //
  // The underlying inconsistency is that cons(IB) (advanced by the HLLD flux divergence) and
  // the projection of the CT face field are two independent estimates of the same quantity
  // whose difference is systematically signed, not zero-mean. Neither dumping it in e (default)
  // nor in E (this option) is correct; the principled fix is a dual-energy formulation that
  // evolves the internal energy directly. Kept as an opt-in strictly for A/B experiments.
  const bool ct_energy_neutral =
      use_ct && pin->GetOrAddBoolean("hydro", "ct_energy_neutral_projection", false);
  pkg->AddParam<>("ct_energy_neutral_projection", ct_energy_neutral);

  // Opt-in per-cell diagnostic of the projection's internal-energy transfer. See
  // CT_ProjectBfToCC. Allocates one extra derived cell field; default off.
  const bool ct_proj_diag =
      use_ct && pin->GetOrAddBoolean("hydro", "ct_proj_diag", false);
  pkg->AddParam<>("ct_proj_diag", ct_proj_diag);

  const bool ct_legacy_glm_source =
      pin->GetOrAddBoolean("hydro", "ct_legacy_glm_source", false);
  const bool ct_glm_inert = use_ct && !ct_legacy_glm_source;
  pkg->AddParam<>("ct_glm_inert", ct_glm_inert);
  if (ct_glm_inert && parthenon::Globals::my_rank == 0) {
    std::cout << "## CT: GLM/Dedner machinery held inert (psi == 0, no Powell/extended "
                 "source terms)."
              << std::endl;
  }
  // Following params should (currently) be present independent of solver because
  // they're all used in the main loop.
  // TODO(pgrete) think about which approach (selective versus always is preferable)
  pkg->AddParam<Real>(
      "c_h", 0.0, Params::Mutability::Mutable); // hyperbolic divergence cleaning speed
  // global minimum dx (used to calc c_h)
  pkg->AddParam<Real>("mindx", std::numeric_limits<Real>::max(),
                      Params::Mutability::Mutable);
  // hyperbolic timestep constraint
  pkg->AddParam<Real>("dt_hyp", std::numeric_limits<Real>::max(),
                      Params::Mutability::Mutable);

  const auto recon_str = pin->GetString("hydro", "reconstruction");
  int recon_need_nghost = 3; // largest number for the choices below
  auto recon = Reconstruction::undefined;
  if (recon_str == "dc") {
    recon = Reconstruction::dc;
    recon_need_nghost = 1;
  } else if (recon_str == "plm") {
    recon = Reconstruction::plm;
    recon_need_nghost = 2;
  } else if (recon_str == "ppm") {
    recon = Reconstruction::ppm;
    recon_need_nghost = 3;
  } else if (recon_str == "limo3") {
    recon = Reconstruction::limo3;
    recon_need_nghost = 2;
  } else if (recon_str == "weno3") {
    recon = Reconstruction::weno3;
    recon_need_nghost = 2;
  } else if (recon_str == "wenoz") {
    recon = Reconstruction::wenoz;
    recon_need_nghost = 3;
  } else {
    PARTHENON_FAIL("AthenaPK hydro: Unknown reconstruction method.");
  }
  // Adding recon independently of flux function pointer as it's used in 3D flux func.
  pkg->AddParam<>("reconstruction", recon);

  // Use hyperbolic timestep constraint by default
  bool calc_dt_hyp = true;
  const auto riemann_str = pin->GetString("hydro", "riemann");
  auto riemann = RiemannSolver::undefined;
  if (riemann_str == "llf") {
    riemann = RiemannSolver::llf;
    PARTHENON_REQUIRE(recon == Reconstruction::dc,
                      "LLF Riemann solver only implemented with DC reconstruction.")
  } else if (riemann_str == "hlle") {
    riemann = RiemannSolver::hlle;
  } else if (riemann_str == "hllc") {
    riemann = RiemannSolver::hllc;
  } else if (riemann_str == "hlld") {
    riemann = RiemannSolver::hlld;
  } else if (riemann_str == "none") {
    riemann = RiemannSolver::none;
    // If hyperbolic fluxes are disabled, there's no restriction from those
    // on the timestep
    calc_dt_hyp = false;
    PARTHENON_REQUIRE(recon == Reconstruction::dc,
                      "Disabling hyperbolic fluxes via 'none' Riemann solver only "
                      "supported in comination with DC reconstruction.")
  } else {
    PARTHENON_FAIL("AthenaPK hydro: Unknown riemann solver.");
  }
  pkg->AddParam<>("riemann", riemann);

  // Set calculation of hyperbolic timestep. Input file option takes precedence.
  if (pin->DoesParameterExist("hydro", "calc_dt_hyp")) {
    calc_dt_hyp = pin->GetBoolean("hydro", "calc_dt_hyp");
  }
  pkg->AddParam<>("calc_dt_hyp", calc_dt_hyp);

  // Maximum dt. Useful for debugging.
  const auto max_dt = pin->GetOrAddReal("hydro", "max_dt", -1.0);
  pkg->AddParam<>("max_dt", max_dt);

  // Map contaning all compiled in flux functions
  std::map<std::tuple<Fluid, Reconstruction, RiemannSolver>, FluxFun_t *>
      flux_functions{};
  // TODO(?) The following line could potentially be set by configure-time options
  // so that the resulting binary can only contain a subset of included flux functions
  // to reduce size.
  add_flux_fun<Fluid::euler, Reconstruction::dc, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::dc, RiemannSolver::none>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::plm, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::ppm, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::weno3, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::limo3, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::wenoz, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::dc, RiemannSolver::hllc>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::plm, RiemannSolver::hllc>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::ppm, RiemannSolver::hllc>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::weno3, RiemannSolver::hllc>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::limo3, RiemannSolver::hllc>(flux_functions);
  add_flux_fun<Fluid::euler, Reconstruction::wenoz, RiemannSolver::hllc>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::dc, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::dc, RiemannSolver::none>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::plm, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::ppm, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::weno3, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::limo3, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::wenoz, RiemannSolver::hlle>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::dc, RiemannSolver::hlld>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::plm, RiemannSolver::hlld>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::ppm, RiemannSolver::hlld>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::weno3, RiemannSolver::hlld>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::limo3, RiemannSolver::hlld>(flux_functions);
  add_flux_fun<Fluid::glmmhd, Reconstruction::wenoz, RiemannSolver::hlld>(flux_functions);
  // Add first order recon with LLF fluxes (implemented for testing as tight loop)
  flux_functions[std::make_tuple(Fluid::euler, Reconstruction::dc, RiemannSolver::llf)] =
      Hydro::CalculateFluxesTight<Fluid::euler>;
  flux_functions[std::make_tuple(Fluid::glmmhd, Reconstruction::dc, RiemannSolver::llf)] =
      Hydro::CalculateFluxesTight<Fluid::glmmhd>;

  // flux used in all stages expect the first. First stage is set below based on integr.
  FluxFun_t *flux_other_stage = nullptr;
  flux_other_stage = flux_functions.at(std::make_tuple(fluid, recon, riemann));

  parthenon::HstVar_list hst_vars = {};
  hst_vars.emplace_back(HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                         HydroHst<Hst::idx, IDN>, "mass"));
  hst_vars.emplace_back(HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                         HydroHst<Hst::idx, IM1>, "1-mom"));
  hst_vars.emplace_back(HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                         HydroHst<Hst::idx, IM2>, "2-mom"));
  hst_vars.emplace_back(HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                         HydroHst<Hst::idx, IM3>, "3-mom"));
  hst_vars.emplace_back(
      HistoryOutputVar(parthenon::UserHistoryOperation::sum, HydroHst<Hst::ekin>, "KE"));
  hst_vars.emplace_back(HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                         HydroHst<Hst::idx, IEN>, "tot-E"));
  if (fluid == Fluid::glmmhd) {
    hst_vars.emplace_back(HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                           HydroHst<Hst::emag>, "ME"));
    hst_vars.emplace_back(HistoryOutputVar(parthenon::UserHistoryOperation::sum,
                                           HydroHst<Hst::divb>, "relDivB"));
    // WS-5c: peak relative div B error (max over cells), for the div B fidelity audit.
    hst_vars.emplace_back(HistoryOutputVar(parthenon::UserHistoryOperation::max,
                                           HydroHstMaxDivB, "maxRelDivB"));
  }
  pkg->AddParam<>(parthenon::hist_param_key, hst_vars, true);

  // not using GetOrAdd here until there's a reasonable default
  const auto nghost = pin->GetInteger("parthenon/mesh", "nghost");
  if (nghost < recon_need_nghost) {
    PARTHENON_FAIL("AthenaPK hydro: Need more ghost zones for chosen reconstruction.");
  }

  const auto integrator_str = pin->GetString("parthenon/time", "integrator");
  auto integrator = Integrator::undefined;
  FluxFun_t *flux_first_stage = flux_other_stage;

  if (integrator_str == "rk1") {
    integrator = Integrator::rk1;
  } else if (integrator_str == "rk2") {
    integrator = Integrator::rk2;
  } else if (integrator_str == "rk3") {
    integrator = Integrator::rk3;
  } else if (integrator_str == "vl2") {
    integrator = Integrator::vl2;
    // override first stage (predictor) to first order
    flux_first_stage =
        flux_functions.at(std::make_tuple(fluid, Reconstruction::dc, riemann));
  }
  pkg->AddParam<>("integrator", integrator);
  pkg->AddParam<FluxFun_t *>("flux_first_stage", flux_first_stage);
  pkg->AddParam<FluxFun_t *>("flux_other_stage", flux_other_stage);

  auto first_order_flux_correct =
      pin->GetOrAddBoolean("hydro", "first_order_flux_correct", false);
  pkg->AddParam<>("first_order_flux_correct", first_order_flux_correct);
  if (first_order_flux_correct) {
    if (fluid == Fluid::euler) {
      pkg->AddParam<FirstOrderFluxCorrectFun_t *>("first_order_flux_correct_fun",
                                                  FirstOrderFluxCorrect<Fluid::euler>);
    } else if (fluid == Fluid::glmmhd) {
      pkg->AddParam<FirstOrderFluxCorrectFun_t *>("first_order_flux_correct_fun",
                                                  FirstOrderFluxCorrect<Fluid::glmmhd>);
    }
  }

  if (pin->DoesBlockExist("units")) {
    Units units(pin, pkg);
  }

  auto eos_str = pin->GetString("hydro", "eos");
  if (eos_str == "adiabatic" || eos_str == "hydrogen") {
    Real gamma = pin->GetReal("hydro", "gamma");
    pkg->AddParam<>("AdiabaticIndex", gamma);

    // Optional tabulated protostellar EOS (second-core physics). eos=hydrogen loads the
    // offline full-Saha table (src/eos/gen_eos_table.py: H2 dissociation + H ionization +
    // He/He+ ionization + H2 rot/vib + inert He, in FHC code units) and interpolates it for
    // pressure/sound-speed/internal-energy. Requires riemann=hlld.
    const bool use_h2diss = (eos_str == "hydrogen");
    EOSTable::EosTable eos_tab;
    if (use_h2diss) {
      // WP-11 (2026-08-02): the fallback default below is the COARSE 180x220x200 table, and
      // silently taking it is a real accuracy hazard -- measured against eos_table_hires.bin it
      // carries 8.3% p99 and 37.6% worst-case error in the SOUND SPEED (vs <=0.51% p99 for the
      // hi-res table). Production decks all set eos_table_file explicitly to eos_table_hires.bin
      // (fhc_flagship.in, fhc_rootladder.in), so a missing key means a deck typo or an
      // unmaintained deck, not an intentional choice. Warn loudly rather than degrade quietly.
      // Not a hard failure: the coarse table is legitimate for smoke tests and cheap gates.
      const bool tabfile_given = pin->DoesParameterExist("hydro", "eos_table_file");
      const auto tabfile = pin->GetOrAddString(
          "hydro", "eos_table_file",
          "/beegfs/u/bbg6470/athenapk/src/eos/eos_table.bin");
      if (!tabfile_given) {
        PARTHENON_WARN(
            "<hydro> eos_table_file was not set, so eos=hydrogen is falling back to the COARSE "
            "built-in table (src/eos/eos_table.bin, 180x220x200). Measured interpolation error "
            "vs eos_table_hires.bin: 8.3% p99 / 37.6% max in cs^2. This is FINE for smoke tests "
            "and NOT fine for production -- set eos_table_file=.../src/eos/eos_table_hires.bin.");
      }
      eos_tab.Load(tabfile);
      const auto riem = pin->GetOrAddString("hydro", "riemann", "hlld");
      PARTHENON_REQUIRE(riem == "hlld",
                        "Tabulated protostellar EOS (eos=hydrogen) requires riemann=hlld "
                        "(HLLE's Roe average bakes in a constant gamma).");
      // The table EOS owns the gas thermodynamics; without the M1 RT package the pgen's
      // barotropic source overwrites e_th with an ideal-gamma law inconsistent with the
      // table. Require radiation on (RT then owns the energy; barotropic is skipped).
      PARTHENON_REQUIRE(pin->GetOrAddBoolean("physics", "radiation", false),
                        "eos=hydrogen requires <physics> radiation=true (else the "
                        "barotropic cooling overwrites the table EOS thermal energy).");
      if (parthenon::Globals::my_rank == 0) {
        std::cout << "## EOS: tabulated protostellar EOS ACTIVE (second-core). table="
                  << tabfile << " grid " << eos_tab.nr_ << "x" << eos_tab.ne_
                  << " (rho,esp), " << eos_tab.nr_ << "x" << eos_tab.nT_ << " (rho,T)"
                  << std::endl;
      }
    }

    if (pin->DoesParameterExist("hydro", "He_mass_fraction") &&
        pkg->AllParams().hasKey("units")) {
      auto units = pkg->Param<Units>("units");
      const auto He_mass_fraction = pin->GetReal("hydro", "He_mass_fraction");
      const auto mu = 1 / (He_mass_fraction * 3. / 4. + (1 - He_mass_fraction) * 2);
      const auto mu_e = 1 / (He_mass_fraction * 2. / 4. + (1 - He_mass_fraction));
      pkg->AddParam<>("mu", mu);
      pkg->AddParam<>("mu_e", mu_e);
      pkg->AddParam<>("He_mass_fraction", He_mass_fraction);
      pkg->AddParam<>("mbar", mu * units.atomic_mass_unit());
      // Following convention in the astro community, we're using mh as unit for the mean
      // molecular weight
      pkg->AddParam<>("mbar_over_kb", mu * units.mh() / units.k_boltzmann());
    }

    // By default disable floors by setting a negative value
    Real dfloor = pin->GetOrAddReal("hydro", "dfloor", -1.0);
    // WP-13: expose the density floor as a package Param. SelfGravity::FillPoissonRHS needs
    // it to reproduce prim's floored density exactly while reading `cons` (it must read cons,
    // not prim, because prim is another package's FillDerived output and the ordering between
    // packages is not guaranteed -- on restart prim was still zero when the Poisson RHS was
    // assembled). Registered unconditionally so it does not depend on any diagnostic gate.
    pkg->AddParam<Real>("grav_rho_floor", dfloor);
    Real pfloor = pin->GetOrAddReal("hydro", "pfloor", -1.0);
    Real Tfloor = pin->GetOrAddReal("hydro", "Tfloor", -1.0);
    Real efloor = Tfloor;
    if (efloor > 0.0) {
      if (!pkg->AllParams().hasKey("mbar_over_kb")) {
        PARTHENON_FAIL("Temperature floor requires units and gas composition. "
                       "Either set a 'units' block and the 'hydro/He_mass_fraction' in "
                       "input file or use a pressure floor "
                       "(defined code units) instead.");
      }
      auto mbar_over_kb = pkg->Param<Real>("mbar_over_kb");
      efloor = Tfloor / mbar_over_kb / (gamma - 1.0);
    }

    // By default disable ceilings by setting to infinity
    Real vceil =
        pin->GetOrAddReal("hydro", "vceil", std::numeric_limits<Real>::infinity());
    Real Tceil =
        pin->GetOrAddReal("hydro", "Tceil", std::numeric_limits<Real>::infinity());
    Real eceil = Tceil;
    if (eceil < std::numeric_limits<Real>::infinity()) {
      if (!pkg->AllParams().hasKey("mbar_over_kb")) {
        PARTHENON_FAIL("Temperature ceiling requires units and gas composition. "
                       "Either set a 'units' block and the 'hydro/He_mass_fraction' in "
                       "input file or use a pressure floor "
                       "(defined code units) instead.");
      }
      auto mbar_over_kb = pkg->Param<Real>("mbar_over_kb");
      eceil = Tceil / mbar_over_kb / (gamma - 1.0);
    }

    auto conduction = Conduction::none;
    auto conduction_str = pin->GetOrAddString("diffusion", "conduction", "none");
    if (conduction_str == "isotropic") {
      conduction = Conduction::isotropic;
    } else if (conduction_str == "anisotropic") {
      conduction = Conduction::anisotropic;
    } else if (conduction_str != "none") {
      PARTHENON_FAIL(
          "Unknown conduction method. Options are: none, isotropic, anisotropic");
    }
    // If conduction is enabled, process supported coefficients
    if (conduction != Conduction::none) {
      PARTHENON_REQUIRE_THROWS(
          typeid(parthenon::Coordinates_t) == typeid(parthenon::UniformCartesian),
          "Probably need to update derivative calc (with respect to spacing between the "
          "faces when calculating fluxes from cell centered quantities average to cell "
          "faces) in thermal conduction to support non-Cartesian coordiates.");
      auto conduction_coeff_str =
          pin->GetOrAddString("diffusion", "conduction_coeff", "none");
      auto conduction_coeff = ConductionCoeff::none;

      // Saturated conduction factor to account for "uncertainty", see
      // Cowie & McKee 77 and a value of 0.3 is typical chosen (though using "weak
      // evidence", see Balbus & MacKee 1982 and Max, McKee, and Mead 1980).
      const auto conduction_sat_phi =
          pin->GetOrAddReal("diffusion", "conduction_sat_phi", 0.3);
      Real conduction_sat_prefac = 0.0;

      if (conduction_coeff_str == "spitzer") {
        if (!pkg->AllParams().hasKey("mbar")) {
          PARTHENON_FAIL("Spitzer thermal conduction requires units and gas composition. "
                         "Please set a 'units' block and the 'hydro/He_mass_fraction' in "
                         "the input file.");
        }
        conduction_coeff = ConductionCoeff::spitzer;

        // Default value assume fully ionized hydrogen plasma with Coulomb logarithm of 40
        // to approximate ICM conditions, i.e., 1.84e-5/ln Lambda = 4.6e-7.
        Real spitzer_coeff =
            pin->GetOrAddReal("diffusion", "spitzer_cond_in_erg_by_s_K_cm", 4.6e-7);
        // Convert to code units. No temp conversion as [T_phys] = [T_code].
        auto units = pkg->Param<Units>("units");
        spitzer_coeff *= units.erg() / (units.s() * units.cm());

        const auto mbar = pkg->Param<Real>("mbar");
        auto thermal_diff =
            ThermalDiffusivity(conduction, conduction_coeff, spitzer_coeff, mbar,
                               units.electron_mass(), units.k_boltzmann());
        pkg->AddParam<>("thermal_diff", thermal_diff);

        const auto mu = pkg->Param<Real>("mu");
        // 6.86 again assumes a fully ionized hydrogen plasma in agreement with
        // the assumptions above (technically this means mu = 0.5) and can be derived
        // from eq (7) in CM77 assuming T_e = T_i.
        conduction_sat_prefac = 6.86 * std::sqrt(mu) * conduction_sat_phi;

      } else if (conduction_coeff_str == "fixed") {
        conduction_coeff = ConductionCoeff::fixed;
        Real thermal_diff_coeff_code =
            pin->GetReal("diffusion", "thermal_diff_coeff_code");
        auto thermal_diff = ThermalDiffusivity(conduction, conduction_coeff,
                                               thermal_diff_coeff_code, 0.0, 0.0, 0.0);
        pkg->AddParam<>("thermal_diff", thermal_diff);
        // 5.0 prefactor comes from eq (8) in Cowie & McKee 1977
        // https://doi.org/10.1086/154911
        conduction_sat_prefac = 5.0 * conduction_sat_phi;

      } else {
        PARTHENON_FAIL("Thermal conduction is enabled but no coefficient is set. Please "
                       "set diffusion/conduction_coeff to either 'spitzer' or 'fixed'");
      }
      PARTHENON_REQUIRE(conduction_sat_prefac != 0.0,
                        "Saturated thermal conduction prefactor uninitialized.");
      pkg->AddParam<>("conduction_sat_prefac", conduction_sat_prefac);
    }
    pkg->AddParam<>("conduction", conduction);

    auto viscosity = Viscosity::none;
    auto viscosity_str = pin->GetOrAddString("diffusion", "viscosity", "none");
    if (viscosity_str == "isotropic") {
      viscosity = Viscosity::isotropic;
    } else if (viscosity_str != "none") {
      PARTHENON_FAIL("Unknown viscosity method. Options are: none, isotropic");
    }
    // If viscosity is enabled, process supported coefficients
    if (viscosity != Viscosity::none) {
      PARTHENON_REQUIRE_THROWS(
          typeid(parthenon::Coordinates_t) == typeid(parthenon::UniformCartesian),
          "Probably need to update derivative calc (with respect to spacing between the "
          "faces when calculating fluxes from cell centered quantities average to cell "
          "faces) in viscosity to support non-Cartesian coordiates.");
      auto viscosity_coeff_str =
          pin->GetOrAddString("diffusion", "viscosity_coeff", "none");
      auto viscosity_coeff = ViscosityCoeff::none;

      if (viscosity_coeff_str == "fixed") {
        viscosity_coeff = ViscosityCoeff::fixed;
        Real mom_diff_coeff_code = pin->GetReal("diffusion", "mom_diff_coeff_code");
        auto mom_diff = MomentumDiffusivity(viscosity, viscosity_coeff,
                                            mom_diff_coeff_code, 0.0, 0.0, 0.0);
        pkg->AddParam<>("mom_diff", mom_diff);

      } else {
        PARTHENON_FAIL("Viscosity is enabled but no coefficient is set. Please "
                       "set diffusion/viscosity_coeff to 'fixed' and "
                       "diffusion/mom_diff_coeff_code to the desired value.");
      }
    }
    pkg->AddParam<>("viscosity", viscosity);

    auto resistivity = Resistivity::none;
    auto resistivity_str = pin->GetOrAddString("diffusion", "resistivity", "none");
    if (resistivity_str == "ohmic") {
      resistivity = Resistivity::ohmic;
    } else if (resistivity_str != "none") {
      PARTHENON_FAIL("Unknown resistivity method. Options are: none, ohmic");
    }
    // If resistivity is enabled, process supported coefficients
    if (resistivity != Resistivity::none) {
      auto resistivity_coeff_str =
          pin->GetOrAddString("diffusion", "resistivity_coeff", "none");
      auto resistivity_coeff = ResistivityCoeff::none;
      // Ceiling on the ionization-model eta_O (code units; default disabled). Bounds the
      // parabolic dt in the diffusion-decoupled first-core interior; ignored for
      // resistivity_coeff=fixed. See OhmicDiffusivity::eta_cap_.
      const Real eta_ohm_cap_code = pin->GetOrAddReal(
          "diffusion", "eta_ohm_cap_code", std::numeric_limits<Real>::max());

      if (resistivity_coeff_str == "spitzer") {
        // If this is implemented, check how the Spitzer coeff for thermal conduction is
        // handled.
        PARTHENON_FAIL("needs impl");

      } else if (resistivity_coeff_str == "fixed") {
        resistivity_coeff = ResistivityCoeff::fixed;
        Real ohm_diff_coeff_code = pin->GetReal("diffusion", "ohm_diff_coeff_code");
        auto ohm_diff = OhmicDiffusivity(resistivity, resistivity_coeff,
                                         ohm_diff_coeff_code, 0.0, 0.0, 0.0);
        pkg->AddParam<>("ohm_diff", ohm_diff);

      } else if (resistivity_coeff_str == "ionization") {
        // Self-consistent Ohmic resistivity eta_O = c^2/(4 pi sigma_O) from the reduced
        // CR+thermal+grain ionization model (no ohm_diff_coeff_code needed).
        resistivity_coeff = ResistivityCoeff::ionization;
        Ionization::IonizationModel ion = BuildIonizationModel(pin);
        auto ohm_diff = OhmicDiffusivity(resistivity, resistivity_coeff, 0.0, 0.0, 0.0,
                                         0.0, ion, -1, eta_ohm_cap_code);
        pkg->AddParam<>("ohm_diff", ohm_diff);
        if (parthenon::Globals::my_rank == 0) {
          std::cout << "## Ohmic resistivity: self-consistent ionization model "
                       "(conductivity-tensor sigma_O)"
                    << std::endl;
          if (eta_ohm_cap_code < std::numeric_limits<Real>::max()) {
            std::cout << "## Ohmic eta_O cap: " << eta_ohm_cap_code << " (code units)"
                      << std::endl;
          }
        }

      } else if (resistivity_coeff_str == "ionization_chem") {
        // Ohmic eta_O = c^2/(4 pi sigma_O) from the conductivity tensor built on the
        // TIME-DEPENDENT electron abundance evolved by the chemistry package
        // (network=gow17_reduced). x_e is prim component nhydro + xe_scalar_index.
        resistivity_coeff = ResistivityCoeff::ionization_chem;
        Ionization::IonizationModel ion = BuildIonizationModel(pin);
        const int xe_scalar = pin->GetOrAddInteger("diffusion", "xe_scalar_index", 4);
        const int i_xe_prim = nhydro + xe_scalar;
        auto ohm_diff = OhmicDiffusivity(resistivity, resistivity_coeff, 0.0, 0.0, 0.0, 0.0,
                                         ion, i_xe_prim, eta_ohm_cap_code);
        pkg->AddParam<>("ohm_diff", ohm_diff);
        if (parthenon::Globals::my_rank == 0) {
          std::cout << "## Ohmic resistivity: chemistry-coupled x_e (conductivity-tensor "
                       "sigma_O, x_e from prim scalar index "
                    << i_xe_prim << ")" << std::endl;
          if (eta_ohm_cap_code < std::numeric_limits<Real>::max()) {
            std::cout << "## Ohmic eta_O cap: " << eta_ohm_cap_code << " (code units)"
                      << std::endl;
          }
        }

      } else {
        PARTHENON_FAIL("Resistivity is enabled but no valid coefficient is set. Set "
                       "diffusion/resistivity_coeff to 'fixed' (+ohm_diff_coeff_code), "
                       "'ionization', or 'ionization_chem'.");
      }
    }
    pkg->AddParam<>("resistivity", resistivity);

    auto ambipolar = Ambipolar::none;
    auto ambipolar_str = pin->GetOrAddString("diffusion", "ambipolar", "none");
    if (ambipolar_str == "ambipolar") {
      ambipolar = Ambipolar::ambipolar;
    } else if (ambipolar_str != "none") {
      PARTHENON_FAIL("Unknown ambipolar method. Options are: none, ambipolar");
    }
    // If ambipolar diffusion is enabled, process supported coefficients
    if (ambipolar != Ambipolar::none) {
      PARTHENON_REQUIRE_THROWS(fluid == Fluid::glmmhd,
                               "Ambipolar diffusion requires MHD (hydro/fluid=glmmhd).");
      auto ambipolar_coeff_str =
          pin->GetOrAddString("diffusion", "ambipolar_coeff", "none");
      auto ambipolar_coeff = AmbipolarCoeff::none;
      // Absolute ceiling on eta_A for the ionization[_chem] coefficient family (code units;
      // ignored for ambipolar_coeff=fixed). The single-fluid/tensor eta_A ~ B^2/(rho_i rho)
      // runs away (~1/rho^2) in low-density magnetized cells and sets the parabolic dt to
      // ~0 there; this bounds that runaway in the diffusion-decoupled regime, the ambipolar
      // analogue of diffusion/eta_ohm_cap_code. See AmbipolarDiffusivity::eta_cap_.
      const Real eta_ad_cap_code = pin->GetOrAddReal(
          "diffusion", "eta_ad_cap_code", std::numeric_limits<Real>::max());

      if (ambipolar_coeff_str == "fixed") {
        ambipolar_coeff = AmbipolarCoeff::fixed;
        Real ambipolar_coeff_code = pin->GetReal("diffusion", "ambipolar_coeff_code");
        auto ad_diff = AmbipolarDiffusivity(ambipolar, ambipolar_coeff,
                                            ambipolar_coeff_code, 0.0, 0.0, 0.0);
        pkg->AddParam<>("ad_diff", ad_diff);

      } else if (ambipolar_coeff_str == "ionization") {
        // Self-consistent eta_A = B^2/(4 pi gamma_AD rho_i rho) with rho_i from the
        // reduced CR+thermal+grain ionization model (replaces the constant Q_A). Defaults
        // are the shared FHC code-unit calibration; override per <diffusion> ion_* keys.
        ambipolar_coeff = AmbipolarCoeff::ionization;
        Ionization::IonizationModel ion = BuildIonizationModel(pin);
        auto ad_diff = AmbipolarDiffusivity(ambipolar, ambipolar_coeff, 0.0, 0.0, 0.0, 0.0,
                                            ion, -1, eta_ad_cap_code);
        pkg->AddParam<>("ad_diff", ad_diff);
        if (parthenon::Globals::my_rank == 0) {
          const bool tensor = ion.ad_closure == Ionization::ADClosure::tensor;
          std::cout << "## Ambipolar diffusion: self-consistent ionization model"
                    << " (CR zeta=" << ion.zeta << " s^-1, MRN dust f_dg=" << ion.f_dg
                    << ", thermal K x_K=" << ion.x_K
                    << ", AD closure=" << (tensor ? "tensor" : "single_fluid") << ")"
                    << std::endl;
          if (eta_ad_cap_code < std::numeric_limits<Real>::max()) {
            std::cout << "## Ambipolar eta_A cap: " << eta_ad_cap_code << " (code units)"
                      << std::endl;
          }
        }

      } else if (ambipolar_coeff_str == "ionization_chem") {
        // Single-fluid eta_A = B^2/(4 pi gamma_AD rho_i rho) with rho_i set by the
        // TIME-DEPENDENT electron abundance evolved by the chemistry package
        // (network=gow17_reduced). The x_e scalar is prim component nhydro + xe_scalar_index
        // (the gow17_reduced electron species sits at scalar index 4 by default).
        ambipolar_coeff = AmbipolarCoeff::ionization_chem;
        Ionization::IonizationModel ion = BuildIonizationModel(pin);
        const int xe_scalar = pin->GetOrAddInteger("diffusion", "xe_scalar_index", 4);
        const int i_xe_prim = nhydro + xe_scalar;
        auto ad_diff = AmbipolarDiffusivity(ambipolar, ambipolar_coeff, 0.0, 0.0, 0.0, 0.0,
                                            ion, i_xe_prim, eta_ad_cap_code);
        pkg->AddParam<>("ad_diff", ad_diff);
        if (parthenon::Globals::my_rank == 0) {
          std::cout << "## Ambipolar diffusion: chemistry-coupled x_e (single_fluid eta_A,"
                    << " x_e from prim scalar index " << i_xe_prim << ", gamma_AD="
                    << ion.gamma_AD << ")" << std::endl;
          if (eta_ad_cap_code < std::numeric_limits<Real>::max()) {
            std::cout << "## Ambipolar eta_A cap: " << eta_ad_cap_code << " (code units)"
                      << std::endl;
          }
        }

      } else {
        PARTHENON_FAIL("Ambipolar diffusion is enabled but no valid coefficient is set. "
                       "Set diffusion/ambipolar_coeff to 'fixed' (+ambipolar_coeff_code), "
                       "'ionization', or 'ionization_chem'.");
      }
    }
    pkg->AddParam<>("ambipolar", ambipolar);

    auto diffint_str = pin->GetOrAddString("diffusion", "integrator", "none");
    auto diffint = DiffInt::none;
    if (diffint_str == "unsplit") {
      diffint = DiffInt::unsplit;
    } else if (diffint_str == "rkl2") {
      diffint = DiffInt::rkl2;
      auto rkl2_dt_ratio = pin->GetOrAddReal("diffusion", "rkl2_max_dt_ratio", -1.0);
      pkg->AddParam<>("rkl2_max_dt_ratio", rkl2_dt_ratio);
    } else if (diffint_str != "none") {
      PARTHENON_FAIL("AthenaPK unknown integration method for diffusion processes. "
                     "Options are: none, unsplit, rkl2");
    }
    if (diffint != DiffInt::none) {
      // As in Athena++ a cfl safety factor is also applied to the theoretical limit.
      // By default it is equal to the hyperbolic cfl.
      auto cfl_diff = pin->GetOrAddReal("diffusion", "cfl", pkg->Param<Real>("cfl"));
      pkg->AddParam<>("cfl_diff", cfl_diff);
    }
    pkg->AddParam<Real>("dt_diff", std::numeric_limits<Real>::max(),
                        Params::Mutability::Mutable); // diffusive timestep constraint
    pkg->AddParam<>("diffint", diffint);

    auto hall = Hall::none;
    auto hall_str = pin->GetOrAddString("diffusion", "hall", "none");
    if (hall_str == "hall") {
      hall = Hall::hall;
    } else if (hall_str != "none") {
      PARTHENON_FAIL("Unknown hall method. Options are: none, hall");
    }
    // If the Hall effect is enabled, process supported coefficients. Hall is dispersive
    // and is never super-time-stepped: with integrator=unsplit everything is unsplit;
    // with integrator=rkl2 the parabolic terms go into RKL2 while the Hall EMF is applied
    // unsplit under its own (whistler) strict dt -- the mixed mode.
    if (hall != Hall::none) {
      PARTHENON_REQUIRE_THROWS(fluid == Fluid::glmmhd,
                               "The Hall effect requires MHD (hydro/fluid=glmmhd).");
      PARTHENON_REQUIRE_THROWS(
          diffint != DiffInt::none,
          "The Hall effect needs a diffusion integrator. Set diffusion/integrator to "
          "'unsplit' (everything unsplit) or 'rkl2' (mixed mode: parabolic terms "
          "super-time-stepped, dispersive Hall EMF unsplit at its whistler dt).");
      // Placement of the Ohmic stabilizer floor in the mixed rkl2+Hall mode:
      // 'unsplit' (default) keeps eta_floor*J with the unsplit Hall EMF, which retains
      // the strict parabolic ceiling dt <~ cfl_diff*dx^2/(6 eta_floor) on the step;
      // 'rkl2' moves the floor into the RKL2 parabolic set, lifting that ceiling so the
      // step is whistler-limited (validate against the whistler eigenmode before use).
      const auto hall_floor_int_str =
          pin->GetOrAddString("diffusion", "hall_floor_integrator", "unsplit");
      PARTHENON_REQUIRE_THROWS(hall_floor_int_str == "unsplit" ||
                                   hall_floor_int_str == "rkl2",
                               "Unknown diffusion/hall_floor_integrator. Options are: "
                               "unsplit, rkl2 (the latter only with integrator=rkl2).");
      pkg->AddParam<>("hall_floor_int_rkl2",
                      diffint == DiffInt::rkl2 && hall_floor_int_str == "rkl2");
      if (diffint == DiffInt::rkl2 && parthenon::Globals::my_rank == 0) {
        std::cout << "## Mixed diffusion integrator: parabolic terms via RKL2 STS, "
                     "dispersive Hall EMF unsplit (whistler dt strict); Ohmic floor "
                  << (hall_floor_int_str == "rkl2" ? "in the RKL2 set" : "unsplit")
                  << std::endl;
      }
      auto hall_coeff_str = pin->GetOrAddString("diffusion", "hall_coeff", "none");
      auto hall_coeff = HallCoeff::none;

      // Ceiling on |eta_H| (code units; default disabled), the Hall analogue of
      // diffusion/eta_ohm_cap_code. Bounds the strict whistler dt (~dx^2/|eta_H|) where
      // the ionization-model |eta_H| ~ B/n_e runs away in the Ohmic-dominated interior;
      // signed clamp, ignored for hall_coeff=fixed. See HallDiffusivity::eta_cap_.
      const Real eta_hall_cap_code = pin->GetOrAddReal(
          "diffusion", "eta_hall_cap_code", std::numeric_limits<Real>::max());

      if (hall_coeff_str == "fixed") {
        hall_coeff = HallCoeff::fixed;
        Real hall_coeff_code = pin->GetReal("diffusion", "hall_coeff_code");
        Real hall_ohmic_floor_code =
            pin->GetOrAddReal("diffusion", "hall_ohmic_floor_code", 0.0);
        auto hall_diff = HallDiffusivity(hall, hall_coeff, hall_coeff_code,
                                         hall_ohmic_floor_code, 0.0, 0.0, 0.0);
        pkg->AddParam<>("hall_diff", hall_diff);

        if ((resistivity == Resistivity::none) && (hall_ohmic_floor_code <= 0.0) &&
            (parthenon::Globals::my_rank == 0)) {
          std::cout << "### WARNING: The Hall effect is enabled without Ohmic resistivity "
                       "or an Ohmic floor (diffusion/hall_ohmic_floor_code). The Hall term "
                       "may be numerically unstable on the cell-centered grid."
                    << std::endl;
        }
      } else if (hall_coeff_str == "ionization") {
        // Self-consistent (signed) Hall diffusivity from the conductivity tensor of the
        // reduced CR+thermal+grain ionization model. An Ohmic floor is still recommended
        // for stability on the cell-centered grid (Hall is dispersive).
        hall_coeff = HallCoeff::ionization;
        Real hall_ohmic_floor_code =
            pin->GetOrAddReal("diffusion", "hall_ohmic_floor_code", 0.0);
        Ionization::IonizationModel ion = BuildIonizationModel(pin);
        auto hall_diff = HallDiffusivity(hall, hall_coeff, 0.0, hall_ohmic_floor_code, 0.0,
                                         0.0, 0.0, ion, -1, eta_hall_cap_code);
        pkg->AddParam<>("hall_diff", hall_diff);
        if (parthenon::Globals::my_rank == 0) {
          std::cout << "## Hall effect: self-consistent ionization model "
                       "(conductivity-tensor sigma_H, signed)"
                    << std::endl;
          if (eta_hall_cap_code < std::numeric_limits<Real>::max()) {
            std::cout << "## Hall |eta_H| cap: " << eta_hall_cap_code << " (code units)"
                      << std::endl;
          }
          if ((resistivity == Resistivity::none) && (hall_ohmic_floor_code <= 0.0)) {
            std::cout << "### WARNING: Hall enabled without Ohmic resistivity or an Ohmic "
                         "floor (diffusion/hall_ohmic_floor_code); may be unstable."
                      << std::endl;
          }
        }
      } else if (hall_coeff_str == "ionization_chem") {
        // Signed Hall diffusivity from the conductivity tensor built on the TIME-DEPENDENT
        // electron abundance evolved by the chemistry package. x_e is prim component
        // nhydro + xe_scalar_index. An Ohmic floor is still recommended for stability.
        hall_coeff = HallCoeff::ionization_chem;
        Real hall_ohmic_floor_code =
            pin->GetOrAddReal("diffusion", "hall_ohmic_floor_code", 0.0);
        Ionization::IonizationModel ion = BuildIonizationModel(pin);
        const int xe_scalar = pin->GetOrAddInteger("diffusion", "xe_scalar_index", 4);
        const int i_xe_prim = nhydro + xe_scalar;
        auto hall_diff = HallDiffusivity(hall, hall_coeff, 0.0, hall_ohmic_floor_code, 0.0,
                                         0.0, 0.0, ion, i_xe_prim, eta_hall_cap_code);
        pkg->AddParam<>("hall_diff", hall_diff);
        if (parthenon::Globals::my_rank == 0) {
          std::cout << "## Hall effect: chemistry-coupled x_e (conductivity-tensor sigma_H, "
                       "signed, x_e from prim scalar index "
                    << i_xe_prim << ")" << std::endl;
          if (eta_hall_cap_code < std::numeric_limits<Real>::max()) {
            std::cout << "## Hall |eta_H| cap: " << eta_hall_cap_code << " (code units)"
                      << std::endl;
          }
          if ((resistivity == Resistivity::none) && (hall_ohmic_floor_code <= 0.0)) {
            std::cout << "### WARNING: Hall enabled without Ohmic resistivity or an Ohmic "
                         "floor (diffusion/hall_ohmic_floor_code); may be unstable."
                      << std::endl;
          }
        }
      } else {
        PARTHENON_FAIL("The Hall effect is enabled but no valid coefficient is set. Set "
                       "diffusion/hall_coeff to 'fixed' (+hall_coeff_code), 'ionization', "
                       "or 'ionization_chem'.");
      }
    }
    pkg->AddParam<>("hall", hall);

    // Constrained Transport (divergence_control=ct) support for non-ideal MHD. Ohmic
    // resistivity (CT_AddOhmicEMF) and ambipolar diffusion (CT_AddAmbipolarEMF) run in BOTH
    // the unsplit integrator (edge EMF added to the ideal EMF each stage) and the RKL2 super-
    // time-stepping integrator (the divergence-free induction is super-time-stepped on the
    // face field Bf with M(Bf)=-curl(E_diff); see AddSTSTasks) -- both validated to analytic.
    // The Hall effect (CT_AddHallEMF) uses the COMPACT edge-current formulation (tight edge
    // curls; transverse J interpolated only in transverse directions), applied unsplit in the
    // main step. Div-free, stable (with the usual Ohmic floor), and now validated: the CT
    // whistler dispersion matches the GLM reference to 3 sig figs across Q_H (fast branch,
    // 4.7% at Q_H=0.05, 1.9% at Q_H=0.2 -- identical to GLM). The earlier "dispersion deficit"
    // was a task-graph race (CT_UpdateBf did not depend on the Hall EMF task and curled a pre-
    // Hall edge EMF); fixed by gating CT_UpdateBf on `emf` in hydro_driver.cpp. No warning.

    // Fused non-ideal diffusive-timestep path. The per-cell Wardle conductivity tensor
    // (Ionization::Diffusivities) is the dominant GPU cost, and the three per-term dt
    // reductions each recompute the full tensor but keep only one eta. When EVERY active
    // non-ideal term uses an ionization-family coefficient ("ionization" or the
    // chem-coupled "ionization_chem"), the reductions can share the minimal set of
    // tensor evaluations per cell (FusedNonidealEval: one equilibrium solve, plus a
    // chem-x_e solve only where Ohm/Hall are chem-coupled; the tier-4 production mix of
    // equilibrium Ohm/Hall + chem-capped AD needs just ONE solve where the per-term path
    // ran three). Fixed-coeff (cheap) terms keep their exact per-term estimators.
    {
      bool any_nonideal = false;
      bool all_ionization = true;
      if (resistivity != Resistivity::none) {
        any_nonideal = true;
        const auto ct = pkg->Param<OhmicDiffusivity>("ohm_diff").GetCoeffType();
        if (ct != ResistivityCoeff::ionization && ct != ResistivityCoeff::ionization_chem)
          all_ionization = false;
      }
      if (ambipolar != Ambipolar::none) {
        any_nonideal = true;
        const auto ct = pkg->Param<AmbipolarDiffusivity>("ad_diff").GetCoeffType();
        if (ct != AmbipolarCoeff::ionization && ct != AmbipolarCoeff::ionization_chem)
          all_ionization = false;
      }
      if (hall != Hall::none) {
        any_nonideal = true;
        const auto ct = pkg->Param<HallDiffusivity>("hall_diff").GetCoeffType();
        if (ct != HallCoeff::ionization && ct != HallCoeff::ionization_chem)
          all_ionization = false;
      }
      // In the mixed rkl2+Hall mode the whistler limit must stay separate from the
      // parabolic aggregate; the two-reducer FusedMixed estimator handles that (one
      // tensor evaluation per cell, two minima), so fusion is available there too.
      // NOTE (profiled 2026-07-11, job 2321477): the per-term estimators were ~15% of
      // all GPU kernel time in the mixed production run (2x843ms + 2x421ms per cycle).
      const bool fusable = any_nonideal && all_ionization;
      // Optional override (default = auto-on when fusable). Setting
      // diffusion/fused_nonideal_dt=false forces the exact per-term estimators (used for
      // A/B validation and as an escape hatch); it can never enable fusion where the
      // coefficients are not all ionization-family.
      const bool fused = pin->GetOrAddBoolean("diffusion", "fused_nonideal_dt", fusable);
      pkg->AddParam<>("nonideal_dt_fused", fused && fusable);
    }

    // Cell-centered non-ideal diffusivity cache. When enabled, CalcDiffFluxes first
    // fills a 3-component derived field (eta_O, eta_H, eta_A) once per stage
    // (PrecomputeNonidealEta) and the Ohmic/ambipolar/Hall flux kernels face-average
    // the cached per-cell values instead of re-running the (expensive) ionization
    // solve per face -- cutting up to 9 tensor evaluations per cell per stage to 1.
    // This is the Athena++ pattern (cell-centered CalcDiffusivity, then face/edge
    // averaging). Face eta becomes the arithmetic mean of the two adjacent cell
    // values rather than eta(face-averaged state): same order of accuracy, but not
    // bit-identical -- diffusion/eta_cache=false restores the per-face evaluation
    // for A/B validation. Default: on when any active term uses an ionization-type
    // coefficient (where it pays), off for purely fixed coefficients (where the
    // per-face evaluation is trivially cheap and bit-identical history matters).
    {
      const bool any_nonideal = (resistivity != Resistivity::none) ||
                                (ambipolar != Ambipolar::none) || (hall != Hall::none);
      bool any_ionization = false;
      if (resistivity != Resistivity::none) {
        const auto ct = pkg->Param<OhmicDiffusivity>("ohm_diff").GetCoeffType();
        any_ionization |= (ct == ResistivityCoeff::ionization ||
                           ct == ResistivityCoeff::ionization_chem);
      }
      if (ambipolar != Ambipolar::none) {
        const auto ct = pkg->Param<AmbipolarDiffusivity>("ad_diff").GetCoeffType();
        any_ionization |= (ct == AmbipolarCoeff::ionization ||
                           ct == AmbipolarCoeff::ionization_chem);
      }
      if (hall != Hall::none) {
        const auto ct = pkg->Param<HallDiffusivity>("hall_diff").GetCoeffType();
        any_ionization |=
            (ct == HallCoeff::ionization || ct == HallCoeff::ionization_chem);
      }
      const bool eta_cache =
          pin->GetOrAddBoolean("diffusion", "eta_cache", any_ionization) && any_nonideal;
      // Under eos=hydrogen the correct code temperature for the ionization model comes from
      // the EOS table; that plumbing exists on the cached path (PrecomputeNonidealEta + the
      // fused dt estimator) but not the per-face cache-off flux kernels, so require the
      // cache there (it is the default when any term uses an ionization coefficient).
      PARTHENON_REQUIRE(!(use_h2diss && any_ionization) || eta_cache,
                        "eos=hydrogen with ionization-coefficient non-ideal MHD requires "
                        "diffusion/eta_cache=true (the default).");
      pkg->AddParam<>("nonideal_eta_cache", eta_cache);
      if (eta_cache) {
        Metadata m_eta({Metadata::Cell, Metadata::Derived, Metadata::OneCopy},
                       std::vector<int>({3}),
                       std::vector<std::string>{"eta_ohmic", "eta_hall", "eta_ambipolar"});
        pkg->AddField("nonideal_eta", m_eta);
        if (parthenon::Globals::my_rank == 0) {
          std::cout << "## Non-ideal eta cache: cell-centered (eta_O, eta_H, eta_A) "
                       "precomputed once per stage; flux kernels use face-averaged "
                       "cached values (diffusion/eta_cache=false to disable)."
                    << std::endl;
        }
      }
      // FLAGSHIP AUDIT ITEM 5 -- WS-4 DUST -> CONDUCTIVITY COUPLING.
      // Grains dominate recombination over most of the collapse, so they set x_e and hence
      // all three non-ideal diffusivities. Until now the ionization model used a FROZEN ISM
      // MRN population while the dust package evolved (f_dg, a_c) for the OPACITY only --
      // i.e. the run's own grain evolution never reached its own magnetic microphysics.
      // With this on, the conductivity grain population is rescaled per cell to the evolved
      // scalars (see IonizationModel::FDG/Ak and FusedNonidealEval::CellIon).
      // OFF by default; when off nothing is copied and results are bit-identical.
      const bool dust_couple =
          pin->GetOrAddBoolean("diffusion", "dust_coupling", false);
      pkg->AddParam<>("nonideal_dust_on", dust_couple);
      if (dust_couple) {
        PARTHENON_REQUIRE(pin->GetOrAddBoolean("dust", "evolve", false),
                          "diffusion/dust_coupling=true requires the dust package to be "
                          "evolving (<dust> evolve=true) -- otherwise the grain population "
                          "is static and the coupling is a no-op with extra cost.");
        PARTHENON_REQUIRE(any_ionization,
                          "diffusion/dust_coupling=true only affects ionization-family "
                          "coefficients (resistivity/ambipolar/hall _coeff = ionization*).");
        // Dust stores f_dg at scalar_index and a_c at scalar_index+1 (see dust.cpp);
        // prim component = nhydro + scalar index, the same convention as x_e.
        const int dsi = pin->GetOrAddInteger("dust", "scalar_index", 0);
        pkg->AddParam<>("nonideal_dust_i_fdg", nhydro + dsi);
        pkg->AddParam<>("nonideal_dust_i_ac", nhydro + dsi + 1);
        // Reference population the STATIC MRN bins represent. a_ref/f_dg_ref must be the
        // values the bins were built from, or the scaling is biased at t=0.
        // NOTE the key is "a_ref_cm", not "a_ref" (see dust.cpp:38). Reading the wrong key
        // would silently fall back to 1e-5 and bias a_scale for any deck that sets it.
        pkg->AddParam<>("nonideal_dust_a_ref",
                        pin->GetOrAddReal("dust", "a_ref_cm", 1.0e-5));
        pkg->AddParam<>("nonideal_dust_fdg_ref",
                        pin->GetOrAddReal("dust", "f_dg_ref", 0.01));
        if (parthenon::Globals::my_rank == 0) {
          std::cout << "## Dust->conductivity coupling ON (diffusion/dust_coupling): the "
                       "ionization grain population is rescaled per cell from the evolved "
                       "(f_dg, a_c) at prim indices " << nhydro + dsi << ", "
                    << nhydro + dsi + 1 << "." << std::endl;
        }
      }
      // FLAGSHIP AUDIT ITEM 1 -- non-ideal cap-activation diagnostic.
      // eta_ohm_cap_code / eta_ad_cap_code / eta_hall_cap_code are numerical stabilizers
      // that MODIFY the induction equation, so a flux-retention number is uninterpretable
      // without knowing where and how hard they bit. Filled inside PrecomputeNonidealEta
      // from the same tensor solve (no extra solve, no extra task); OFF by default, and
      // bit-identical when off (it only writes a Derived field, never feeds back).
      const bool cap_diag = pin->GetOrAddBoolean("diffusion", "cap_diag", false);
      PARTHENON_REQUIRE(!cap_diag || eta_cache,
                        "diffusion/cap_diag=true requires diffusion/eta_cache=true.");
      pkg->AddParam<>("nonideal_cap_diag", cap_diag);
      if (cap_diag) {
        Metadata m_cap({Metadata::Cell, Metadata::Derived, Metadata::OneCopy},
                       std::vector<int>({CapDiagIdx::NCOMP}),
                       std::vector<std::string>{"cap_flag_O", "cap_flag_H", "cap_flag_A",
                                                "cap_dec_O", "cap_dec_H", "cap_dec_A"});
        pkg->AddField("diff.capdiag", m_cap);
        auto hst_cap = pkg->Param<parthenon::HstVar_list>(parthenon::hist_param_key);
        using parthenon::HistoryOutputVar;
        using parthenon::UserHistoryOperation;
        // Normalizer: the summed volume of exactly the cells the flags are summed over.
        // Volume fraction capped = cap-V<term> / cap-Vtot; mass fraction = cap-M<term>/mass;
        // volume-weighted mean clipped decades over the capped set = cap-D<term>/cap-V<term>.
        hst_cap.emplace_back(HistoryOutputVar(UserHistoryOperation::sum,
                                              NonidealCapHstVol, "cap-Vtot"));
        const std::pair<int, const char *> terms[3] = {{CapDiagIdx::flagO, "O"},
                                                       {CapDiagIdx::flagH, "H"},
                                                       {CapDiagIdx::flagA, "A"}};
        for (const auto &tm : terms) {
          const int c = tm.first;
          hst_cap.emplace_back(HistoryOutputVar(
              UserHistoryOperation::sum,
              [c](MeshData<Real> *md) { return NonidealCapHstFlag(md, c); },
              std::string("cap-V") + tm.second));
          hst_cap.emplace_back(HistoryOutputVar(
              UserHistoryOperation::sum,
              [c](MeshData<Real> *md) { return NonidealCapHstMassFlag(md, c); },
              std::string("cap-M") + tm.second));
          hst_cap.emplace_back(HistoryOutputVar(
              UserHistoryOperation::sum,
              [c](MeshData<Real> *md) { return NonidealCapHstDecade(md, c); },
              std::string("cap-D") + tm.second));
        }
        pkg->UpdateParam(parthenon::hist_param_key, hst_cap);
        if (parthenon::Globals::my_rank == 0) {
          std::cout << "## Non-ideal cap diagnostic ON (diffusion/cap_diag): hst cap-V*/"
                       "cap-M*/cap-D* + field diff.capdiag (add to an output block for "
                       "the spatial picture)."
                    << std::endl;
        }
      }
      // Freeze the cached eta across the RKL2 stages of each Strang half: the stage
      // tasks skip the PrecomputeNonidealEta refresh, so eta stays at its value from
      // the half's first stage (refreshed 2x per cycle instead of once per stage).
      // RKL2's stability polynomial (Meyer+2012) is derived for an operator held fixed
      // across the super-step, and the coefficient lag is O(dt/2) on eta(rho,T,B),
      // which evolves on the collapse timescale >> dt. Requires the eta cache; only
      // meaningful with integrator=rkl2. Default off (bit-identical history matters
      // for A/B validation).
      const bool freeze_eta =
          pin->GetOrAddBoolean("diffusion", "rkl2_freeze_eta", false);
      PARTHENON_REQUIRE(!freeze_eta || (eta_cache && diffint == DiffInt::rkl2),
                        "diffusion/rkl2_freeze_eta=true requires diffusion/eta_cache=true "
                        "and diffusion/integrator=rkl2.");
      pkg->AddParam<>("rkl2_freeze_eta", freeze_eta && eta_cache);
      if (freeze_eta && parthenon::Globals::my_rank == 0) {
        std::cout << "## RKL2 eta freeze: non-ideal eta cache refreshed at the first "
                     "stage of each Strang half only (diffusion/rkl2_freeze_eta)."
                  << std::endl;
      }
    }

    if (fluid == Fluid::euler) {
      PARTHENON_REQUIRE(!use_h2diss,
                        "H2-dissociation EOS (eos=hydrogen) is wired for fluid=glmmhd only.");
      AdiabaticHydroEOS eos(pfloor, dfloor, efloor, vceil, eceil, gamma);
      pkg->AddParam<>("eos", eos);
      pkg->FillDerivedMesh = ConsToPrim<AdiabaticHydroEOS>;
      pkg->EstimateTimestepMesh = EstimateTimestep<Fluid::euler>;
    } else if (fluid == Fluid::glmmhd) {
      // Boris / semi-relativistic Alfven-speed limiter: cap c_b on the Alfven speed used in
      // the Maxwell-stress FORCE + signal speeds (via the EOS + a per-interface LLF-Boris
      // branch in the HLLD solver). Relaxes the timestep / tames the vA runaway in evacuated
      // strongly-magnetized cells (the first-core magnetic wall). 0 (default) = disabled,
      // bit-identical.
      const Real boris_ca_max =
          pin->GetOrAddReal("hydro", "boris_ca_max_code", 0.0);
      AdiabaticGLMMHDEOS eos(pfloor, dfloor, efloor, vceil, eceil, gamma, use_h2diss,
                             eos_tab, boris_ca_max, ct_glm_inert);
      pkg->AddParam<>("eos", eos);
      pkg->FillDerivedMesh = ConsToPrim<AdiabaticGLMMHDEOS>;
      pkg->EstimateTimestepMesh = EstimateTimestep<Fluid::glmmhd>;
      if (boris_ca_max > 0.0 && parthenon::Globals::my_rank == 0) {
        std::cout << "## Boris Alfven-speed limiter: c_b = " << boris_ca_max
                  << " (code units); HLLD->LLF-Boris where vA>0.33 c_b." << std::endl;
      }
    }
  } else {
    PARTHENON_FAIL("AthenaPK hydro: Unknown EOS");
  }

  /************************************************************
   * Read Tabular Cooling
   ************************************************************/

  const auto enable_cooling_str =
      pin->GetOrAddString("cooling", "enable_cooling", "none");

  auto cooling = Cooling::none;
  if (enable_cooling_str == "tabular") {
    cooling = Cooling::tabular;
  } else if (enable_cooling_str != "none") {
    PARTHENON_FAIL("AthenaPK hydro: Unknown cooling string. Supported options are "
                   "'none' and 'tabular'");
  }
  pkg->AddParam<>("enable_cooling", cooling);

  if (cooling == Cooling::tabular) {
    TabularCooling tabular_cooling(pin, pkg);
    pkg->AddParam<>("tabular_cooling", tabular_cooling);
  }

  auto scratch_level = pin->GetOrAddInteger("hydro", "scratch_level", 0);
  pkg->AddParam("scratch_level", scratch_level);

  auto nscalars = pin->GetOrAddInteger("hydro", "nscalars", 0);
  pkg->AddParam("nscalars", nscalars);

  std::vector<std::string> cons_labels(nhydro);
  cons_labels[IDN] = "density";
  cons_labels[IM1] = "momentum_density_1";
  cons_labels[IM2] = "momentum_density_2";
  cons_labels[IM3] = "momentum_density_3";
  cons_labels[IEN] = "total_energy_density";
  if (fluid == Fluid::glmmhd) {
    cons_labels[IB1] = "magnetic_field_1";
    cons_labels[IB2] = "magnetic_field_2";
    cons_labels[IB3] = "magnetic_field_3";
    cons_labels[IPS] = "magnetic_psi";
  }

  // TODO(pgrete) check if this could be "one-copy" for two stage SSP integrators
  std::vector<std::string> prim_labels(nhydro);
  prim_labels[IDN] = "density";
  prim_labels[IV1] = "velocity_1";
  prim_labels[IV2] = "velocity_2";
  prim_labels[IV3] = "velocity_3";
  prim_labels[IPR] = "pressure";
  if (fluid == Fluid::glmmhd) {
    prim_labels[IB1] = "magnetic_field_1";
    prim_labels[IB2] = "magnetic_field_2";
    prim_labels[IB3] = "magnetic_field_3";
    prim_labels[IPS] = "magnetic_psi";
  }
  for (auto i = 0; i < nscalars; i++) {
    cons_labels.emplace_back("scalar_density_" + std::to_string(i));
    prim_labels.emplace_back("scalar_" + std::to_string(i));
  }

  Metadata m(
      {Metadata::Cell, Metadata::Independent, Metadata::FillGhost, Metadata::WithFluxes},
      std::vector<int>({nhydro + nscalars}), cons_labels);
  m.RegisterRefinementOps<refinement_ops::ProlongateCellMinModMultiD,
                          parthenon::refinement_ops::RestrictAverage>();
  pkg->AddField("cons", m);

  m = Metadata({Metadata::Cell, Metadata::Derived}, std::vector<int>({nhydro + nscalars}),
               prim_labels);
  pkg->AddField("prim", m);

  // CT (Phase 2): face-centered primary magnetic field. Registered with the div-free
  // AMR refinement ops (Toth & Roe internal prolongation) even though increment 1 is
  // single level, per docs/CT_DESIGN.md sec 3.1. The edge-flux slots (WithFluxes) hold
  // the CT EMF; FillGhost exchanges the staggered ghosts.
  if (use_ct) {
    auto m_bf = Metadata({Metadata::Face, Metadata::Independent, Metadata::WithFluxes,
                          Metadata::FillGhost});
    m_bf.RegisterRefinementOps<parthenon::refinement_ops::ProlongateSharedMinMod,
                               parthenon::refinement_ops::RestrictAverage,
                               parthenon::refinement_ops::ProlongateInternalTothAndRoe>();
    pkg->AddField<Hydro::CT::Bf>(m_bf);

    // Face-divergence history diagnostic (max |div B|*dx/|B|); ~round-off on the CT path.
    auto hst_vars = pkg->Param<parthenon::HstVar_list>(parthenon::hist_param_key);
    // Projection diagnostic (hydro/ct_proj_diag): per-cell relative internal-energy change
    // imposed by the face->cell projection, i.e. the direct measurement of the E-vs-ME
    // bookkeeping mismatch. Derived => not carried in restarts; opt-in so production runs
    // pay no memory. Dump "ct.dEint" in an output block to localize it spatially.
    if (pkg->Param<bool>("ct_proj_diag")) {
      auto m_diag = Metadata({Metadata::Cell, Metadata::Derived, Metadata::OneCopy});
      pkg->AddField("ct.dEint", m_diag);
      hst_vars.emplace_back(parthenon::HistoryOutputVar(
          parthenon::UserHistoryOperation::min, Hydro::CT::CT_ProjEintMin,
          "ct_projEintMin"));
      hst_vars.emplace_back(parthenon::HistoryOutputVar(
          parthenon::UserHistoryOperation::max, Hydro::CT::CT_ProjEintMaxAbs,
          "ct_projEintMaxAbs"));
    }
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::max, Hydro::CT::CT_MaxRelFaceDivB, "ct_maxRelDivB"));
    hst_vars.emplace_back(parthenon::HistoryOutputVar(
        parthenon::UserHistoryOperation::max, Hydro::CT::CT_MaxAbsFaceDivB, "ct_maxAbsDivB"));
    pkg->UpdateParam(parthenon::hist_param_key, hst_vars);
  }

  // FLAGSHIP AUDIT ITEM 2 -- magnetic-transport (fossil-field) history diagnostics.
  // See src/diagnostics/mag_diag.hpp. Pure read-only volume reductions => bit-identical
  // when off, and off by default. MHD only (they all need B).
  const bool mag_diag =
      pin->GetOrAddBoolean("hydro", "mag_diag", false) && fluid == Fluid::glmmhd;
  pkg->AddParam<>("mag_diag", mag_diag);
  // WP-8: density threshold (CODE density) splitting the dissipation integrals into a core
  // and an envelope budget. 0 = OFF, which registers no extra columns, so the OFF state and
  // the original column set stay bit-identical. Always added because MagDiagReduce reads it
  // unconditionally. Sensible value: ~1e2-1e4 code, i.e. well above the t=0 peak (~5 code)
  // and below rhocrit (1.83e5 code) -- see mag_diag.hpp for why a single global integral of
  // eta|J|^2 is not a convergent diagnostic.
  pkg->AddParam<Real>("mag_diag_rho_split",
                      pin->GetOrAddReal("hydro", "mag_diag_rho_split", 0.0));
  if (mag_diag) {
    using Diagnostics::MagDiag;
    using Diagnostics::MagDiagReduce;
    auto hst_mag = pkg->Param<parthenon::HstVar_list>(parthenon::hist_param_key);
    auto add_mag = [&hst_mag](MagDiag w, const std::string &label) {
      hst_mag.emplace_back(parthenon::HistoryOutputVar(
          parthenon::UserHistoryOperation::sum,
          [w](MeshData<Real> *md) { return MagDiagReduce(md, w); }, label));
    };
    add_mag(MagDiag::Jsq, "mag-Jsq");
    add_mag(MagDiag::Hc, "mag-Hc");
    add_mag(MagDiag::MEtor, "mag-MEtor");
    add_mag(MagDiag::MEpol, "mag-MEpol");
    // eta-weighted dissipation needs the cell-centered eta cache; the cap-restricted
    // variants additionally need the cap flags. Both params only exist when a non-ideal
    // term is active, so probe rather than assume (mag_diag is legal on an ideal run).
    const bool have_eta = pkg->AllParams().hasKey("nonideal_eta_cache") &&
                          pkg->Param<bool>("nonideal_eta_cache");
    const bool have_cap = pkg->AllParams().hasKey("nonideal_cap_diag") &&
                          pkg->Param<bool>("nonideal_cap_diag");
    const Real rho_split = pkg->Param<Real>("mag_diag_rho_split");
    if (have_eta) {
      add_mag(MagDiag::dissO, "mag-dissO");
      add_mag(MagDiag::dissA, "mag-dissA");
      // "how much of the dissipation is a numerical stabilizer" needs the cap flags.
      if (have_cap) {
        add_mag(MagDiag::dissOcap, "mag-dissOcap");
        add_mag(MagDiag::dissAcap, "mag-dissAcap");
      }
      // WP-8 columns are ALL behind rho_split>0, deliberately. Appending columns shifts
      // every later column index (maxRelDivB would move from 31 to 33), which would break
      // both the existing analysis scripts and any comparison against the history files
      // already on disk from the njeans ladder. With rho_split=0 (the default) the column
      // set is byte-for-byte what it was.
      if (rho_split > 0.0) {
        // The carrying volume -- what makes the ill-conditioning visible: 90% of the global
        // integral was measured to come from a volume fraction of ~1e-7.
        add_mag(MagDiag::dissOsq, "mag-dissOsq");
        add_mag(MagDiag::dissAsq, "mag-dissAsq");
        add_mag(MagDiag::dissOhi, "mag-dissO-hi");
        add_mag(MagDiag::dissOlo, "mag-dissO-lo");
        add_mag(MagDiag::dissAhi, "mag-dissA-hi");
        add_mag(MagDiag::dissAlo, "mag-dissA-lo");
        add_mag(MagDiag::Vhi, "mag-Vhi");
      }
    }
    pkg->UpdateParam(parthenon::hist_param_key, hst_mag);
    if (parthenon::Globals::my_rank == 0) {
      std::cout << "## Magnetic-transport diagnostics ON (hydro/mag_diag): hst mag-Jsq, "
                   "mag-Hc (current helicity), mag-MEtor/MEpol"
                << (have_eta ? ", mag-dissO/dissA" : "")
                << (have_eta && have_cap ? ", mag-dissOcap/dissAcap" : "")
                << "." << std::endl;
      if (have_eta && rho_split > 0.0) {
        std::cout << "## mag_diag density split ON (hydro/mag_diag_rho_split=" << rho_split
                  << " code): hst mag-dissO-hi/lo, mag-dissA-hi/lo, mag-Vhi. Use these,"
                     " NOT the global mag-dissO/dissA, for convergence work -- WP-8 measured"
                     " 90% of the global integral coming from ~1e-7 of the volume."
                  << std::endl;
      }
    }
  }

  // WP-4 -- angular-momentum history diagnostics. See src/diagnostics/angmom_diag.hpp.
  // The production history had NO angular-momentum column at all, which for a
  // magnetic-braking result is the central quantity. Pure read-only reductions => the OFF
  // state is bit-identical, and OFF is the default.
  const bool angmom_diag = pin->GetOrAddBoolean("hydro", "angmom_diag", false);
  pkg->AddParam<>("angmom_diag", angmom_diag);
  const Real angmom_rho_split =
      pin->GetOrAddReal("hydro", "angmom_diag_rho_split", 0.0);
  if (angmom_diag) {
    using Diagnostics::AngMom;
    using Diagnostics::AngMomReduce;
    auto hst_am = pkg->Param<parthenon::HstVar_list>(parthenon::hist_param_key);
    auto add_am = [&hst_am](AngMom w, const std::string &label) {
      hst_am.emplace_back(parthenon::HistoryOutputVar(
          parthenon::UserHistoryOperation::sum,
          [w](MeshData<Real> *md) { return AngMomReduce(md, w); }, label));
    };
    pkg->AddParam<bool>("angmom_diag_mhd", fluid == Fluid::glmmhd);
    pkg->AddParam<Real>("angmom_diag_rho_split", angmom_rho_split);
    add_am(AngMom::Lx, "am-Lx");
    add_am(AngMom::Ly, "am-Ly");
    add_am(AngMom::Lz, "am-Lz");
    // THE budget term: the TOTAL stress flux, including the pressure and Maxwell parts.
    // Measured 2026-07-31: with the advective part alone the budget misses by >2 orders of
    // magnitude, because the surface pressure torque dominates under outflow BCs.
    add_am(AngMom::FTx, "am-FTx");
    add_am(AngMom::FTy, "am-FTy");
    add_am(AngMom::FTz, "am-FTz");
    // Advective part alone -- kept because it separates "L was carried off by outflowing
    // gas" from "L was removed by a surface stress", which are different physical claims.
    add_am(AngMom::FLx, "am-FLx");
    add_am(AngMom::FLy, "am-FLy");
    add_am(AngMom::FLz, "am-FLz");
    // Gravitational torque -- zero in the continuum for an isolated system, so its measured
    // size probes the Poisson solve rather than the physics. Needs a potential to exist.
    if (pin->GetOrAddBoolean("physics", "self_gravity", false)) {
      add_am(AngMom::TgravX, "am-Tgravx");
      add_am(AngMom::TgravY, "am-Tgravy");
      add_am(AngMom::TgravZ, "am-Tgravz");
    }
    // The Lorentz torque needs a field. INTERPRETIVE ONLY -- do not add to am-FT*, which
    // already contains the Maxwell stress; see angmom_diag.hpp on the double-count.
    if (fluid == Fluid::glmmhd) {
      add_am(AngMom::TmagX, "am-Tmagx");
      add_am(AngMom::TmagY, "am-Tmagy");
      add_am(AngMom::TmagZ, "am-Tmagz");
    }
    // Density-split L (WP-4 follow-up). ALL of these sit behind rho_split > 0, deliberately
    // and for the same reason as mag_diag's WP-8 columns: appending history columns shifts
    // every downstream index and would silently break parsers reading .hst files already on
    // disk. With the default rho_split = 0 the column set is bit-identical.
    if (angmom_rho_split > 0.0) {
      add_am(AngMom::Lxhi, "am-Lx-hi");
      add_am(AngMom::Lxlo, "am-Lx-lo");
      add_am(AngMom::Lyhi, "am-Ly-hi");
      add_am(AngMom::Lylo, "am-Ly-lo");
      add_am(AngMom::Lzhi, "am-Lz-hi");
      add_am(AngMom::Lzlo, "am-Lz-lo");
      add_am(AngMom::Mhi, "am-Mhi");
    }
    pkg->UpdateParam(parthenon::hist_param_key, hst_am);
    if (parthenon::Globals::my_rank == 0) {
      const auto &ms = pin->GetOrAddReal("parthenon/mesh", "x1min", 0.0);
      const auto &me = pin->GetOrAddReal("parthenon/mesh", "x1max", 0.0);
      std::cout << "## Angular-momentum diagnostics ON (hydro/angmom_diag): hst am-Lx/Ly/Lz,"
                   " am-FLx/FLy/FLz (advective flux OUT through the domain faces,"
                   " outward-positive)"
                << (fluid == Fluid::glmmhd ? ", am-Tmagx/y/z (Lorentz torque)" : "")
                << ". Moments about the BOX CENTRE, x1 centre = " << 0.5 * (ms + me)
                << " code. Budget: d(am-L)/dt = am-Tmag - am-FL + (pressure and gravity"
                   " torques, NOT instrumented -- see angmom_diag.hpp)."
                << std::endl;
    }
  }

  // WP-6 -- conservation-budget history diagnostics. See src/diagnostics/cons_diag.hpp.
  // These exist because tot-E excludes the gravitational potential energy (so energy
  // conservation is currently UNVERIFIABLE, not merely unverified), mass RISES under
  // outflow BCs, and 1-mom drifts ~10%. All read-only => bit-identical when off.
  const bool cons_diag = pin->GetOrAddBoolean("hydro", "cons_diag", false);
  pkg->AddParam<>("cons_diag", cons_diag);
  if (cons_diag) {
    using Diagnostics::ConsDiag;
    using Diagnostics::ConsDiagReduce;
    // Cached for the reduction kernel: the floor values and whether B exists. Read from
    // pin with the SAME defaults hydro.cpp uses above, so a deck that omits them gets the
    // same -1.0 "disabled" sentinel here and the floor columns stay identically zero.
    pkg->AddParam<Real>("cons_diag_dfloor", pin->GetOrAddReal("hydro", "dfloor", -1.0));
    pkg->AddParam<Real>("cons_diag_pfloor", pin->GetOrAddReal("hydro", "pfloor", -1.0));
    pkg->AddParam<bool>("cons_diag_mhd", fluid == Fluid::glmmhd);

    auto hst_c = pkg->Param<parthenon::HstVar_list>(parthenon::hist_param_key);
    auto add_c = [&hst_c](ConsDiag w, const std::string &label) {
      hst_c.emplace_back(parthenon::HistoryOutputVar(
          parthenon::UserHistoryOperation::sum,
          [w](MeshData<Real> *md) { return ConsDiagReduce(md, w); }, label));
    };
    // W = 1/2 int rho Phi dV needs the solver's potential field to exist.
    const bool have_phi = pin->GetOrAddBoolean("physics", "self_gravity", false);
    if (have_phi) add_c(ConsDiag::Wgrav, "cons-W");
    add_c(ConsDiag::Mout, "cons-Mout");
    // The factor-2 experiment: the same faces read from the SOLVER's own flux array.
    // Agreement with cons-Mout localises the discrepancy to the time integration;
    // a factor of 2 between them localises it to the surface integral. See cons_diag.hpp.
    add_c(ConsDiag::MoutSolver, "cons-Mout-solver");
    add_c(ConsDiag::MoutSolverInner, "cons-Mout-solver-in");
    add_c(ConsDiag::MoutSolverOuter, "cons-Mout-solver-out");
    add_c(ConsDiag::MoutInner, "cons-Mout-in");
    add_c(ConsDiag::MoutOuter, "cons-Mout-out");
    add_c(ConsDiag::PoutX, "cons-Poutx");
    add_c(ConsDiag::PoutY, "cons-Pouty");
    add_c(ConsDiag::PoutZ, "cons-Poutz");
    add_c(ConsDiag::nfloor, "cons-nfloor");
    add_c(ConsDiag::Mfloor, "cons-Mfloor");
    add_c(ConsDiag::npfloor, "cons-npfloor");
    pkg->UpdateParam(parthenon::hist_param_key, hst_c);
    if (parthenon::Globals::my_rank == 0) {
      std::cout << "## Conservation-budget diagnostics ON (hydro/cons_diag): hst "
                << (have_phi ? "cons-W (grav. PE), " : "")
                << "cons-Mout, cons-Poutx/y/z (TOTAL momentum flux incl. pressure and"
                   " Maxwell stress), cons-nfloor/Mfloor/npfloor. All *out columns are"
                   " OUTWARD-POSITIVE, so a NEGATIVE cons-Mout is inflow through an"
                   " outflow face. Budgets: d(mass)/dt = -cons-Mout, d(i-mom)/dt ="
                   " -cons-Pouti. Floor columns are a PROXY for floor activity, NOT the"
                   " mass injected -- see cons_diag.hpp."
                << std::endl;
      if (!have_phi) {
        std::cout << "## cons_diag: self_gravity is off, so cons-W is not registered."
                  << std::endl;
      }
    }
  }

  const auto refine_str = pin->GetOrAddString("refinement", "type", "unset");
  if (refine_str == "pressure_gradient") {
    pkg->CheckRefinementBlock = refinement::gradient::PressureGradient;
    const auto thr = pin->GetOrAddReal("refinement", "threshold_pressure_gradient", 0.0);
    PARTHENON_REQUIRE(thr > 0.,
                      "Make sure to set refinement/threshold_pressure_gradient >0.");
    pkg->AddParam<Real>("refinement/threshold_pressure_gradient", thr);
  } else if (refine_str == "xyvelocity_gradient") {
    pkg->CheckRefinementBlock = refinement::gradient::VelocityGradient;
    const auto thr =
        pin->GetOrAddReal("refinement", "threshold_xyvelocity_gradient", 0.0);
    PARTHENON_REQUIRE(thr > 0.,
                      "Make sure to set refinement/threshold_xyvelocity_gradient >0.");
    pkg->AddParam<Real>("refinement/threshold_xyvelocity_gradient", thr);
  } else if (refine_str == "maxdensity") {
    pkg->CheckRefinementBlock = refinement::other::MaxDensity;
    const auto deref_below =
        pin->GetOrAddReal("refinement", "maxdensity_deref_below", 0.0);
    const auto refine_above =
        pin->GetOrAddReal("refinement", "maxdensity_refine_above", 0.0);
    PARTHENON_REQUIRE(deref_below > 0.,
                      "Make sure to set refinement/maxdensity_deref_below > 0.");
    PARTHENON_REQUIRE(refine_above > 0.,
                      "Make sure to set refinement/maxdensity_refine_above > 0.");
    PARTHENON_REQUIRE(deref_below < refine_above,
                      "Make sure to set refinement/maxdensity_deref_below < "
                      "refinement/maxdensity_refine_above");
    pkg->AddParam<Real>("refinement/maxdensity_deref_below", deref_below);
    pkg->AddParam<Real>("refinement/maxdensity_refine_above", refine_above);
  } else if (refine_str == "jeans") {
    pkg->CheckRefinementBlock = refinement::jeans::Jeans;
    const auto njeans = pin->GetOrAddReal("refinement", "njeans", 0.0);
    PARTHENON_REQUIRE(njeans > 0.,
                      "Make sure to set refinement/njeans > 0 (typically 8-16).");
    pkg->AddParam<Real>("refinement/njeans", njeans);
  } else if (refine_str == "jeans_nonideal") {
    // Flagship Phase 7 physics-based AMR: Jeans + current-sheet resolution (opt-in).
    pkg->CheckRefinementBlock = refinement::nonideal::JeansNonideal;
    const auto njeans = pin->GetOrAddReal("refinement", "njeans", 0.0);
    PARTHENON_REQUIRE(njeans > 0.,
                      "Make sure to set refinement/njeans > 0 (typically 8-16).");
    // cells per current-sheet field-reversal scale L_B=|B|/|curlB| below which we refine.
    const auto curr_nsheet = pin->GetOrAddReal("refinement", "curr_nsheet", 4.0);
    PARTHENON_REQUIRE(curr_nsheet > 0.,
                      "Make sure to set refinement/curr_nsheet > 0 (typically 4-8).");
    // Density gate (code units): current-sheet refinement only applies to cells denser than this,
    // so it targets the collapsing region, NOT the diffuse turbulent envelope (which is full of
    // irrelevant numerical current sheets -> block runaway/OOM if ungated). Default 10 code =
    // ~2x the canonical BE-sphere t=0 peak (~5 code), so it activates only once gas is clearly
    // collapsing; tune per IC.
    const auto curr_rho_thresh =
        pin->GetOrAddReal("refinement", "curr_rho_thresh", 10.0);
    PARTHENON_REQUIRE(curr_rho_thresh >= 0.,
                      "Make sure to set refinement/curr_rho_thresh >= 0 (code density).");
    // Max refinement level the current-sheet criterion may drive (Jeans still refines deeper).
    // Bounds the deep-core block explosion; default 99 = effectively uncapped (prior behavior).
    const auto curr_max_level =
        pin->GetOrAddInteger("refinement", "curr_max_level", 99);
    pkg->AddParam<Real>("refinement/njeans", njeans);
    pkg->AddParam<Real>("refinement/curr_nsheet", curr_nsheet);
    pkg->AddParam<Real>("refinement/curr_rho_thresh", curr_rho_thresh);
    pkg->AddParam<int>("refinement/curr_max_level", curr_max_level);
  } else if (refine_str == "user") {
    pkg->CheckRefinementBlock = Hydro::ProblemCheckRefinementBlock;
  }

  if (ProblemInitPackageData != nullptr) {
    ProblemInitPackageData(pin, pkg.get());
  }

  return pkg;
}

template <Fluid fluid>
Real EstimateHyperbolicTimestep(MeshData<Real> *md) {
  // get to package via first block in Meshdata (which exists by construction)
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  const auto &cfl_hyp = hydro_pkg->Param<Real>("cfl");
  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});
  const auto &eos_ =
      hydro_pkg->Param<typename std::conditional<fluid == Fluid::euler, AdiabaticHydroEOS,
                                                 AdiabaticGLMMHDEOS>::type>("eos");

  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  Real min_dt_hyperbolic = std::numeric_limits<Real>::max();

  const auto ndim_ = prim_pack.GetNdim();
  Kokkos::parallel_reduce(
      "EstimateHyperbolicTimestep",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExecSpace(), {0, kb.s, jb.s, ib.s},
          {prim_pack.GetDim(5), kb.e + 1, jb.e + 1, ib.e + 1},
          {1, 1, 1, ib.e + 1 - ib.s}),
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, Real &min_dt) {
        const auto &prim = prim_pack(b);
        const auto &coords = prim_pack.GetCoords(b);
        // Need to reference variables here so that they are properly caught by
        // nvcc, which cannot determine captured variables only used within constexpr if.
        const auto &ndim = ndim_;
        const auto &eos = eos_;

        Real w[(NHYDRO)];
        w[IDN] = prim(IDN, k, j, i);
        w[IV1] = prim(IV1, k, j, i);
        w[IV2] = prim(IV2, k, j, i);
        w[IV3] = prim(IV3, k, j, i);
        w[IPR] = prim(IPR, k, j, i);
        Real lambda_max_x, lambda_max_y, lambda_max_z;
        if constexpr (fluid == Fluid::euler) {
          lambda_max_x = eos.SoundSpeed(w);
          lambda_max_y = lambda_max_x;
          lambda_max_z = lambda_max_x;

        } else if constexpr (fluid == Fluid::glmmhd) {
          lambda_max_x = eos.FastMagnetosonicSpeed(
              w[IDN], w[IPR], prim(IB1, k, j, i), prim(IB2, k, j, i), prim(IB3, k, j, i));
          if (ndim > 1) {
            lambda_max_y =
                eos.FastMagnetosonicSpeed(w[IDN], w[IPR], prim(IB2, k, j, i),
                                          prim(IB3, k, j, i), prim(IB1, k, j, i));
          }
          if (ndim > 2) {
            lambda_max_z =
                eos.FastMagnetosonicSpeed(w[IDN], w[IPR], prim(IB3, k, j, i),
                                          prim(IB1, k, j, i), prim(IB2, k, j, i));
          }
        } else {
          PARTHENON_FAIL("Unknown fluid in EstimateTimestep");
        }
        min_dt = fmin(min_dt, coords.Dxc<1>(k, j, i) / (fabs(w[IV1]) + lambda_max_x));
        if (ndim > 1) {
          min_dt = fmin(min_dt, coords.Dxc<2>(k, j, i) / (fabs(w[IV2]) + lambda_max_y));
        }
        if (ndim > 2) {
          min_dt = fmin(min_dt, coords.Dxc<3>(k, j, i) / (fabs(w[IV3]) + lambda_max_z));
        }
      },
      Kokkos::Min<Real>(min_dt_hyperbolic));

  // TODO(pgrete) THIS WORKAROUND IS NOT THREAD SAFE (though this will only become
  // relevant once parthenon uses host-multithreading in the driver).
  // We need to save the the hyperbolic part to recover it later as
  // the divergence cleaning speed is only limited in relation to the other
  // hyperbolic signal speeds and not by (potentially more restrictive) diffusive
  // processes.
  if constexpr (fluid == Fluid::glmmhd) {
    auto dt_hyp_pkg = hydro_pkg->Param<Real>("dt_hyp");
    if (cfl_hyp * min_dt_hyperbolic < dt_hyp_pkg) {
      hydro_pkg->UpdateParam("dt_hyp", cfl_hyp * min_dt_hyperbolic);
    }
  }
  return cfl_hyp * min_dt_hyperbolic;
}

// provide the routine that estimates a stable timestep for this package
template <Fluid fluid>
Real EstimateTimestep(MeshData<Real> *md) {
  // get to package via first block in Meshdata (which exists by construction)
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");
  auto min_dt = std::numeric_limits<Real>::max();
  auto dt_hyp = std::numeric_limits<Real>::max();

  const auto calc_dt_hyp = hydro_pkg->Param<bool>("calc_dt_hyp");
  if (calc_dt_hyp) {
    dt_hyp = EstimateHyperbolicTimestep<fluid>(md);
    min_dt = std::min(min_dt, dt_hyp);
  }

  const auto &enable_cooling = hydro_pkg->Param<Cooling>("enable_cooling");

  if (enable_cooling == Cooling::tabular) {
    const TabularCooling &tabular_cooling =
        hydro_pkg->Param<TabularCooling>("tabular_cooling");

    min_dt = std::min(min_dt, tabular_cooling.EstimateTimeStep(md));
  }

  auto dt_diff = std::numeric_limits<Real>::max();
  if (hydro_pkg->Param<DiffInt>("diffint") != DiffInt::none) {
    const auto diffint = hydro_pkg->Param<DiffInt>("diffint");
    // Strict dt constraint of the unsplit-in-rkl2 Hall part (mixed mode): the dispersive
    // whistler limit plus, if diffusion/hall_floor_integrator=unsplit, the floor's
    // parabolic limit. Kept OUT of dt_diff, which under rkl2 only sets the STS stage
    // count (and the optional ratio cap).
    auto dt_strict_hall = std::numeric_limits<Real>::max();
    if (hydro_pkg->Param<Conduction>("conduction") != Conduction::none) {
      dt_diff = std::min(dt_diff, EstimateConductionTimestep(md));
    }
    if (hydro_pkg->Param<Viscosity>("viscosity") != Viscosity::none) {
      dt_diff = std::min(dt_diff, EstimateViscosityTimestep(md));
    }
    // When every active non-ideal term uses the "ionization" coefficient, evaluate the
    // (expensive) Wardle conductivity tensor ONCE per cell and reduce the Ohmic/ambipolar/
    // Hall diffusive dt together. In the mixed rkl2+Hall mode the two-reducer FusedMixed
    // variant keeps the whistler limit separate from the parabolic aggregate (replacing
    // up to four full tensor sweeps per cycle with one); otherwise the single-reducer
    // fused estimator applies. Fall back to the exact per-term estimators when fusion is
    // off (diffusion/fused_nonideal_dt=false) or coefficients are not all "ionization".
    const bool mixed_hall_dt = (diffint == DiffInt::rkl2) &&
                               (hydro_pkg->Param<Hall>("hall") != Hall::none);
    if (hydro_pkg->Param<bool>("nonideal_dt_fused") && mixed_hall_dt) {
      const bool floor_rkl2 = hydro_pkg->Param<bool>("hall_floor_int_rkl2");
      Real dt_par, dt_strict;
      EstimateNonidealTimestepIonizationFusedMixed(md, !floor_rkl2, dt_par, dt_strict);
      dt_diff = std::min(dt_diff, dt_par);
      dt_strict_hall = std::min(dt_strict_hall, dt_strict);
    } else if (hydro_pkg->Param<bool>("nonideal_dt_fused")) {
      dt_diff = std::min(dt_diff, EstimateNonidealTimestepIonizationFused(md));
    } else {
      if (hydro_pkg->Param<Resistivity>("resistivity") != Resistivity::none) {
        dt_diff = std::min(dt_diff, EstimateResistivityTimestep(md));
      }
      if (hydro_pkg->Param<Ambipolar>("ambipolar") != Ambipolar::none) {
        dt_diff = std::min(dt_diff, EstimateAmbipolarTimestep(md));
      }
      if (hydro_pkg->Param<Hall>("hall") != Hall::none) {
        if (diffint == DiffInt::unsplit) {
          // whistler + floor both strict via the unsplit branch below
          dt_diff = std::min(dt_diff, EstimateHallTimestep(md, true, true));
        } else {
          // Mixed mode: whistler (+ floor if unsplit) strict; floor joins the parabolic
          // aggregate instead when it is RKL2-integrated.
          const bool floor_rkl2 = hydro_pkg->Param<bool>("hall_floor_int_rkl2");
          dt_strict_hall =
              std::min(dt_strict_hall, EstimateHallTimestep(md, true, !floor_rkl2));
          if (floor_rkl2) {
            dt_diff = std::min(dt_diff, EstimateHallTimestep(md, false, true));
          }
        }
      }
    }

    // For unsplit ingegration use strict limit
    if (diffint == DiffInt::unsplit) {
      min_dt = std::min(min_dt, dt_diff);
      // and for RKL2 integration use limit taking into account the maxium ratio
      // or not constrain limit further (which is why RKL2 is there in first place)
    } else if (diffint == DiffInt::rkl2) {
      min_dt = std::min(min_dt, dt_strict_hall); // no-op unless mixed rkl2+Hall mode
      const auto max_dt_ratio = hydro_pkg->Param<Real>("rkl2_max_dt_ratio");
      if (max_dt_ratio > 0.0 && dt_hyp / dt_diff > max_dt_ratio) {
        min_dt = std::min(min_dt, max_dt_ratio * dt_diff);
      }
    } else {
      PARTHENON_THROW("Looks like a a new diffusion integrator was implemented without "
                      "taking into accout timestep contstraints yet.");
    }
    auto dt_diff_param = hydro_pkg->Param<Real>("dt_diff");
    hydro_pkg->UpdateParam("dt_diff", std::min(dt_diff, dt_diff_param));
  }

  if (ProblemEstimateTimestep != nullptr) {
    min_dt = std::min(min_dt, ProblemEstimateTimestep(md));
  }

  // maximum user dt
  const auto max_dt = hydro_pkg->Param<Real>("max_dt");
  if (max_dt > 0.0) {
    min_dt = std::min(min_dt, max_dt);
  }

  return min_dt;
}

// Calculate fluxes using a tightly nested 3D loop over the entire block.
// Currently only used for testing the LLF Riemann solver used in first-order flux corr.
template <Fluid fluid>
TaskStatus CalculateFluxesTight(std::shared_ptr<MeshData<Real>> &md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  std::vector<parthenon::MetadataFlag> flags_ind({Metadata::Independent});
  auto cons_in = md->PackVariablesAndFluxes(flags_ind);
  auto pkg = pmb->packages.Get("Hydro");

  const auto &eos =
      pkg->Param<typename std::conditional<fluid == Fluid::euler, AdiabaticHydroEOS,
                                           AdiabaticGLMMHDEOS>::type>("eos");

  // Hyperbolic divergence cleaning speed for GLM MHD
  Real c_h = 0.0;
  if (fluid == Fluid::glmmhd) {
    c_h = pkg->Param<Real>("c_h");
  }
  // TODO(pgrete) fix scalar fluxes, too
  auto const &prim_in = md->PackVariables(std::vector<std::string>{"prim"});

  const int ndim = pmb->pmy_mesh->ndim;
  auto riemann = Riemann<fluid, RiemannSolver::llf>();
  // loop bounds are chosen so that all active fluxes are calculated
  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "DC LLF fluxes", parthenon::DevExecSpace(), 0,
      cons_in.GetDim(5) - 1, kb.s, kb.e + 1, jb.s, jb.e + 1, ib.s, ib.e + 1,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        auto &cons = cons_in(b);
        const auto &prim = prim_in(b);
        riemann.Solve(eos, k, j, i, IV1, prim, cons, c_h);
        if (ndim >= 2) {
          riemann.Solve(eos, k, j, i, IV2, prim, cons, c_h);
        }
        if (ndim >= 3) {
          riemann.Solve(eos, k, j, i, IV3, prim, cons, c_h);
        }
      });

  return TaskStatus::complete;
}

// Calculate fluxes using scratch pad memory, i.e., over cached pencils in i-dir.
template <Fluid fluid, Reconstruction recon, RiemannSolver rsolver>
TaskStatus CalculateFluxes(std::shared_ptr<MeshData<Real>> &md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
  int il, iu, jl, ju, kl, ku;
  jl = jb.s, ju = jb.e, kl = kb.s, ku = kb.e;
  // TODO(pgrete): are these looop limits are likely too large for 2nd order
  if (pmb->block_size.nx(X2DIR) > 1) {
    if (pmb->block_size.nx(X3DIR) == 1) // 2D
      jl = jb.s - 1, ju = jb.e + 1, kl = kb.s, ku = kb.e;
    else // 3D
      jl = jb.s - 1, ju = jb.e + 1, kl = kb.s - 1, ku = kb.e + 1;
  }

  std::vector<parthenon::MetadataFlag> flags_ind({Metadata::Independent});
  auto cons_in = md->PackVariablesAndFluxes(flags_ind);
  auto pkg = pmb->packages.Get("Hydro");
  const auto nhydro = pkg->Param<int>("nhydro");
  const auto nscalars = pkg->Param<int>("nscalars");

  const auto &eos =
      pkg->Param<typename std::conditional<fluid == Fluid::euler, AdiabaticHydroEOS,
                                           AdiabaticGLMMHDEOS>::type>("eos");

  auto num_scratch_vars = nhydro + nscalars;

  // Hyperbolic divergence cleaning speed for GLM MHD
  Real c_h = 0.0;
  if (fluid == Fluid::glmmhd) {
    c_h = pkg->Param<Real>("c_h");
  }

  auto const &prim_in = md->PackVariables(std::vector<std::string>{"prim"});

  const int scratch_level =
      pkg->Param<int>("scratch_level"); // 0 is actual scratch (tiny); 1 is HBM
  const int nx1 = pmb->cellbounds.ncellsi(IndexDomain::entire);

  size_t scratch_size_in_bytes =
      parthenon::ScratchPad2D<Real>::shmem_size(num_scratch_vars, nx1) * 2;

  auto riemann = Riemann<fluid, rsolver>();

  parthenon::par_for_outer(
      DEFAULT_OUTER_LOOP_PATTERN, "x1 flux", DevExecSpace(), scratch_size_in_bytes,
      scratch_level, 0, cons_in.GetDim(5) - 1, kl, ku, jl, ju,
      KOKKOS_LAMBDA(parthenon::team_mbr_t member, const int b, const int k, const int j) {
        const auto &prim = prim_in(b);
        auto &cons = cons_in(b);
        parthenon::ScratchPad2D<Real> wl(member.team_scratch(scratch_level),
                                         num_scratch_vars, nx1);
        parthenon::ScratchPad2D<Real> wr(member.team_scratch(scratch_level),
                                         num_scratch_vars, nx1);
        // get reconstructed state on faces
        Reconstruct<recon, X1DIR>(member, k, j, ib.s - 1, ib.e + 1, prim, wl, wr);
        // Sync all threads in the team so that scratch memory is consistent
        member.team_barrier();

        riemann.Solve(member, k, j, ib.s, ib.e + 1, IV1, wl, wr, cons, eos, c_h);
        member.team_barrier();

        // Passive scalar fluxes
        for (auto n = nhydro; n < nhydro + nscalars; ++n) {
          parthenon::par_for_inner(member, ib.s, ib.e + 1, [&](const int i) {
            if (cons.flux(IV1, IDN, k, j, i) >= 0.0) {
              cons.flux(IV1, n, k, j, i) = cons.flux(IV1, IDN, k, j, i) * wl(n, i);
            } else {
              cons.flux(IV1, n, k, j, i) = cons.flux(IV1, IDN, k, j, i) * wr(n, i);
            }
          });
        }
      });

  //--------------------------------------------------------------------------------------
  // j-direction
  if (pmb->pmy_mesh->ndim >= 2) {
    scratch_size_in_bytes =
        parthenon::ScratchPad2D<Real>::shmem_size(num_scratch_vars, nx1) * 3;
    // set the loop limits
    il = ib.s - 1, iu = ib.e + 1, kl = kb.s, ku = kb.e;
    if (pmb->block_size.nx(X3DIR) == 1) // 2D
      kl = kb.s, ku = kb.e;
    else // 3D
      kl = kb.s - 1, ku = kb.e + 1;

    parthenon::par_for_outer(
        DEFAULT_OUTER_LOOP_PATTERN, "x2 flux", DevExecSpace(), scratch_size_in_bytes,
        scratch_level, 0, cons_in.GetDim(5) - 1, kl, ku,
        KOKKOS_LAMBDA(parthenon::team_mbr_t member, const int b, const int k) {
          const auto &prim = prim_in(b);
          auto &cons = cons_in(b);
          parthenon::ScratchPad2D<Real> wl(member.team_scratch(scratch_level),
                                           num_scratch_vars, nx1);
          parthenon::ScratchPad2D<Real> wr(member.team_scratch(scratch_level),
                                           num_scratch_vars, nx1);
          parthenon::ScratchPad2D<Real> wlb(member.team_scratch(scratch_level),
                                            num_scratch_vars, nx1);
          for (int j = jb.s - 1; j <= jb.e + 1; ++j) {
            // reconstruct L/R states at j
            Reconstruct<recon, X2DIR>(member, k, j, il, iu, prim, wlb, wr);
            // Sync all threads in the team so that scratch memory is consistent
            member.team_barrier();

            if (j > jb.s - 1) {
              riemann.Solve(member, k, j, il, iu, IV2, wl, wr, cons, eos, c_h);
              member.team_barrier();

              // Passive scalar fluxes
              for (auto n = nhydro; n < nhydro + nscalars; ++n) {
                parthenon::par_for_inner(member, il, iu, [&](const int i) {
                  if (cons.flux(IV2, IDN, k, j, i) >= 0.0) {
                    cons.flux(IV2, n, k, j, i) = cons.flux(IV2, IDN, k, j, i) * wl(n, i);
                  } else {
                    cons.flux(IV2, n, k, j, i) = cons.flux(IV2, IDN, k, j, i) * wr(n, i);
                  }
                });
              }
              member.team_barrier();
            }

            // swap the arrays for the next step
            auto *tmp = wl.data();
            wl.assign_data(wlb.data());
            wlb.assign_data(tmp);
          }
        });
  }
  //--------------------------------------------------------------------------------------
  // k-direction
  if (pmb->pmy_mesh->ndim >= 3) {
    // set the loop limits
    il = ib.s - 1, iu = ib.e + 1, jl = jb.s - 1, ju = jb.e + 1;

    parthenon::par_for_outer(
        DEFAULT_OUTER_LOOP_PATTERN, "x3 flux", DevExecSpace(), scratch_size_in_bytes,
        scratch_level, 0, cons_in.GetDim(5) - 1, jl, ju,
        KOKKOS_LAMBDA(parthenon::team_mbr_t member, const int b, const int j) {
          const auto &prim = prim_in(b);
          auto &cons = cons_in(b);
          parthenon::ScratchPad2D<Real> wl(member.team_scratch(scratch_level),
                                           num_scratch_vars, nx1);
          parthenon::ScratchPad2D<Real> wr(member.team_scratch(scratch_level),
                                           num_scratch_vars, nx1);
          parthenon::ScratchPad2D<Real> wlb(member.team_scratch(scratch_level),
                                            num_scratch_vars, nx1);
          for (int k = kb.s - 1; k <= kb.e + 1; ++k) {
            // reconstruct L/R states at j
            Reconstruct<recon, X3DIR>(member, k, j, il, iu, prim, wlb, wr);
            // Sync all threads in the team so that scratch memory is consistent
            member.team_barrier();

            if (k > kb.s - 1) {
              riemann.Solve(member, k, j, il, iu, IV3, wl, wr, cons, eos, c_h);
              member.team_barrier();

              // Passive scalar fluxes
              for (auto n = nhydro; n < nhydro + nscalars; ++n) {
                parthenon::par_for_inner(member, il, iu, [&](const int i) {
                  if (cons.flux(IV3, IDN, k, j, i) >= 0.0) {
                    cons.flux(IV3, n, k, j, i) = cons.flux(IV3, IDN, k, j, i) * wl(n, i);
                  } else {
                    cons.flux(IV3, n, k, j, i) = cons.flux(IV3, IDN, k, j, i) * wr(n, i);
                  }
                });
              }
              member.team_barrier();
            }
            // swap the arrays for the next step
            auto *tmp = wl.data();
            wl.assign_data(wlb.data());
            wlb.assign_data(tmp);
          }
        });
  }

  const auto &diffint = pkg->Param<DiffInt>("diffint");
  if (diffint == DiffInt::unsplit) {
    CalcDiffFluxes(pkg.get(), md.get(), DiffTermSet::all);
  } else if (diffint == DiffInt::rkl2 && pkg->Param<Hall>("hall") != Hall::none) {
    // Mixed mode: the dispersive Hall EMF (plus the Ohmic floor unless
    // hall_floor_integrator=rkl2) rides on the hyperbolic fluxes every stage under the
    // strict whistler dt; the parabolic terms are applied by the RKL2 tasks.
    CalcDiffFluxes(pkg.get(), md.get(), DiffTermSet::dispersive);
  }

  return TaskStatus::complete;
}

// Apply first order flux correction, i.e., use first order reconstruction and a
// diffusive LLF Riemann solver if a negative density or energy density is expected.
// The current implementation is computationally not the most efficient one, but works
// for all standard integrators (rk1, rk2, rk3, and vl) and with and without AMR.
// In principle, without AMR one could directly use the results from the actual
// flux divergence call.
// However, with AMR (and coarse/fine flux correction) we need to correct the local
// fluxes first before calling coarse/fine flux correction.
// In addition, it may be enough to call first order flux correction once at the
// final stage (rather than at every stage as right row).
// However, this'd require an additional register to store the initial state and
// we should first evaluate where the tradeoff between extra computational costs
// (multiple calls) versus extra memory usage is.
template <Fluid fluid>
TaskStatus FirstOrderFluxCorrect(MeshData<Real> *u0_data, MeshData<Real> *u1_data,
                                 const Real gam0_, const Real gam1_,
                                 const Real beta_dt_) {
  // Work around for CUDA <=11.6
  const Real gam0 = gam0_;
  const Real gam1 = gam1_;
  const Real beta_dt = beta_dt_;

  auto pmb = u0_data->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  // Pack "cons" by NAME rather than by the {Independent} flag. Under Constrained Transport
  // the {Independent} flag also matches the FACE field Bf (WithFluxes) whose EDGE fluxes
  // (E1/E2/E3) have no X1/X2/X3 cell-flux slots -- indexing them in this cell-centered
  // kernel segfaults (cudaErrorIllegalAddress). FOFC only needs the cell-centered hydro cons
  // (rho, mom, E, cell-B, psi); CT evolves the staggered Bf via its own EMF, so the
  // cell-B-flux redo here is inert (the cell-centered B is re-projected from Bf each stage).
  // Mirrors the cons-by-name pack the STS/driver use for the same CT reason (hydro_driver.cpp).
  const std::vector<std::string> cons_name{"cons"};
  auto u0_cons_pack = u0_data->PackVariablesAndFluxes(cons_name, cons_name);
  auto const &u0_prim_pack = u0_data->PackVariables(std::vector<std::string>{"prim"});
  auto u1_cons_pack = u1_data->PackVariablesAndFluxes(cons_name, cons_name);
  auto pkg = pmb->packages.Get("Hydro");

  const auto &eos =
      pkg->Param<typename std::conditional<fluid == Fluid::euler, AdiabaticHydroEOS,
                                           AdiabaticGLMMHDEOS>::type>("eos");

  // Hyperbolic divergence cleaning speed for GLM MHD
  Real c_h = 0.0;
  if (fluid == Fluid::glmmhd) {
    c_h = pkg->Param<Real>("c_h");
  }

  const int ndim = pmb->pmy_mesh->ndim;

  constexpr auto NVAR = GetNVars<fluid>();

  auto riemann = Riemann<fluid, RiemannSolver::llf>();

  std::int64_t num_corrected, num_need_floor;
  // Potentially need multiple attempts as flux correction corrects 6 (in 3D) fluxes
  // of a single cell at the same time. So the neighboring cells need to be rechecked with
  // the corrected fluxes as the corrected fluxes in one cell may result in the need to
  // correct all the fluxes of an originally "good" neighboring cell.
  size_t num_attempts = 0;
  do {
    num_corrected = 0;

    Kokkos::parallel_reduce(
        "FirstOrderFluxCorrect",
        Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
            DevExecSpace(), {0, kb.s, jb.s, ib.s},
            {u0_cons_pack.GetDim(5), kb.e + 1, jb.e + 1, ib.e + 1},
            {1, 1, 1, ib.e + 1 - ib.s}),
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i,
                      std::int64_t &lnum_corrected, std::int64_t &lnum_need_floor) {
          const auto &coords = u0_cons_pack.GetCoords(b);
          const auto &u0_prim = u0_prim_pack(b);
          auto &u0_cons = u0_cons_pack(b);

          // In principle, the u_cons.fluxes could be updated in parallel by a
          // different thread resulting in a race conditon here. However, if the
          // fluxes of a cell have been updated (anywhere) then the entire kernel will
          // be called again anyway, and, at that point the already fixed
          // u0_cons.fluxes will automaticlly be used here.
          Real new_cons[NVAR];
          for (auto v = 0; v < NVAR; v++) {
            new_cons[v] =
                gam0 * u0_cons(v, k, j, i) + gam1 * u1_cons_pack(b, v, k, j, i) +
                beta_dt *
                    parthenon::Update::FluxDivHelper(v, k, j, i, ndim, coords, u0_cons);
          }

          // no need to include gamma - 1 as we only care for negative values
          auto new_p =
              new_cons[IEN] -
              0.5 * (SQR(new_cons[IM1]) + SQR(new_cons[IM2]) + SQR(new_cons[IM3])) /
                  new_cons[IDN];
          if constexpr (fluid == Fluid::glmmhd) {
            new_p -= 0.5 * (SQR(new_cons[IB1]) + SQR(new_cons[IB2]) + SQR(new_cons[IB3]));
          }
          // no correction required
          if (new_cons[IDN] > 0.0 && new_p > 0.0) {
            return;
          }
          // if already tried 3 times and only pressure is negative, then we'll rely
          // on the pressure floor during ConsToPrim conversion
          if (num_attempts > 2 && new_cons[IDN] > 0.0 && new_p < 0.0) {
            lnum_need_floor += 1;
            return;
          }
          // In principle, there could be a racecondion as this loop goes over all
          // k,j,i and we updating the i+1 flux here. However, the results are
          // idential because u0_prim is never updated in this kernel so we don't
          // worry about it.
          // TODO(pgrete) as we need to keep the function signature idential for now
          // (due to Cuda compiler bug) we could potentially template these function
          // and get rid of the `if constexpr`
          riemann.Solve(eos, k, j, i, IV1, u0_prim, u0_cons, c_h);
          riemann.Solve(eos, k, j, i + 1, IV1, u0_prim, u0_cons, c_h);

          if (ndim >= 2) {
            riemann.Solve(eos, k, j, i, IV2, u0_prim, u0_cons, c_h);
            riemann.Solve(eos, k, j + 1, i, IV2, u0_prim, u0_cons, c_h);
          }
          if (ndim >= 3) {
            riemann.Solve(eos, k, j, i, IV3, u0_prim, u0_cons, c_h);
            riemann.Solve(eos, k + 1, j, i, IV3, u0_prim, u0_cons, c_h);
          }
          lnum_corrected += 1;
        },
        Kokkos::Sum<std::int64_t>(num_corrected),
        Kokkos::Sum<std::int64_t>(num_need_floor));
    // TODO(pgrete) make this optional and global (potentially store values in Params)
    // std::cout << "[" << parthenon::Globals::my_rank << "] Attempt: " <<
    // num_attempts
    //           << " Corrected (center): " << num_corrected
    //           << " Failed (will rely on floor): " << num_need_floor << std::endl;
    num_attempts += 1;
  } while (num_corrected > 0 && num_attempts < 4);

  return TaskStatus::complete;
}

} // namespace Hydro

