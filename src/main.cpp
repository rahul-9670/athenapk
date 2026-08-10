// AthenaPK - a performance portable block structured AMR MHD code
// Copyright (c) 2020-2021, Athena Parthenon Collaboration. All rights reserved.
// Licensed under the 3-Clause License (the "LICENSE");

#include <sstream>

// Parthenon headers
#include "bvals/boundary_conditions_generic.hpp"
#include "defs.hpp"
#include "globals.hpp"
#include "parthenon_manager.hpp"

// AthenaPK headers
#include "bvals/boundary_conditions_apk.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_driver.hpp"
#include "main.hpp"

#include "pgen/pgen.hpp"
#include "tracers/tracers.hpp"
// Initialize defaults for package specific callback functions
namespace Hydro {
InitPackageDataFun_t ProblemInitPackageData = nullptr;
SourceFun_t ProblemSourceFirstOrder = nullptr;
SourceFun_t ProblemSourceStrangSplit = nullptr;
SourceFun_t ProblemSourceUnsplit = nullptr;
EstimateTimestepFun_t ProblemEstimateTimestep = nullptr;
std::function<AmrTag(MeshBlockData<Real> *mbd)> ProblemCheckRefinementBlock = nullptr;
} // namespace Hydro

namespace Tracers {
InitPackageDataFun_t ProblemInitTracerData = nullptr;
SeedInitialFun_t ProblemSeedInitialTracers = nullptr;
FillTracersFun_t ProblemFillTracers = nullptr;
} // namespace Tracers

int main(int argc, char *argv[]) {
  using parthenon::ParthenonManager;
  using parthenon::ParthenonStatus;
  ParthenonManager pman;

  // call ParthenonInit to initialize MPI and Kokkos, parse the input deck, and set up
  auto manager_status = pman.ParthenonInitEnv(argc, argv);
  if (manager_status == ParthenonStatus::complete) {
    pman.ParthenonFinalize();
    return 0;
  }
  if (manager_status == ParthenonStatus::error) {
    pman.ParthenonFinalize();
    return 1;
  }
  // Now that ParthenonInit has been called and setup succeeded, the code can now
  // make use of MPI and Kokkos

  // Redefine defaults
  pman.app_input->ProcessPackages = Hydro::ProcessPackages;
  pman.app_input->PreStepMeshUserWorkInLoop = Hydro::PreStepMeshUserWorkInLoop;
  const auto problem = pman.pinput->GetOrAddString("job", "problem_id", "unset");

  // Flagship build: `collapse_be` is the only problem generator compiled in.
  // The other 23 (upstream demos + the validation problems) live in
  // /beegfs/u/bbg6470/validation/code/pgen/ as of 2026-08-10, and in git at tag
  // `validation-complete-2026-08-10`.  Re-enabling one means restoring its .cpp,
  // its namespace declaration in pgen/pgen.hpp, its entry in pgen/CMakeLists.txt,
  // and an `else if` branch here.
  if (problem == "collapse_be") {
    pman.app_input->ProblemGenerator = collapse_be::ProblemGenerator;
    // WP-13: restart-safe Param registration. ProblemGenerator is NOT called on restart, so
    // this hook (run from Hydro::Initialize on both fresh start and resume) is what keeps
    // collapse_be_rhocrit alive -- and with it the ApplyBarotropicCooling task, which carries
    // the outside-sphere momentum BC. Without it every restart silently dropped that task.
    Hydro::ProblemInitPackageData = collapse_be::ProblemInitPackageData;
  } else {
    std::stringstream msg;
    msg << "Problem ID '" << problem
        << "' is not available in this build. This is the flagship build, which "
           "compiles only 'collapse_be'; the other problem generators are in "
           "/beegfs/u/bbg6470/validation/code/pgen/ (git tag "
           "validation-complete-2026-08-10).";
    PARTHENON_THROW(msg);
  }

  // The "reflecting" BC was removed on 2026-08-10 with the pure-hydro path. Its
  // implementation opened with PARTHENON_REQUIRE_THROWS(fluid == Fluid::euler,
  // "Reflecting boundary conditions for MHD need special treatment.") -- so in an MHD-only
  // build it could only ever abort. Production uses `diode` (registered below) or
  // Parthenon's built-in `outflow`/`periodic`.
  using BF = parthenon::BoundaryFace;
  using parthenon::BoundaryFunction::BCSide;

  // VALIDATION B1: inflow-suppressing outflow ("diode"). Parthenon's `outflow` is a plain
  // zero-gradient copy and admits inflow -- measured as ~-178 per face in WP-6's cons-Mout,
  // the source of the production mass rise. Registering it as a SEPARATE named BC leaves
  // `outflow` untouched, so every existing deck stays bit-identical; switching production to
  // it is result-changing and therefore an explicit, deliberate deck edit.
  const std::string DIODE = "diode";
  using Hydro::BoundaryFunction::DiodeBC;
  pman.app_input->RegisterBoundaryCondition(BF::inner_x1, DIODE,
                                            DiodeBC<X1DIR, BCSide::Inner>);
  pman.app_input->RegisterBoundaryCondition(BF::outer_x1, DIODE,
                                            DiodeBC<X1DIR, BCSide::Outer>);
  pman.app_input->RegisterBoundaryCondition(BF::inner_x2, DIODE,
                                            DiodeBC<X2DIR, BCSide::Inner>);
  pman.app_input->RegisterBoundaryCondition(BF::outer_x2, DIODE,
                                            DiodeBC<X2DIR, BCSide::Outer>);
  pman.app_input->RegisterBoundaryCondition(BF::inner_x3, DIODE,
                                            DiodeBC<X3DIR, BCSide::Inner>);
  pman.app_input->RegisterBoundaryCondition(BF::outer_x3, DIODE,
                                            DiodeBC<X3DIR, BCSide::Outer>);

  pman.ParthenonInitPackagesAndMesh();

  if (parthenon::Globals::my_rank == 0) {
    std::cout << "Starting up hydro driver" << std::endl;
  }

  {
    Hydro::HydroDriver driver(pman.pinput.get(), pman.app_input.get(), pman.pmesh.get());

    driver.Execute();
  }

  pman.ParthenonFinalize();


  return (0);
}
