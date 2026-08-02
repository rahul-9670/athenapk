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

  if (problem == "linear_wave") {
    pman.app_input->InitUserMeshData = linear_wave::InitUserMeshData;
    pman.app_input->ProblemGenerator = linear_wave::ProblemGenerator;
    pman.app_input->UserWorkAfterLoop = linear_wave::UserWorkAfterLoop;
  } else if (problem == "linear_wave_mhd") {
    pman.app_input->InitUserMeshData = linear_wave_mhd::InitUserMeshData;
    pman.app_input->ProblemGenerator = linear_wave_mhd::ProblemGenerator;
    pman.app_input->UserWorkAfterLoop = linear_wave_mhd::UserWorkAfterLoop;
    Hydro::ProblemInitPackageData = linear_wave_mhd::ProblemInitPackageData;
  } else if (problem == "jeans") {
    pman.app_input->ProblemGenerator = jeans::ProblemGenerator;
  } else if (problem == "polytrope") {
    pman.app_input->ProblemGenerator = polytrope::ProblemGenerator;
  } else if (problem == "collapse_be") {
    pman.app_input->ProblemGenerator = collapse_be::ProblemGenerator;
    // WP-13: restart-safe Param registration. ProblemGenerator is NOT called on restart, so
    // this hook (run from Hydro::Initialize on both fresh start and resume) is what keeps
    // collapse_be_rhocrit alive -- and with it the ApplyBarotropicCooling task, which carries
    // the outside-sphere momentum BC. Without it every restart silently dropped that task.
    Hydro::ProblemInitPackageData = collapse_be::ProblemInitPackageData;
  } else if (problem == "poisson_test") {
    pman.app_input->ProblemGenerator = poisson_test::ProblemGenerator;
  } else if (problem == "rad_pulse") {
    pman.app_input->ProblemGenerator = rad_pulse::ProblemGenerator;
  } else if (problem == "rad_shadow") {
    pman.app_input->InitUserMeshData = rad_shadow::InitUserMeshData;
    pman.app_input->ProblemGenerator = rad_shadow::ProblemGenerator;
    pman.app_input->RegisterBoundaryCondition(parthenon::BoundaryFace::inner_x1,
                                              "rad_shadow_beam_x1",
                                              rad_shadow::InflowBeamX1);
  } else if (problem == "rad_shock") {
    pman.app_input->ProblemGenerator = rad_shock::ProblemGenerator;
  } else if (problem == "cpaw") {
    pman.app_input->InitUserMeshData = cpaw::InitUserMeshData;
    pman.app_input->ProblemGenerator = cpaw::ProblemGenerator;
    pman.app_input->UserWorkAfterLoop = cpaw::UserWorkAfterLoop;
  } else if (problem == "cloud") {
    pman.app_input->InitUserMeshData = cloud::InitUserMeshData;
    pman.app_input->ProblemGenerator = cloud::ProblemGenerator;
    pman.app_input->RegisterBoundaryCondition(parthenon::BoundaryFace::inner_x2,
                                              "cloud_inflow_x2", cloud::InflowWindX2);
    Hydro::ProblemCheckRefinementBlock = cloud::ProblemCheckRefinementBlock;
  } else if (problem == "blast") {
    pman.app_input->InitUserMeshData = blast::InitUserMeshData;
    pman.app_input->ProblemGenerator = blast::ProblemGenerator;
    pman.app_input->UserWorkAfterLoop = blast::UserWorkAfterLoop;
  } else if (problem == "advection") {
    pman.app_input->InitUserMeshData = advection::InitUserMeshData;
    pman.app_input->ProblemGenerator = advection::ProblemGenerator;
  } else if (problem == "orszag_tang") {
    pman.app_input->ProblemGenerator = orszag_tang::ProblemGenerator;
  } else if (problem == "diffusion") {
    pman.app_input->ProblemGenerator = diffusion::ProblemGenerator;
    pman.app_input->UserWorkAfterLoop = diffusion::UserWorkAfterLoop;
  } else if (problem == "cshock") {
    pman.app_input->ProblemGenerator = cshock::ProblemGenerator;
  } else if (problem == "field_loop") {
    pman.app_input->ProblemGenerator = field_loop::ProblemGenerator;
    Hydro::ProblemInitPackageData = field_loop::ProblemInitPackageData;
  } else if (problem == "kh") {
    pman.app_input->MeshProblemGenerator = kh::ProblemGenerator;
  } else if (problem == "lw_implode") {
    pman.app_input->ProblemGenerator = lw_implode::ProblemGenerator;
  } else if (problem == "rand_blast") {
    pman.app_input->ProblemGenerator = rand_blast::ProblemGenerator;
    Hydro::ProblemInitPackageData = rand_blast::ProblemInitPackageData;
    Hydro::ProblemSourceFirstOrder = rand_blast::RandomBlasts;
  } else if (problem == "cluster") {
    pman.app_input->MeshProblemGenerator = cluster::ProblemGenerator;
    pman.app_input->MeshBlockUserWorkBeforeOutput = cluster::UserWorkBeforeOutput;
    Hydro::ProblemInitPackageData = cluster::ProblemInitPackageData;
    Hydro::ProblemSourceUnsplit = cluster::ClusterUnsplitSrcTerm;
    Hydro::ProblemSourceFirstOrder = cluster::ClusterSplitSrcTerm;
    Hydro::ProblemEstimateTimestep = cluster::ClusterEstimateTimestep;
  } else if (problem == "shattering") {
    pman.app_input->InitUserMeshData = shattering::InitUserMeshData;
    pman.app_input->ProblemGenerator = shattering::ProblemGenerator;
  } else if (problem == "sod") {
    pman.app_input->ProblemGenerator = sod::ProblemGenerator;
  } else if (problem == "turbulence") {
    pman.app_input->MeshProblemGenerator = turbulence::ProblemGenerator;
    Hydro::ProblemInitPackageData = turbulence::ProblemInitPackageData;
    Tracers::ProblemInitTracerData = turbulence::ProblemInitTracerData;
    Tracers::ProblemFillTracers = turbulence::ProblemFillTracers;
    Hydro::ProblemSourceFirstOrder = turbulence::Driving;
    pman.app_input->InitMeshBlockUserData = turbulence::SetPhases;
    pman.app_input->MeshBlockUserWorkBeforeOutput = turbulence::UserWorkBeforeOutput;
  } else {
    // parthenon throw error message for the invalid problem
    std::stringstream msg;
    msg << "Problem ID '" << problem << "' is not implemented yet.";
    PARTHENON_THROW(msg);
  }

  const std::string REFLECTING = "reflecting";
  using BF = parthenon::BoundaryFace;
  using Hydro::BoundaryFunction::ReflectBC;
  using parthenon::BoundaryFunction::BCSide;
  pman.app_input->RegisterBoundaryCondition(BF::inner_x1, REFLECTING,
                                            ReflectBC<X1DIR, BCSide::Inner>);
  pman.app_input->RegisterBoundaryCondition(BF::outer_x1, REFLECTING,
                                            ReflectBC<X1DIR, BCSide::Outer>);
  pman.app_input->RegisterBoundaryCondition(BF::inner_x2, REFLECTING,
                                            ReflectBC<X2DIR, BCSide::Inner>);
  pman.app_input->RegisterBoundaryCondition(BF::outer_x2, REFLECTING,
                                            ReflectBC<X2DIR, BCSide::Outer>);
  pman.app_input->RegisterBoundaryCondition(BF::inner_x3, REFLECTING,
                                            ReflectBC<X3DIR, BCSide::Inner>);
  pman.app_input->RegisterBoundaryCondition(BF::outer_x3, REFLECTING,
                                            ReflectBC<X3DIR, BCSide::Outer>);

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
