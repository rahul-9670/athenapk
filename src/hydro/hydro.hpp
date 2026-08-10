#ifndef HYDRO_HYDRO_HPP_
#define HYDRO_HYDRO_HPP_
//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2020-2021, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

// Parthenon headers
#include <parthenon/package.hpp>


using namespace parthenon::package::prelude;

namespace Hydro {

parthenon::Packages_t ProcessPackages(std::unique_ptr<ParameterInput> &pin);
void PreStepMeshUserWorkInLoop(Mesh *pmesh, ParameterInput *pin, parthenon::SimTime &tm);
std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin);

template <Fluid fluid>
Real EstimateTimestep(MeshData<Real> *md);

using parthenon::SimTime;
//! B1-remainder (2026-08-03): clamp the DOMAIN-BOUNDARY face fluxes so no gas can enter through
//! a nominally-outflow face. Runs after SetFluxCorrections and before the flux divergence.
//!
//! WHY THE `diode` BC IS NOT ENOUGH. `diode` zeroes the face-normal momentum in the GHOST cell
//! when it points inward. Measured, that removes ~51 % of the spurious mass gain (0.212 -> 0.104 %
//! at 128^3, 0.106 -> 0.052 % at 256^3, and independently 53 % on the 32^3 smoke deck via
//! cons-Mout-solver). It cannot remove the rest, because zeroing a ghost VELOCITY does not stop
//! the Riemann solver returning an inward mass flux when the INTERIOR cell is being accelerated
//! inward by gravity -- the ghost still carries density, and that density is advected in.
//! Clamping the flux itself is the only place the remaining inflow can be caught.
//!
//! SCOPE, deliberately narrow. Only the HYDRO conserved fluxes (IDN, IM1..3, IEN) and the passive
//! scalars are zeroed, and only on faces that lie on the physical domain boundary, and only when
//! the mass flux points INTO the domain. The MAGNETIC fluxes are left alone: under GLM
//! cons.flux(IB*) is the induction flux, and zeroing it at a boundary face would inject a
//! div(B) error rather than remove a mass error; under CT the field is advanced from edge EMFs
//! and these fluxes do not carry it at all. B1 is a mass/momentum/energy problem and the fix is
//! kept to those equations.
//!
//! Gated on `hydro/boundary_flux_clamp` (default false => the task is not even added, so the
//! OFF state is bit-identical).
TaskStatus ClampBoundaryFluxes(MeshData<Real> *md);

TaskStatus AddUnsplitSources(MeshData<Real> *md, const SimTime &tm, const Real beta_dt);
TaskStatus AddSplitSourcesFirstOrder(MeshData<Real> *md, const SimTime &tm);
TaskStatus AddSplitSourcesStrang(MeshData<Real> *md, const SimTime &tm);

using SourceFun_t =
    std::function<void(MeshData<Real> *md, const SimTime &tm, const Real dt)>;
using EstimateTimestepFun_t = std::function<Real(MeshData<Real> *md)>;

extern SourceFun_t ProblemSourceFirstOrder;
extern SourceFun_t ProblemSourceUnsplit;
extern SourceFun_t ProblemSourceStrangSplit;
extern EstimateTimestepFun_t ProblemEstimateTimestep;
extern InitPackageDataFun_t ProblemInitPackageData;
extern std::function<AmrTag(MeshBlockData<Real> *mbd)> ProblemCheckRefinementBlock;

template <Fluid fluid>
TaskStatus CalculateFluxesTight(std::shared_ptr<MeshData<Real>> &md);
template <Fluid fluid, Reconstruction recon, RiemannSolver rsolver>
TaskStatus CalculateFluxes(std::shared_ptr<MeshData<Real>> &md);
using FluxFun_t =
    decltype(CalculateFluxes<Fluid::glmmhd, Reconstruction::dc, RiemannSolver::hlle>);

template <Fluid fluid>
TaskStatus FirstOrderFluxCorrect(MeshData<Real> *u0_data, MeshData<Real> *u1_data,
                                 const Real gam0, const Real gam1, const Real beta_dt);
using FirstOrderFluxCorrectFun_t = decltype(FirstOrderFluxCorrect<Fluid::glmmhd>);

using FluxFunKey_t = std::tuple<Fluid, Reconstruction, RiemannSolver>;

// Add flux function pointer to map containing all compiled in flux functions
template <Fluid fluid, Reconstruction recon, RiemannSolver rsolver>
void add_flux_fun(std::map<FluxFunKey_t, FluxFun_t *> &flux_functions) {
  flux_functions[std::make_tuple(fluid, recon, rsolver)] =
      Hydro::CalculateFluxes<fluid, recon, rsolver>;
}

// Get number of "fluid" variable used
template <Fluid fluid>
constexpr size_t GetNVars();

template <>
constexpr size_t GetNVars<Fluid::glmmhd>() {
  return 9; // above plus B_x, B_y, B_z, psi
}

} // namespace Hydro

#endif // HYDRO_HYDRO_HPP_
