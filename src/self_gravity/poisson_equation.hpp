//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Poisson equation for self-gravity solver.
// Adapted from parthenon/example/poisson_gmg (C) Triad National Security, LLC.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
#ifndef SELF_GRAVITY_POISSON_EQUATION_HPP_
#define SELF_GRAVITY_POISSON_EQUATION_HPP_

#include <memory>
#include <string>
#include <utility>

#include <coordinates/coordinates.hpp>
#include <kokkos_abstraction.hpp>
#include <parthenon/package.hpp>

namespace SelfGravity {

constexpr parthenon::TopologicalElement te = parthenon::TopologicalElement::CC;

// Implements A.x = y and diag(A) for A = -grad^2 in conservative flux form.
// Templated only on var_t (the solution variable type). D = 1 is hard-coded
// (no variable-coefficient diffusion); alpha = 0 is hard-coded (no Helmholtz shift).
// This is a stripped Parthenon poisson_gmg PoissonEquation.
template <class var_t>
class PoissonEquation {
 public:
  using IndependentVars = parthenon::TypeList<var_t>;

  // Flux correction: on for AMR base grid, off inside MG two-level composite grids.
  // Matches the Artemis gating exactly.
  PoissonEquation(parthenon::ParameterInput *pin, const std::string &label) {}

  // y = A.x
  parthenon::TaskID Ax(parthenon::TaskList &tl, parthenon::TaskID depends_on,
                       std::shared_ptr<parthenon::MeshData<Real>> &md_mat,
                       std::shared_ptr<parthenon::MeshData<Real>> &md_in,
                       std::shared_ptr<parthenon::MeshData<Real>> &md_out) {
    auto flux_res = tl.AddTask(depends_on, CalculateFluxes, md_in);
    // Flux correction at coarse-fine boundaries, but NOT inside MG composite grids
    // (MG has its own coarse-fine handling).
    if (!(md_mat->grid.type == parthenon::GridType::two_level_composite)) {
      auto start_flxcor =
          tl.AddTask(flux_res, parthenon::StartReceiveFluxCorrections, md_in);
      auto send_flxcor =
          tl.AddTask(flux_res, parthenon::LoadAndSendFluxCorrections, md_in);
      auto recv_flxcor =
          tl.AddTask(start_flxcor, parthenon::ReceiveFluxCorrections, md_in);
      flux_res = tl.AddTask(recv_flxcor, parthenon::SetFluxCorrections, md_in);
    }
    return tl.AddTask(flux_res, FluxMultiplyMatrix, md_in, md_out);
  }

  // Approximate diagonal of A (exact on uniform grid). Used by Jacobi smoother.
  parthenon::TaskStatus SetDiagonal(std::shared_ptr<parthenon::MeshData<Real>> &md_mat,
                                    std::shared_ptr<parthenon::MeshData<Real>> &md_diag) {
    using namespace parthenon;
    const int ndim = md_mat->GetMeshPointer()->ndim;
    IndexRange ib = md_mat->GetBoundsI(IndexDomain::interior, te);
    IndexRange jb = md_mat->GetBoundsJ(IndexDomain::interior, te);
    IndexRange kb = md_mat->GetBoundsK(IndexDomain::interior, te);

    auto desc_diag = parthenon::MakePackDescriptor<var_t>(md_diag.get());
    auto pack = desc_diag.GetPack(md_diag.get());
    parthenon::par_for(
        "SG::StoreDiagonal", 0, pack.GetNBlocks() - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &coords = pack.GetCoordinates(b);
          Real diag_elem = 0.0;
          // X1
          {
            const Real dxc = coords.template Dxc<parthenon::X1DIR>(k, j, i);
            diag_elem -= 2.0 / (dxc * dxc);
          }
          if (ndim > 1) {
            const Real dxc = coords.template Dxc<parthenon::X2DIR>(k, j, i);
            diag_elem -= 2.0 / (dxc * dxc);
          }
          if (ndim > 2) {
            const Real dxc = coords.template Dxc<parthenon::X3DIR>(k, j, i);
            diag_elem -= 2.0 / (dxc * dxc);
          }
          pack(b, te, var_t(), k, j, i) = diag_elem;
        });
    return TaskStatus::complete;
  }

  // flux = -grad(phi) on each face.
  // Sign convention: A returns -grad^2(phi).  Combined with rhs = +4piG*rho,
  // the system  A.phi = rhs  is equivalent to  -grad^2(phi) = 4piG*rho,
  // i.e. phi here satisfies grad^2(phi) = -4piG*rho.
  // This is Artemis's convention. Acceleration = -grad(phi) in ApplyGravitySource
  // is then correct (pulls toward overdensities). DO NOT change this without
  // understanding what you're doing.
  static parthenon::TaskStatus
  CalculateFluxes(std::shared_ptr<parthenon::MeshData<Real>> &md) {
    using namespace parthenon;
    const int ndim = md->GetMeshPointer()->ndim;
    IndexRange ib = md->GetBoundsI(IndexDomain::interior, te);
    IndexRange jb = md->GetBoundsJ(IndexDomain::interior, te);
    IndexRange kb = md->GetBoundsK(IndexDomain::interior, te);

    auto desc = parthenon::MakePackDescriptor<var_t>(md.get(), {}, {PDOpt::WithFluxes});
    auto pack = desc.GetPack(md.get());
    parthenon::par_for(
        "SG::CalculateFluxes", 0, pack.GetNBlocks() - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &coords = pack.GetCoordinates(b);
          // X1: flux at left face of cell (k,j,i)
          pack.flux(b, X1DIR, var_t(), k, j, i) =
              (pack(b, te, var_t(), k, j, i - 1) - pack(b, te, var_t(), k, j, i)) /
              coords.template Dxc<X1DIR>(k, j, i);
          if (i == ib.e) {
            pack.flux(b, X1DIR, var_t(), k, j, i + 1) =
                (pack(b, te, var_t(), k, j, i) - pack(b, te, var_t(), k, j, i + 1)) /
                coords.template Dxc<X1DIR>(k, j, i + 1);
          }
          if (ndim > 1) {
            pack.flux(b, X2DIR, var_t(), k, j, i) =
                (pack(b, te, var_t(), k, j - 1, i) - pack(b, te, var_t(), k, j, i)) /
                coords.template Dxc<X2DIR>(k, j, i);
            if (j == jb.e) {
              pack.flux(b, X2DIR, var_t(), k, j + 1, i) =
                  (pack(b, te, var_t(), k, j, i) - pack(b, te, var_t(), k, j + 1, i)) /
                  coords.template Dxc<X2DIR>(k, j + 1, i);
            }
          }
          if (ndim > 2) {
            pack.flux(b, X3DIR, var_t(), k, j, i) =
                (pack(b, te, var_t(), k - 1, j, i) - pack(b, te, var_t(), k, j, i)) /
                coords.template Dxc<X3DIR>(k, j, i);
            if (k == kb.e) {
              pack.flux(b, X3DIR, var_t(), k + 1, j, i) =
                  (pack(b, te, var_t(), k, j, i) - pack(b, te, var_t(), k + 1, j, i)) /
                  coords.template Dxc<X3DIR>(k + 1, j, i);
            }
          }
        });
    return TaskStatus::complete;
  }

  // A.x = divergence of fluxes (from CalculateFluxes, possibly corrected at c/f)
  static parthenon::TaskStatus
  FluxMultiplyMatrix(std::shared_ptr<parthenon::MeshData<Real>> &md,
                     std::shared_ptr<parthenon::MeshData<Real>> &md_out) {
    using namespace parthenon;
    const int ndim = md->GetMeshPointer()->ndim;
    IndexRange ib = md->GetBoundsI(IndexDomain::interior, te);
    IndexRange jb = md->GetBoundsJ(IndexDomain::interior, te);
    IndexRange kb = md->GetBoundsK(IndexDomain::interior, te);

    auto desc = parthenon::MakePackDescriptor<var_t>(md.get(), {}, {PDOpt::WithFluxes});
    auto desc_out = parthenon::MakePackDescriptor<var_t>(md_out.get());
    auto pack = desc.GetPack(md.get());
    auto pack_out = desc_out.GetPack(md_out.get());
    parthenon::par_for(
        "SG::FluxMultiplyMatrix", 0, pack.GetNBlocks() - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          const auto &coords = pack.GetCoordinates(b);
          const Real Vol = coords.template Volume<parthenon::TopologicalElement::CC>(k, j, i);
          Real out = 0.0;
          {
            const Real Al = coords.template Volume<parthenon::TopologicalElement::F1>(k, j, i);
            const Real Ar = coords.template Volume<parthenon::TopologicalElement::F1>(k, j, i + 1);
            out += (pack.flux(b, X1DIR, var_t(), k, j, i) * Al -
                    pack.flux(b, X1DIR, var_t(), k, j, i + 1) * Ar) / Vol;
          }
          if (ndim > 1) {
            const Real Al = coords.template Volume<parthenon::TopologicalElement::F2>(k, j, i);
            const Real Ar = coords.template Volume<parthenon::TopologicalElement::F2>(k, j + 1, i);
            out += (pack.flux(b, X2DIR, var_t(), k, j, i) * Al -
                    pack.flux(b, X2DIR, var_t(), k, j + 1, i) * Ar) / Vol;
          }
          if (ndim > 2) {
            const Real Al = coords.template Volume<parthenon::TopologicalElement::F3>(k, j, i);
            const Real Ar = coords.template Volume<parthenon::TopologicalElement::F3>(k + 1, j, i);
            out += (pack.flux(b, X3DIR, var_t(), k, j, i) * Al -
                    pack.flux(b, X3DIR, var_t(), k + 1, j, i) * Ar) / Vol;
          }
          pack_out(b, te, var_t(), k, j, i) = out;
        });
    return TaskStatus::complete;
  }
};

} // namespace SelfGravity

#endif // SELF_GRAVITY_POISSON_EQUATION_HPP_
