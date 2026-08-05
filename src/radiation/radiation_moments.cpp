//========================================================================================
// AthenaPK - M1 moment radiation transport package.
// Increment 2b: explicit, operator-split M1 transport (free-streaming + diffusion).
// Donor-cell HLL Riemann fluxes for the two-moment system
//     d_t Er   + d_j  Fr_j            = 0   (matter coupling added in increment 3)
//     d_t Fr_i + d_j (chat^2 P_ij)    = 0 ,   P_ij = D_ij(f) Er  (M1 closure)
// with reduced speed of light chat and causal flux limit |Fr| <= chat Er.
//========================================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <interface/update.hpp>
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

#include "../eos/eos_table.hpp" // tabulated protostellar EOS (RT owns gas energy)
#include "../main.hpp" // IDN, IM1..IM3, IEN, NHYDRO, IB1..IB3 (cons component indices)
#include "radiation.hpp"
#include "radiation_closure.hpp"
#include "radiation_groups.hpp"
#include "radiation_opacity.hpp"
#include "../dust/dust.hpp" // WS-4 DustFactor / DustModel (opacity consumer)

namespace Radiation {

using namespace parthenon::package::prelude;
using parthenon::X1DIR;
using parthenon::X2DIR;
using parthenon::X3DIR;

//----------------------------------------------------------------------------------------
//! Physical M1 flux 4-vector G[Er,Frx,Fry,Frz] across a face normal to `dir` (1,2,3),
//! plus the M1 min/max signal speeds for the given cell state. (M1 closure only.)
KOKKOS_INLINE_FUNCTION
void M1FaceFlux(const int dir, const Real E, const Real Fx, const Real Fy, const Real Fz,
                const Real chat, Real G[4], Real &smin, Real &smax) {
  const Real iceE = 1.0 / (chat * E + RadFuzz());
  const std::array<Real, 3> fred{Fx * iceE, Fy * iceE, Fz * iceE};
  const Real fraw =
      std::sqrt(RadSqr(fred[0]) + RadSqr(fred[1]) + RadSqr(fred[2])) + RadFuzz();
  const Real fmag = std::min(1.0, fraw);
  const std::array<Real, 6> D = EddingtonTensor<Closure::M1>(fred); // {xx,yy,zz,yz,xz,xy}
  const Real c2 = chat * chat;

  Real nd;
  if (dir == 1) {
    nd = fred[0] / fraw;
    G[0] = Fx;
    G[1] = c2 * E * D[0];
    G[2] = c2 * E * D[5];
    G[3] = c2 * E * D[4];
  } else if (dir == 2) {
    nd = fred[1] / fraw;
    G[0] = Fy;
    G[1] = c2 * E * D[5];
    G[2] = c2 * E * D[1];
    G[3] = c2 * E * D[3];
  } else {
    nd = fred[2] / fraw;
    G[0] = Fz;
    G[1] = c2 * E * D[4];
    G[2] = c2 * E * D[3];
    G[3] = c2 * E * D[2];
  }
  const Kokkos::pair<Real, Real> lam = WaveSpeed<Closure::M1>(nd, fmag);
  smin = chat * lam.first;
  smax = chat * lam.second;
}

//----------------------------------------------------------------------------------------
//! Donor-cell HLL fluxes for (Er,Fr1,Fr2,Fr3) into the WithFluxes flux arrays.
//! Minmod slope limiter for the WS-3b PLM reconstruction (2nd order, TVD).
KOKKOS_FORCEINLINE_FUNCTION
Real RadMinmod(const Real a, const Real b) {
  return (a * b <= 0.0) ? 0.0 : ((std::abs(a) < std::abs(b)) ? a : b);
}

TaskStatus CalculateRadFluxes(MeshData<Real> *md) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
  const int ndim = pmb->pmy_mesh->ndim;

  auto pkg = pmb->pmy_mesh->packages.Get("radiation");
  PARTHENON_REQUIRE(pkg->Param<std::string>("closure") == "M1",
                    "Increment 2b transport supports closure=M1 only.");
  const Real chat = pkg->Param<Real>("chat");
  const bool plm = pkg->Param<bool>("rad_recon_plm"); // WS-3b: dc (false) or plm (true)
  const int n_group = pkg->Param<int>("n_group");

  // Multigroup: each group's (Er_g,Fr_g) is an INDEPENDENT M1 system (same chat, same
  // frequency-independent closure); they couple only through the matter temperature in
  // MatterCoupling. Loop the identical HLL sweep per group. n_group=1 runs it once over the
  // gray field names -> bit-identical to the pre-multigroup transport.
  for (int g = 0; g < n_group; ++g) {
  const std::vector<std::string> names = GroupFieldNames(g);
  parthenon::PackIndexMap imap;
  auto pack = md->PackVariablesAndFluxes(names, imap);
  const int iE = imap[names[0]].first;
  const int iX = imap[names[1]].first;
  const int iY = imap[names[2]].first;
  const int iZ = imap[names[3]].first;

  auto hll = KOKKOS_LAMBDA(const int dir, const int b, const int k, const int j,
                           const int i, const int dk, const int dj, const int di) {
    Real GL[4], GR[4], smnL, smxL, smnR, smxR;
    // Reconstruct the left/right face states of (Er, Fr1, Fr2, Fr3). The face at (i)
    // sits between the "minus" cell (i-di) and the "plus" cell (i). Donor-cell uses the
    // cell centers directly (1st order); PLM adds a minmod-limited half-slope (2nd order,
    // +-2 stencil). plm=false is bit-identical to the pre-WS-3b donor-cell flux.
    const int comp[4] = {iE, iX, iY, iZ};
    Real UL[4], UR[4];
    for (int c = 0; c < 4; ++c) {
      const int ic = comp[c];
      const Real cm1 = pack(b, ic, k - dk, j - dj, i - di); // minus cell (i-di)
      const Real c0 = pack(b, ic, k, j, i);                 // plus cell  (i)
      if (plm) {
        const Real cm2 = pack(b, ic, k - 2 * dk, j - 2 * dj, i - 2 * di); // (i-2di)
        const Real cp1 = pack(b, ic, k + dk, j + dj, i + di);            // (i+di)
        UL[c] = cm1 + 0.5 * RadMinmod(cm1 - cm2, c0 - cm1); // right face of minus cell
        UR[c] = c0 - 0.5 * RadMinmod(c0 - cm1, cp1 - c0);   // left face of plus cell
      } else {
        UL[c] = cm1;
        UR[c] = c0;
      }
    }
    M1FaceFlux(dir, UL[0], UL[1], UL[2], UL[3], chat, GL, smnL, smxL);
    M1FaceFlux(dir, UR[0], UR[1], UR[2], UR[3], chat, GR, smnR, smxR);
    const Real SL = std::min(std::min(smnL, smnR), 0.0);
    const Real SR = std::max(std::max(smxL, smxR), 0.0);
    const Real idS = 1.0 / (SR - SL + RadFuzz());
    auto &v = pack(b);
    const int dirv = (dir == 1) ? X1DIR : ((dir == 2) ? X2DIR : X3DIR);
    v.flux(dirv, iE, k, j, i) = (SR * GL[0] - SL * GR[0] + SL * SR * (UR[0] - UL[0])) * idS;
    v.flux(dirv, iX, k, j, i) = (SR * GL[1] - SL * GR[1] + SL * SR * (UR[1] - UL[1])) * idS;
    v.flux(dirv, iY, k, j, i) = (SR * GL[2] - SL * GR[2] + SL * SR * (UR[2] - UL[2])) * idS;
    v.flux(dirv, iZ, k, j, i) = (SR * GL[3] - SL * GR[3] + SL * SR * (UR[3] - UL[3])) * idS;
  };

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Rad::FluxX1", parthenon::DevExecSpace(), 0,
      pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e + 1,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        hll(1, b, k, j, i, 0, 0, 1);
      });
  if (ndim >= 2) {
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "Rad::FluxX2", parthenon::DevExecSpace(), 0,
        pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e + 1, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          hll(2, b, k, j, i, 0, 1, 0);
        });
  }
  if (ndim >= 3) {
    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, "Rad::FluxX3", parthenon::DevExecSpace(), 0,
        pack.GetDim(5) - 1, kb.s, kb.e + 1, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
          hll(3, b, k, j, i, 1, 0, 0);
        });
  }
  } // group loop
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Explicit FV update U += dt*(-div F) followed by the causal clamp.
TaskStatus ApplyRadUpdate(MeshData<Real> *md, const Real dt) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
  const int ndim = pmb->pmy_mesh->ndim;

  auto pkg = pmb->pmy_mesh->packages.Get("radiation");
  const Real chat = pkg->Param<Real>("chat");
  const Real efloor = pkg->Param<Real>("efloor");
  const int n_group = pkg->Param<int>("n_group");

  for (int g = 0; g < n_group; ++g) {
  const std::vector<std::string> names = GroupFieldNames(g);
  parthenon::PackIndexMap imap;
  auto pack = md->PackVariablesAndFluxes(names, imap);
  const int iE = imap[names[0]].first;
  const int iX = imap[names[1]].first;
  const int iY = imap[names[2]].first;
  const int iZ = imap[names[3]].first;

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Rad::Update", parthenon::DevExecSpace(), 0,
      pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = pack.GetCoords(b);
        auto &v = pack(b);
        Real E = pack(b, iE, k, j, i) +
                 dt * parthenon::Update::FluxDivHelper(iE, k, j, i, ndim, coords, v);
        Real Fx = pack(b, iX, k, j, i) +
                  dt * parthenon::Update::FluxDivHelper(iX, k, j, i, ndim, coords, v);
        Real Fy = pack(b, iY, k, j, i) +
                  dt * parthenon::Update::FluxDivHelper(iY, k, j, i, ndim, coords, v);
        Real Fz = pack(b, iZ, k, j, i) +
                  dt * parthenon::Update::FluxDivHelper(iZ, k, j, i, ndim, coords, v);
        E = std::max(E, efloor);
        const Real fmag = std::sqrt(Fx * Fx + Fy * Fy + Fz * Fz) / (chat * E + RadFuzz());
        if (fmag > 1.0) {
          const Real s = 1.0 / fmag;
          Fx *= s;
          Fy *= s;
          Fz *= s;
        }
        pack(b, iE, k, j, i) = E;
        pack(b, iX, k, j, i) = Fx;
        pack(b, iY, k, j, i) = Fy;
        pack(b, iZ, k, j, i) = Fz;
      });
  } // group loop
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Gas energy-balance residual R(T) for the multigroup implicit coupling (groups eliminated
//! analytically): R(T) = [e(T)-e0] - (c/chat) sum_g a_g (E_g^0 - B_g(T))/(1+a_g), with
//! B_g = arad T^4 PlanckFraction(g,T), a_g = chat dt rho kappa_P,g. Free function (not a
//! nested lambda) so it is unambiguously device-safe under nvcc. Root -> the coupled (T,E_g).
KOKKOS_INLINE_FUNCTION
Real MGResidual(const Real Tt, const Real rho, const Real kdust, const Real arad,
                const Real chat, const Real c, const Real dt, const Real gm1,
                const Real T_unit, const bool use_h2, const EOSTable::EosTable &eos_tab,
                const OpacityParams &op, const RadGroups &groups,
                const GroupOpacityTable &optab, const int n_group, const Real *Eg0,
                const Real e0_ref) {
  const Real kP = PlanckOpacity(op, rho, Tt) * kdust; // gray Planck opacity at Tt
  const Real e = use_h2 ? eos_tab.EintFromRhoTk(rho, Tt * T_unit) : rho * Tt / gm1;
  Real S = 0.0;
  const Real T4 = Tt * Tt * Tt * Tt;
  for (int g = 0; g < n_group; ++g) {
    const Real Bg = arad * T4 * groups.PlanckFraction(g, Tt * T_unit);
    const Real ag = chat * dt * rho * kP * optab.PlanckMult(g, Tt * T_unit);
    S += ag * (Eg0[g] - Bg) / (1.0 + ag);
  }
  return (e - e0_ref) - c / chat * S;
}

//----------------------------------------------------------------------------------------
//! Implicit MULTIGROUP matter coupling (n_group>1). One matter temperature T couples to all
//! groups: group g emits B_g = arad*T^4*PlanckFraction(g,T) and absorbs at a_g = chat*dt*rho*
//! kappa_P,g. The (n_group) group-energy equations are linear in E_g given T
//!     E_g^new = (E_g^0 + a_g B_g(T)) / (1 + a_g),
//! so E_g^new - B_g = (E_g^0 - B_g)/(1+a_g), and the gas energy balance closes on a SINGLE
//! nonlinear equation in T (Newton):
//!     R(T) = [e(T) - e0] - (c/chat) * sum_g a_g (E_g^0 - B_g(T))/(1+a_g) = 0 .
//! With a nu-independent opacity (beta=0 => band mult=1) and sum_g B_g = arad*T^4, this reduces EXACTLY
//! to the gray solve (same T, same sum_g E_g). Groups couple ONLY through T. Flux attenuation
//! + radiation-force momentum exchange are applied per group (Rosseland mean + scattering),
//! and the kinetic radiation-energy exchange is removed pro-rata across groups.
TaskStatus MatterCouplingMultigroup(MeshData<Real> *md, const Real dt) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  auto pkg = pmb->pmy_mesh->packages.Get("radiation");
  const Real chat = pkg->Param<Real>("chat");
  const Real c = pkg->Param<Real>("c");
  const Real arad = pkg->Param<Real>("arad");
  const OpacityParams op = pkg->Param<OpacityParams>("opacity");
  const Real gam = pkg->Param<Real>("gamma");
  const Real gm1 = gam - 1.0;
  const Real tfloor = pkg->Param<Real>("tfloor");
  const Real efloor = pkg->Param<Real>("efloor");
  const int inner_max = pkg->Param<int>("inner_iteration_max");
  const Real inner_tol = pkg->Param<Real>("inner_iteration_tol");
  const Real coupling_tmax = pkg->Param<Real>("coupling_tmax"); // rtsafe upper bracket (code T)
  const Real Bfloor = arad * tfloor * tfloor * tfloor * tfloor;
  const bool use_h2 = pkg->Param<bool>("use_h2diss");
  const auto eos_tab = pkg->Param<EOSTable::EosTable>("eos_tab");
  const Real T_unit = pkg->Param<Real>("T_unit");
  const RadGroups groups = pkg->Param<RadGroups>("groups");
  const GroupOpacityTable optab = pkg->Param<GroupOpacityTable>("optable");
  const int n_group = groups.n_group;

  // Pack gas cons + ALL groups' moments.
  std::vector<std::string> names{"cons"};
  for (const auto &nm : AllRadFieldNames(n_group)) names.push_back(nm);
  parthenon::PackIndexMap imap;
  auto pack = md->PackVariables(names, imap);
  const int ic = imap["cons"].first;
  const int nhydro = pmb->pmy_mesh->packages.Get("Hydro")->Param<int>("nhydro");
  const bool mhd = (nhydro > NHYDRO);
  // Per-group component indices (Kokkos::Array is captured by value into the device kernel).
  Kokkos::Array<int, MAX_GROUP> iEg, iXg, iYg, iZg;
  for (int g = 0; g < n_group; ++g) {
    const auto gn = GroupFieldNames(g);
    iEg[g] = imap[gn[0]].first;
    iXg[g] = imap[gn[1]].first;
    iYg[g] = imap[gn[2]].first;
    iZg[g] = imap[gn[3]].first;
  }

  // WS-4 dust consumer (same as the gray path).
  auto &pkgs = pmb->pmy_mesh->packages;
  const bool dust_on = pkgs.AllPackages().count("dust") > 0 &&
                       pkgs.Get("dust")->Param<bool>("evolve");
  int dust_sidx = 0;
  Real dust_fref = 0.01, dust_aref = 1.0e-5;
  if (dust_on) {
    auto dpkg = pkgs.Get("dust");
    dust_sidx = dpkg->Param<int>("scalar_index");
    const Dust::DustModel dm = dpkg->Param<Dust::DustModel>("model");
    dust_fref = dm.f_dg_ref;
    dust_aref = dm.a_ref;
  }
  const int idust = ic + nhydro + dust_sidx;

  int nfail = 0;
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "Rad::MatterCouplingMG",
      parthenon::DevExecSpace(), 0, pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, int &lnfail) {
        const Real rho = pack(b, ic + IDN, k, j, i);
        const Real m1 = pack(b, ic + IM1, k, j, i);
        const Real m2 = pack(b, ic + IM2, k, j, i);
        const Real m3 = pack(b, ic + IM3, k, j, i);
        Real kdust = 1.0;
        if (dust_on) {
          const Real irho = 1.0 / rho;
          kdust = Dust::DustFactor(pack(b, idust, k, j, i) * irho,
                                   pack(b, idust + 1, k, j, i) * irho, dust_fref, dust_aref);
        }
        const Real ke = 0.5 * (m1 * m1 + m2 * m2 + m3 * m3) / rho;
        Real me = 0.0;
        if (mhd) {
          const Real b1 = pack(b, ic + IB1, k, j, i);
          const Real b2 = pack(b, ic + IB2, k, j, i);
          const Real b3 = pack(b, ic + IB3, k, j, i);
          me = 0.5 * (b1 * b1 + b2 * b2 + b3 * b3);
        }
        Real e0 = pack(b, ic + IEN, k, j, i) - ke - me;
        Real T = use_h2 ? ((e0 > 0.0) ? eos_tab.TemperatureK(rho, e0) / T_unit : tfloor)
                        : gm1 * e0 / rho;
        T = std::max(T, tfloor);
        e0 = use_h2 ? eos_tab.EintFromRhoTk(rho, T * T_unit) : rho * T / gm1;
        const Real e0_ref = e0;

        // Read the per-group incident moments.
        Real Eg0[MAX_GROUP], Fxg0[MAX_GROUP], Fyg0[MAX_GROUP], Fzg0[MAX_GROUP];
        Real Esum0 = 0.0;
        for (int g = 0; g < n_group; ++g) {
          Eg0[g] = pack(b, iEg[g], k, j, i);
          Fxg0[g] = pack(b, iXg[g], k, j, i);
          Fyg0[g] = pack(b, iYg[g], k, j, i);
          Fzg0[g] = pack(b, iZg[g], k, j, i);
          Esum0 += Eg0[g];
        }
        const Real etot = e0_ref + c / chat * Esum0 + RadFuzz();

        // Solve MGResidual(T) = 0 by SAFEGUARDED Newton+bisection (rtsafe, NR). The residual is
        // MONOTONICALLY INCREASING in T (both e(T) and the emission a_R T^4 rise with T), so a
        // bracketed method is guaranteed to converge and can NEVER write garbage. The plain
        // numerical-derivative Newton previously here overshot near the tabulated-EOS e(T) kinks
        // (H2 dissociation / H ionization) at collapse onset -> non-convergence -> NaN blowup.
        auto Rof = [&](const Real Tt) {
          return MGResidual(Tt, rho, kdust, arad, chat, c, dt, gm1, T_unit, use_h2, eos_tab, op,
                            groups, optab, n_group, Eg0, e0_ref);
        };
        // Bracket [Tlo, Thi]. Thi capped inside the tabulated EOS/opacity range (Tmax_table).
        const Real Tguess = T; // warm start: pre-coupling gas temperature (a good initial guess)
        Real Tlo = tfloor;
        Real Thi = coupling_tmax; // code T; opacity/EOS-table Tmax (stay in range)
        if (Thi <= Tlo) Thi = 1.0e5; // fallback (ideal EOS: no table cap)
        // EARLY EXIT: if the gas is already in radiative equilibrium (residual already ~0 at the
        // current T), skip the whole solve. In relaxed regions most cells are near equilibrium, so
        // this avoids the ~30-bisection cold-start over the huge [tfloor, ~1e6 K] bracket.
        const Real fguess = Rof(Tguess);
        bool conv = false;
        if (std::abs(fguess) / etot <= inner_tol) {
          T = Tguess; conv = true;             // already in radiative equilibrium
        } else if (fguess > 0.0 && Rof(Tlo) >= 0.0) {
          T = Tlo; conv = true;                // root at/below floor: gas cools to tfloor
        } else if (fguess < 0.0 && Rof(Thi) <= 0.0) {
          T = Thi; conv = true;                // root above table cap: clamp (unphysical; safe)
        }
        if (conv) {
          // handled above
        } else {
          // rtsafe with a WARM START at Tguess (Newton from a good guess converges in a few steps;
          // the [Tlo,Thi] bracket remains the safety net). Tighten the bracket with fguess's sign.
          if (fguess < 0.0) Tlo = Tguess; else Thi = Tguess;
          T = Tguess; // == clamp(Tguess,[Tlo,Thi]) after the tightening above
          Real dTold = Thi - Tlo, dTstep = dTold;
          Real f = fguess; // reuse: Rof(Tguess) already computed for the early-exit check
          for (int it = 0; it < inner_max; ++it) {
            const Real dT = 1.0e-4 * T + 1.0e-12;
            const Real df = (Rof(T + dT) - f) / dT; // numerical derivative (monotone => df>0)
            // bisect if Newton would leave the bracket or converge too slowly
            if (((T - Thi) * df - f) * ((T - Tlo) * df - f) > 0.0 ||
                std::abs(2.0 * f) > std::abs(dTold * df)) {
              dTold = dTstep;
              dTstep = 0.5 * (Thi - Tlo);
              T = Tlo + dTstep;
            } else {
              dTold = dTstep;
              dTstep = f / (df + RadFuzz());
              T = T - dTstep;
            }
            if (std::abs(dTstep) < inner_tol * T || std::abs(f) / etot <= inner_tol) {
              conv = true;
              break;
            }
            f = Rof(T);
            if (f < 0.0) { Tlo = T; } else { Thi = T; } // maintain the bracket
          }
        }
        T = std::max(T, tfloor);
        if (!conv) lnfail += 1;

        // Converged T -> per-group implicit energies + gas thermal change.
        const Real kP = PlanckOpacity(op, rho, T) * kdust;
        const Real T4 = T * T * T * T;
        Real Egn[MAX_GROUP];
        Real Esum_new = 0.0;
        for (int g = 0; g < n_group; ++g) {
          const Real Bg = arad * T4 * groups.PlanckFraction(g, T * T_unit);
          const Real ag = chat * dt * rho * kP * optab.PlanckMult(g, T * T_unit);
          Egn[g] = std::max((Eg0[g] + ag * Bg) / (1.0 + ag), efloor);
          Esum_new += Egn[g];
        }
        const Real dEg = (use_h2 ? eos_tab.EintFromRhoTk(rho, T * T_unit) : rho * T / gm1) -
                         e0_ref;

        // Flux attenuation + radiation-force momentum exchange, summed over groups.
        Real dFx[MAX_GROUP], dFy[MAX_GROUP], dFz[MAX_GROUP];
        Real sdFx = 0.0, sdFy = 0.0, sdFz = 0.0;
        const Real kR = RosselandOpacity(op, rho, T) * kdust;
        // audit N5 (2026-08-05): scattering here is DUST scattering, so it scales with the
        // dust factor exactly as the two absorption means do. Previously it did not, so a
        // sublimating cell lost its absorption but kept full-strength scattering. Latent, and
        // the reason is the DUST switch, not the scattering coefficient: exactly two decks in
        // the tree set radiation/kappa_s_code non-zero (runs/validation_rt/diff_thick.in =
        // 100.0 and runs/validation_rt/rad_shadow.in = 0.1) and NEITHER enables <physics>
        // dust, so kdust == 1 in both and this line is bit-identical everywhere it has run.
        const Real ks = ScatteringOpacity(op, rho, T) * kdust;
        for (int g = 0; g < n_group; ++g) {
          const Real ag = chat * dt * rho * (kR * optab.RossMult(g, T * T_unit) + ks);
          const Real fac = -ag / (1.0 + ag);
          dFx[g] = fac * Fxg0[g];
          dFy[g] = fac * Fyg0[g];
          dFz[g] = fac * Fzg0[g];
          sdFx += dFx[g];
          sdFy += dFy[g];
          sdFz += dFz[g];
        }
        const Real icc = 1.0 / (c * chat * rho);
        const Real dvx = -icc * sdFx, dvy = -icc * sdFy, dvz = -icc * sdFz;
        const Real vx = m1 / rho, vy = m2 / rho, vz = m3 / rho;
        const Real vnx = vx + dvx, vny = vy + dvy, vnz = vz + dvz;
        const Real dEk = 0.5 * rho * ((vnx * vnx - vx * vx) + (vny * vny - vy * vy) +
                                      (vnz * vnz - vz * vz));

        // Write back. Radiation energy: the emission/absorption part is in Egn; the kinetic
        // radiation-work chat/c*dEk is removed pro-rata across groups (Sum removal exact up
        // to per-group flooring). Then apply flux attenuation and gas momentum/energy.
        const Real kin = chat / c * dEk;
        const Real inv_Esum = 1.0 / (Esum_new + RadFuzz());
        for (int g = 0; g < n_group; ++g) {
          pack(b, iEg[g], k, j, i) = std::max(Egn[g] - kin * Egn[g] * inv_Esum, efloor);
          pack(b, iXg[g], k, j, i) = Fxg0[g] + dFx[g];
          pack(b, iYg[g], k, j, i) = Fyg0[g] + dFy[g];
          pack(b, iZg[g], k, j, i) = Fzg0[g] + dFz[g];
        }
        pack(b, ic + IM1, k, j, i) = m1 + dvx * rho;
        pack(b, ic + IM2, k, j, i) = m2 + dvy * rho;
        pack(b, ic + IM3, k, j, i) = m3 + dvz * rho;
        pack(b, ic + IEN, k, j, i) += dEg + dEk;
      },
      Kokkos::Sum<int>(nfail));

  if (nfail > 0) {
    static bool warned_mg = false;
    if (!warned_mg) {
      warned_mg = true;
      std::cout << "### RADIATION WARNING: multigroup matter coupling did not converge in "
                << nfail << " cell(s) this step (inner_iteration_max=" << inner_max
                << ", tol=" << inner_tol << "); last iterate written. (warn-once)"
                << std::endl;
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Implicit gray matter coupling (RT owns the gas thermal energy).
//!
//! Per cell, solves the operator-split source system  y = U^(0) + dt*S(y)  for the
//! radiation energy E and the Planck function B = arad*T^4 via the linearized Newton
//! step of Artemis' "Simple" coupling (radiation/moments/matter_coupling.hpp), then
//! applies the flux--momentum (radiation force) exchange explicitly with the
//! momentum-conserving sign of Artemis' "Full" coupling (see note in the kernel):
//!     S_E   =  chat * rho*kappa_a * (E - B)            (absorption - emission)
//!     S_e   = -c    * rho*kappa_a * (E - B)            (gas internal energy)
//!     S_Fr  = -chat * rho*(kappa_a+kappa_s) * Fr       (flux attenuation)
//! Gas EOS is ideal: e_int = rho*T/(gamma-1), Cv = rho/(gamma-1) (code units, p=rho*T).
//! No barotropic floor is applied here -- this REPLACES collapse_be's e_th overwrite.
TaskStatus MatterCoupling(MeshData<Real> *md, const Real dt) {
  auto pmb = md->GetBlockData(0)->GetBlockPointer();
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  auto pkg = pmb->pmy_mesh->packages.Get("radiation");
  // Multigroup: dispatch to the group-coupled solve. n_group=1 falls through to the gray
  // path below, byte-for-byte unchanged (production/prod_v9 is gray).
  if (pkg->Param<int>("n_group") > 1) return MatterCouplingMultigroup(md, dt);
  const Real chat = pkg->Param<Real>("chat");
  const Real c = pkg->Param<Real>("c");
  const Real arad = pkg->Param<Real>("arad");
  const OpacityParams op = pkg->Param<OpacityParams>("opacity"); // kappa(rho,T), code units
  const Real gam = pkg->Param<Real>("gamma");
  const Real gm1 = gam - 1.0;
  const Real tfloor = pkg->Param<Real>("tfloor");
  const Real efloor = pkg->Param<Real>("efloor");
  const int inner_max = pkg->Param<int>("inner_iteration_max");
  const Real inner_tol = pkg->Param<Real>("inner_iteration_tol");
  const Real Bfloor = arad * tfloor * tfloor * tfloor * tfloor;
  // Tabulated-EOS coupling: map gas internal-energy DENSITY e <-> code temperature T via
  // the shared table instead of the ideal e=rho*T/gm1. T (code) = T_phys[K]/T_unit.
  const bool use_h2 = pkg->Param<bool>("use_h2diss");
  const auto eos_tab = pkg->Param<EOSTable::EosTable>("eos_tab");
  const Real T_unit = pkg->Param<Real>("T_unit");

  // Pack gas cons + radiation together (cons holds IDN..IEN [+ IB1..IB3 if MHD]).
  static const std::vector<std::string> names{"cons", "rad.Er", "rad.Fr1", "rad.Fr2",
                                               "rad.Fr3"};
  parthenon::PackIndexMap imap;
  auto pack = md->PackVariables(names, imap);
  const int ic = imap["cons"].first; // cons component base: IDN=ic+0, IEN=ic+IEN, ...
  // F2 fix (audit 2026-07-18): detect MHD from the Hydro package's nhydro (glmmhd =>
  // nhydro > NHYDRO), NOT from the packed cons width -- passive scalars widen "cons",
  // so the old (ncons >= NHYDRO+3) test misfired for euler + >=3 scalars and read
  // scalar components as B-field in the magnetic-energy subtraction.
  const int nhydro = pmb->pmy_mesh->packages.Get("Hydro")->Param<int>("nhydro");
  const bool mhd = (nhydro > NHYDRO);
  const int iE = imap["rad.Er"].first;
  const int iX = imap["rad.Fr1"].first;
  const int iY = imap["rad.Fr2"].first;
  const int iZ = imap["rad.Fr3"].first;

  // WS-4 dust consumer: if the dust package is evolving grains, scale kappa per cell by
  // DustFactor = (f_dg/f_ref)(a_ref/a_c) (geometric-area limit). The dust scalars live in
  // the packed "cons" field at [nhydro+scalar_index, +1]. Inert (factor=1) if dust off/frozen.
  auto &pkgs = pmb->pmy_mesh->packages;
  const bool dust_on = pkgs.AllPackages().count("dust") > 0 &&
                       pkgs.Get("dust")->Param<bool>("evolve");
  // (nhydro already fetched above for the F2 mhd test)
  int dust_sidx = 0;
  Real dust_fref = 0.01, dust_aref = 1.0e-5;
  if (dust_on) {
    auto dpkg = pkgs.Get("dust");
    dust_sidx = dpkg->Param<int>("scalar_index");
    const Dust::DustModel dm = dpkg->Param<Dust::DustModel>("model");
    dust_fref = dm.f_dg_ref;
    dust_aref = dm.a_ref;
  }
  const int idust = ic + nhydro + dust_sidx; // f_dg at idust, a_c at idust+1

  int nfail = 0; // audit #6: count cells where the implicit matter coupling did NOT converge
  parthenon::par_reduce(
      parthenon::loop_pattern_mdrange_tag, "Rad::MatterCoupling", parthenon::DevExecSpace(),
      0, pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i, int &lnfail) {
        const Real rho = pack(b, ic + IDN, k, j, i);
        const Real m1 = pack(b, ic + IM1, k, j, i);
        const Real m2 = pack(b, ic + IM2, k, j, i);
        const Real m3 = pack(b, ic + IM3, k, j, i);
        // WS-4 dust opacity factor (1.0 if dust off/frozen at reference => bit-identical).
        Real kdust = 1.0;
        if (dust_on) {
          const Real irho = 1.0 / rho;
          kdust = Dust::DustFactor(pack(b, idust, k, j, i) * irho,
                                   pack(b, idust + 1, k, j, i) * irho, dust_fref, dust_aref);
        }
        const Real ke = 0.5 * (m1 * m1 + m2 * m2 + m3 * m3) / rho;
        Real me = 0.0;
        if (mhd) {
          const Real b1 = pack(b, ic + IB1, k, j, i);
          const Real b2 = pack(b, ic + IB2, k, j, i);
          const Real b3 = pack(b, ic + IB3, k, j, i);
          me = 0.5 * (b1 * b1 + b2 * b2 + b3 * b3); // Heaviside-Lorentz: P_mag = B^2/2
        }

        // Gas internal energy and temperature. Ideal: p=rho*T, e=rho*T/gm1. H2-diss: the
        // Saha EOS maps e <-> T_phys[K]; T (code) = T_phys/T_unit.
        Real e0 = pack(b, ic + IEN, k, j, i) - ke - me;
        // Guard the table path against e0 <= 0: TemperatureK interpolates in log10(e0/rho)
        // and would return NaN, and std::max(NaN, tfloor) propagates the NaN (the ideal
        // branch goes negative and is caught by the tfloor max). Fall through to tfloor.
        Real T = use_h2 ? ((e0 > 0.0) ? eos_tab.TemperatureK(rho, e0) / T_unit : tfloor)
                        : gm1 * e0 / rho;
        T = std::max(T, tfloor);
        // re-sync e0 with T (avoids spurious 0-opacity drift)
        e0 = use_h2 ? eos_tab.EintFromRhoTk(rho, T * T_unit) : rho * T / gm1;
        Real B = arad * T * T * T * T;
        const Real e0_ref = e0;

        const Real Er0 = pack(b, iE, k, j, i);
        const Real Fx0 = pack(b, iX, k, j, i);
        const Real Fy0 = pack(b, iY, k, j, i);
        const Real Fz0 = pack(b, iZ, k, j, i);

        Real E = Er0;
        const Real etot = e0_ref + c / chat * E + RadFuzz();

        // Newton iteration on (E, B): linearized implicit absorption/emission.
        bool conv = false;
        for (int it = 0; it < inner_max; ++it) {
          T = std::pow(B / arad, 0.25);
          const Real e = use_h2 ? eos_tab.EintFromRhoTk(rho, T * T_unit)
                                : rho * T / gm1; // InternalEnergyFromDensityTemperature
          Real Cv;                               // d e_int / dT (code T)
          if (use_h2) {
            const Real dTc = 1.0e-4 * T + 1.0e-12;
            Cv = (eos_tab.EintFromRhoTk(rho, (T + dTc) * T_unit) - e) / dTc;
          } else {
            Cv = rho / gm1;
          }
          // absorption coefficient (1/length); opacity may depend on the current T.
          // Emission/absorption uses the PLANCK mean (kappa_P).
          const Real a = chat * dt * rho * PlanckOpacity(op, rho, T) * kdust;
          const Real fleck = FleckFactor(arad, T, Cv);
          const Real Ri = a * (E - B);
          const Real Fi = (e - e0_ref) - c / chat * Ri;
          const Real Fr = (E - Er0) + Ri;
          const Real idet = 1.0 / (1.0 + a + c / chat * fleck * a);
          const Real dE = ((1.0 + c / chat * fleck * a) * (-Fr) + a * (-fleck * Fi)) * idet;
          const Real dB = ((c / chat * fleck * a) * (-Fr) + (1.0 + a) * (-fleck * Fi)) * idet;
          E = std::max(E + dE, efloor);
          B = std::max(B + dB, Bfloor);
          const Real err =
              std::max(std::abs(Fi) / etot, c / chat * std::abs(Fr) / etot);
          if (err <= inner_tol) {
            conv = true;
            break;
          }
        }
        if (!conv) lnfail += 1; // audit #6: non-silent non-convergence

        // Gas internal-energy change from the converged temperature.
        T = std::pow(B / arad, 0.25);
        const Real dEg = (use_h2 ? eos_tab.EintFromRhoTk(rho, T * T_unit)
                                 : rho * T / gm1) - e0_ref;

        // Flux--momentum (radiation force) exchange, explicit in the attenuation factor.
        // Sign convention follows Artemis' momentum-conserving "Full" coupling
        // (matter_coupling.hpp: icc = +1/(c chat rho), dv = -icc*dF): the gas gains the
        // momentum the radiation flux loses, dp_gas = -dF/(c*chat). NOTE: Artemis'
        // "Simple" coupling (the original source of this port) has icc = -1/(c chat rho),
        // which kicks the gas ANTI-parallel to the attenuated flux -- wrong direction.
        // Flux attenuation uses the ROSSELAND mean (kappa_R) + scattering.
        const Real a = chat * dt * rho *
                       (RosselandOpacity(op, rho, T) + ScatteringOpacity(op, rho, T)) * kdust;
        const Real dFx = -a / (1.0 + a) * Fx0;
        const Real dFy = -a / (1.0 + a) * Fy0;
        const Real dFz = -a / (1.0 + a) * Fz0;
        const Real icc = 1.0 / (c * chat * rho);
        const Real dvx = -icc * dFx, dvy = -icc * dFy, dvz = -icc * dFz;
        const Real vx = m1 / rho, vy = m2 / rho, vz = m3 / rho;
        const Real vnx = vx + dvx, vny = vy + dvy, vnz = vz + dvz;
        const Real dEk = 0.5 * rho *
                         ((vnx * vnx - vx * vx) + (vny * vny - vy * vy) +
                          (vnz * vnz - vz * vz));

        // Conservative update of both gas and radiation (Artemis "Full" pairing):
        // whatever total energy the gas gains (thermal dEg + kinetic dEk), the radiation
        // field loses, scaled by chat/c in the reduced-speed-of-light convention.
        // Floor the conservative write: dEg+dEk can (in pathological floored cells) exceed
        // the radiation budget, and the next CalculateRadFluxes would see a negative Er
        // before ApplyRadUpdate's floor. Matches the flooring already done in the Newton
        // loop (exact conservation is intentionally sacrificed in floored cells only).
        pack(b, iE, k, j, i) = std::max(Er0 - chat / c * (dEg + dEk), efloor);
        pack(b, iX, k, j, i) = Fx0 + dFx;
        pack(b, iY, k, j, i) = Fy0 + dFy;
        pack(b, iZ, k, j, i) = Fz0 + dFz;
        pack(b, ic + IM1, k, j, i) = m1 + dvx * rho;
        pack(b, ic + IM2, k, j, i) = m2 + dvy * rho;
        pack(b, ic + IM3, k, j, i) = m3 + dvz * rho;
        pack(b, ic + IEN, k, j, i) += dEg + dEk;
      },
      Kokkos::Sum<int>(nfail));

  // Audit #6: make non-convergence NON-SILENT. The implicit solve previously wrote the
  // last iterate regardless of residual; now count and report it (warn-once) so it cannot
  // masquerade as physical heating/cooling. (Full per-cell rejection/subcycling is a
  // larger follow-up noted in the audit remediation plan.)
  if (nfail > 0) {
    static bool warned = false;
    if (!warned) {
      warned = true;
      std::cout << "### RADIATION WARNING: implicit matter coupling did not converge in "
                << nfail << " cell(s) this step (inner_iteration_max=" << inner_max
                << ", tol=" << inner_tol
                << "); last iterate written. Raise radiation/inner_iteration_max or reduce"
                   " dt if this persists. (warn-once)"
                << std::endl;
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! Radiation CFL timestep over the whole mesh: cfl * dx_min / chat.
Real RadDtMesh(Mesh *pmesh) {
  auto pkg = pmesh->packages.Get("radiation");
  const Real chat = pkg->Param<Real>("chat");
  const Real cfl = pkg->Param<Real>("cfl");
  Real dxmin = std::numeric_limits<Real>::max();
  for (auto const &pmb : pmesh->block_list) {
    const auto &bs = pmb->block_size;
    for (int d = 0; d < pmesh->ndim; d++) {
      dxmin = std::min(dxmin, (bs.xmax_[d] - bs.xmin_[d]) / bs.nx_[d]);
    }
  }
#ifdef MPI_PARALLEL
  PARTHENON_MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, &dxmin, 1, MPI_PARTHENON_REAL, MPI_MIN,
                                    MPI_COMM_WORLD));
#endif
  return cfl * dxmin / chat;
}

//----------------------------------------------------------------------------------------
//! Operator-split M1 transport: sub-cycle (flux -> update -> ghost exchange) at the
//! radiation CFL over the hydro step dt. Mirrors SelfGravity::SolvePoisson structure.
void AddRadiationTasks(TaskCollection &tc, Mesh *pmesh, const Real dt) {
  using namespace parthenon;
  TaskID none(0);
  if (getenv("RAD_DISABLE_TRANSPORT") != nullptr) return; // diagnostic: no-op transport

  const Real dt_rad = RadDtMesh(pmesh);
  const int nsub = std::max(1, static_cast<int>(std::ceil(dt / (dt_rad + RadFuzz()))));
  const Real dts = dt / static_cast<Real>(nsub);
  if (getenv("RAD_PRINT_NSUB") != nullptr && parthenon::Globals::my_rank == 0)
    printf("[RAD_NSUB] nsub=%d dt=%.3e dt_rad=%.3e\n", nsub, dt, dt_rad), fflush(stdout);

  auto pkg = pmesh->packages.Get("radiation");
  const bool do_coupling = pkg->Param<bool>("matter_coupling");

  // Multigroup: the rad-only sub-container (flux correction + ghost exchange) must carry
  // ALL groups' moments. For n_group=1 this is exactly the four gray names.
  const std::vector<std::string> rad_names = AllRadFieldNames(pkg->Param<int>("n_group"));
  auto partitions = pmesh->GetDefaultBlockPartitions();
  const int num_partitions = partitions.size();
  TaskRegion &region = tc.AddRegion(num_partitions);
  for (int i = 0; i < num_partitions; ++i) {
    TaskList &tl = region[i];
    auto &md = pmesh->mesh_data.Add("base", partitions[i]);
    auto &md_rad = pmesh->mesh_data.Add("rad_sub", md, rad_names);
    TaskID prev = none;
    for (int s = 0; s < nsub; ++s) {
      // Audit fix A1: reflux the radiation moments across coarse-fine interfaces. rad.Er /
      // rad.Fr are registered WithFluxes, but transport previously applied face fluxes
      // with NO flux correction, so coarse and fine faces used different net radiation
      // flux at a refinement boundary -> non-conservative Er/Fr (spurious heating/cooling
      // that MatterCoupling then transfers into the gas). Mirror the hydro reflux
      // sequence, restricted to the rad-only container so only rad fluxes are corrected.
      auto start_flxcor = pmesh->multilevel
                              ? tl.AddTask(prev, parthenon::StartReceiveFluxCorrections, md_rad)
                              : prev;
      auto flx = tl.AddTask(prev, CalculateRadFluxes, md.get());
      TaskID ready = flx;
      if (pmesh->multilevel) {
        auto send = tl.AddTask(flx, parthenon::LoadAndSendFluxCorrections, md_rad);
        auto recv = tl.AddTask(start_flxcor, parthenon::ReceiveFluxCorrections, md_rad);
        ready = tl.AddTask(recv | flx, parthenon::SetFluxCorrections, md_rad);
      }
      auto upd = tl.AddTask(ready, ApplyRadUpdate, md.get(), dts);
      auto bcs = AddBoundaryExchangeTasks(upd, tl, md_rad, pmesh->multilevel);
      prev = bcs;
    }
    // Matter coupling (absorption/emission + radiation force) once over the full hydro dt,
    // after transport. Local per-cell for the update itself, BUT it mutates all four
    // interior radiation moments.
    if (do_coupling) {
      auto mc = tl.AddTask(prev, MatterCoupling, md.get(), dt);
      // Audit fix A2: refresh radiation ghosts AFTER coupling so the next step's first
      // CalculateRadFluxes reconstructs interface states from coupled (not stale
      // pre-coupling) neighbor moments -- otherwise the first flux of every step carries a
      // one-step, decomposition-dependent impulse at block/rank/AMR boundaries.
      AddBoundaryExchangeTasks(mc, tl, md_rad, pmesh->multilevel);
    }
  }
}

} // namespace Radiation
