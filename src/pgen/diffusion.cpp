
// AthenaPK - a performance portable block structured AMR MHD code
// Copyright (c) 2021-2024, Athena Parthenon Collaboration. All rights reserved.
// Licensed under the 3-Clause License (the "LICENSE")

// Parthenon headers
#include "mesh/mesh.hpp"
#include <cmath>
#include <complex>
#include <iostream>
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

// AthenaPK headers
#include "../hydro/ct/ct.hpp"
#include "../main.hpp"
#include "utils/error_checking.hpp"

namespace diffusion {
using namespace parthenon::driver::prelude;

void ProblemGenerator(MeshBlock *pmb, ParameterInput *pin) {
  auto hydro_pkg = pmb->packages.Get("Hydro");
  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  auto &mbd = pmb->meshblock_data.Get();
  auto &u = mbd->Get("cons").data;

  const auto gamma = pin->GetReal("hydro", "gamma");
  const bool mhd_enabled = hydro_pkg->Param<Fluid>("fluid") == Fluid::glmmhd;

  const auto Bx = pin->GetOrAddReal("problem/diffusion", "Bx", 0.0);
  const auto By = pin->GetOrAddReal("problem/diffusion", "By", 0.0);

  const auto iprob = pin->GetInteger("problem/diffusion", "iprob");
  PARTHENON_REQUIRE_THROWS(mhd_enabled || !(iprob == 0 || iprob == 1 || iprob == 2 ||
                                            iprob == 10 || iprob == 20 || iprob == 40),
                           "Selected iprob for diffusion pgen requires MHD enabled.")
  Real t0 = 0.5;
  Real diff_coeff = 0.0;
  Real amp = 1e-6;
  // Get common parameters for Gaussian profile
  if ((iprob == 10) || (iprob == 30) || (iprob == 40)) {
    t0 = pin->GetOrAddReal("problem/diffusion", "t0", t0);
    amp = pin->GetOrAddReal("problem/diffusion", "amp", amp);
  }
  // Heat diffusion of 1D Gaussian
  if (iprob == 10) {
    diff_coeff = pin->GetReal("diffusion", "thermal_diff_coeff_code");
    // Viscous diffusion of 1D Gaussian
  } else if (iprob == 30) {
    diff_coeff = pin->GetReal("diffusion", "mom_diff_coeff_code");
    // Ohmic diffusion of 1D Gaussian
  } else if (iprob == 40) {
    diff_coeff = pin->GetReal("diffusion", "ohm_diff_coeff_code");
  }

  // Ambipolar diffusion of a 1D sinusoidal transverse field (eigenmode decay test).
  // With a uniform guide field B0 along x and By = amp*sin(k x), the current is purely
  // perpendicular to B, so E_A = eta_A J and By decays as exp(-eta_A k^2 t) with
  // eta_A = Q_A B0^2. The analytic decay rate is therefore Q_A * B0^2 * k^2.
  Real kpar = 0.0;
  if (iprob == 50) {
    PARTHENON_REQUIRE_THROWS(mhd_enabled, "iprob=50 (ambipolar decay) requires MHD.");
    amp = pin->GetOrAddReal("problem/diffusion", "amp", amp);
    const auto nwave = pin->GetOrAddReal("problem/diffusion", "nwave", 1.0);
    const auto x1min = pin->GetReal("parthenon/mesh", "x1min");
    const auto x1max = pin->GetReal("parthenon/mesh", "x1max");
    kpar = 2.0 * M_PI * nwave / (x1max - x1min);
  }

  // Hall whistler / ion-cyclotron circularly polarized eigenmode test. A traveling wave
  // b+/b- (set with its matching velocity eigenvector) propagates at omega/k given by the
  // Hall dispersion relation omega = -h*alpha + sqrt(alpha^2 + k^2 v_A^2), with
  // alpha = (Q_H/rho0) B0 k^2 / 2 and helicity h = +/-1. The two helicities split (whistler
  // vs ion-cyclotron) -- the signature of the Hall term.
  Real hall_h = 1.0, hall_cvel = 0.0;
  if (iprob == 60) {
    PARTHENON_REQUIRE_THROWS(mhd_enabled, "iprob=60 (Hall whistler) requires MHD.");
    amp = pin->GetOrAddReal("problem/diffusion", "amp", amp);
    const auto nwave = pin->GetOrAddReal("problem/diffusion", "nwave", 1.0);
    const auto x1min = pin->GetReal("parthenon/mesh", "x1min");
    const auto x1max = pin->GetReal("parthenon/mesh", "x1max");
    kpar = 2.0 * M_PI * nwave / (x1max - x1min);
    hall_h = pin->GetOrAddReal("problem/diffusion", "helicity", 1.0);
    const Real QH = pin->GetReal("diffusion", "hall_coeff_code");
    const Real rho0 = 1.0;
    const Real vA2 = Bx * Bx / rho0;
    const Real alpha = (QH / rho0) * Bx * kpar * kpar / 2.0;
    const Real omega = -hall_h * alpha + std::sqrt(alpha * alpha + kpar * kpar * vA2);
    hall_cvel = (Bx / rho0) * (kpar / omega); // v = -hall_cvel * b (eigenvector)
  }

  auto &coords = pmb->coords;

  pmb->par_for(
      "ProblemGenerator: Diffusion", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        u(IDN, k, j, i) = 1.0;

        u(IM1, k, j, i) = 0.0;
        u(IM2, k, j, i) = 0.0;
        u(IM3, k, j, i) = 0.0;

        if (mhd_enabled) {
          u(IB1, k, j, i) = 0.0;
          u(IB2, k, j, i) = 0.0;
          u(IB3, k, j, i) = 0.0;
        }

        Real eint = -1.0;
        // step function x1
        if (iprob == 0) {
          u(IB1, k, j, i) = Bx;
          u(IB2, k, j, i) = By;
          eint = coords.Xc<1>(i) <= 0.0 ? 10.0 : 12.0;
          // step function x2
        } else if (iprob == 1) {
          u(IB2, k, j, i) = Bx;
          u(IB3, k, j, i) = By;
          eint = coords.Xc<2>(j) <= 0.0 ? 10.0 : 12.0;
          // step function x3
        } else if (iprob == 2) {
          u(IB3, k, j, i) = Bx;
          u(IB1, k, j, i) = By;
          eint = coords.Xc<3>(k) <= 0.0 ? 10.0 : 12.0;
          // Gaussian
        } else if (iprob == 10) {
          u(IB1, k, j, i) = Bx;
          u(IB2, k, j, i) = By;
          // Adjust for anisotropic thermal conduction.
          // If there's no conduction for the setup (because the field is perp.)
          // treat as 1 (also in analysis) to prevent division by 0.
          // Note, this is very constructed and needs to be updated/adjusted for isotropic
          // conduction, other directions, and Bfield configs with |B| != 1
          Real eff_diff_coeff = Bx == 0.0 ? diff_coeff : diff_coeff * Bx * Bx;
          eint = 1 + amp / std::sqrt(4. * M_PI * eff_diff_coeff * t0) *
                         std::exp(-(std::pow(coords.Xc<1>(i), 2.)) /
                                  (4. * eff_diff_coeff * t0));
          // Ring diffusion in x1-x2 plane
        } else if (iprob == 20) {
          const auto x = coords.Xc<1>(i);
          const auto y = coords.Xc<2>(j);
          Real r = std::sqrt(SQR(x) + SQR(y));
          Real phi = std::atan2(y, x);

          u(IB1, k, j, i) = y / r;
          u(IB2, k, j, i) = -x / r;
          eint = std::abs(r - 0.6) < 0.1 && std::abs(phi) < M_PI / 12.0 ? 12.0 : 10.0;
          // Ring diffusion in x2-x3 plane
        } else if (iprob == 21) {
          const auto x = coords.Xc<2>(j);
          const auto y = coords.Xc<3>(k);
          Real r = std::sqrt(SQR(x) + SQR(y));
          Real phi = std::atan2(y, x);

          u(IB2, k, j, i) = y / r;
          u(IB3, k, j, i) = -x / r;
          eint = std::abs(r - 0.6) < 0.1 && std::abs(phi) < M_PI / 12.0 ? 12.0 : 10.0;
          // Ring diffusion in x3-x1 plane
        } else if (iprob == 22) {
          const auto x = coords.Xc<3>(k);
          const auto y = coords.Xc<1>(i);
          Real r = std::sqrt(SQR(x) + SQR(y));
          Real phi = std::atan2(y, x);

          u(IB3, k, j, i) = y / r;
          u(IB1, k, j, i) = -x / r;
          eint = std::abs(r - 0.6) < 0.1 && std::abs(phi) < M_PI / 12.0 ? 12.0 : 10.0;
          // Viscous diffusion of 1D Gaussian
        } else if (iprob == 30) {
          u(IM2, k, j, i) =
              u(IDN, k, j, i) * amp /
              std::pow(std::sqrt(4. * M_PI * diff_coeff * t0), 1.0) *
              std::exp(-(std::pow(coords.Xc<1>(i), 2.)) / (4. * diff_coeff * t0));
          eint = 1.0 / (gamma * (gamma - 1.0)); // c_s = 1 everywhere
          // Ohmic diffusion of 1D Gaussian
        } else if (iprob == 40) {
          u(IB2, k, j, i) =
              amp / std::pow(std::sqrt(4. * M_PI * diff_coeff * t0), 1.0) *
              std::exp(-(std::pow(coords.Xc<1>(i), 2.)) / (4. * diff_coeff * t0));
          eint = 1.0 / (gamma * (gamma - 1.0)); // c_s = 1 everywhere
          // Ambipolar diffusion of a 1D sinusoidal transverse field
        } else if (iprob == 50) {
          u(IB1, k, j, i) = Bx; // uniform guide field B0 along x
          u(IB2, k, j, i) = amp * std::sin(kpar * coords.Xc<1>(i));
          eint = 1.0 / (gamma * (gamma - 1.0)); // c_s = 1 everywhere
          // Hall circularly polarized whistler / ion-cyclotron eigenmode
        } else if (iprob == 60) {
          const Real xx = coords.Xc<1>(i);
          const Real by = amp * std::cos(kpar * xx);
          const Real bz = hall_h * amp * std::sin(kpar * xx);
          u(IB1, k, j, i) = Bx; // guide field B0 along x
          u(IB2, k, j, i) = by;
          u(IB3, k, j, i) = bz;
          // Matching velocity eigenvector v = -(B0/rho0)(k/omega) b
          u(IM2, k, j, i) = -hall_cvel * by * u(IDN, k, j, i);
          u(IM3, k, j, i) = -hall_cvel * bz * u(IDN, k, j, i);
          eint = 1.0 / (gamma * (gamma - 1.0)); // c_s = 1 everywhere
        }

        PARTHENON_REQUIRE(eint > 0.0, "Missing init of eint");
        u(IEN, k, j, i) =
            u(IDN, k, j, i) * eint +
            0.5 * ((SQR(u(IM1, k, j, i)) + SQR(u(IM2, k, j, i)) + SQR(u(IM3, k, j, i))) /
                   u(IDN, k, j, i));

        if (mhd_enabled) {
          u(IEN, k, j, i) +=
              0.5 * (SQR(u(IB1, k, j, i)) + SQR(u(IB2, k, j, i)) + SQR(u(IB3, k, j, i)));
        }
      });

  // ---- Constrained Transport: initialize the face-centered field Bf ----
  // Two diffusion tests are CT-wired, both with a field that varies only in x so the face
  // init is div-B-free by construction and the projection reproduces the cell-centered u:
  //   iprob=40 (Ohmic Gaussian):     B_y = G(x), B_x=B_z=0  -> F2(i,j-1/2)=G(x_i), F1=F3=0.
  //   iprob=50 (ambipolar eigenmode): B_x=B0 (uniform guide), B_y=amp*sin(k x), B_z=0
  //                                   -> F1(i-1/2)=B0, F2(i,j-1/2)=amp*sin(k x_i), F3=0.
  if (hydro_pkg->Param<bool>("use_ct")) {
    PARTHENON_REQUIRE_THROWS(
        (iprob == 40) || (iprob == 50),
        "Constrained Transport (divergence_control=ct) in the diffusion pgen is currently "
        "only wired for iprob=40 (Ohmic Gaussian) and iprob=50 (ambipolar eigenmode).");
    const bool is_ad = (iprob == 50);
    const Real gnorm = is_ad ? 0.0 : amp / std::sqrt(4. * M_PI * diff_coeff * t0);
    const Real inv4dt = is_ad ? 0.0 : 1.0 / (4. * diff_coeff * t0);
    const Real B0 = Bx;   // guide field (iprob=50); 0 for iprob=40 by input
    const Real kx = kpar; // wavenumber (iprob=50); 0 for iprob=40
    const Real amp_l = amp;
    auto desc = parthenon::MakePackDescriptor<Hydro::CT::Bf>(mbd.get());
    auto pack = desc.GetPack(mbd.get());
    const int bidx = 0;
    using TE = parthenon::TopologicalElement;
    // F1 (x-face, B_x): guide field B0 for iprob=50, else 0.
    IndexRange fib = pmb->cellbounds.GetBoundsI(IndexDomain::interior, TE::F1);
    IndexRange fjb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior, TE::F1);
    IndexRange fkb = pmb->cellbounds.GetBoundsK(IndexDomain::interior, TE::F1);
    pmb->par_for(
        "diffusion CT F1", fkb.s, fkb.e, fjb.s, fjb.e, fib.s, fib.e,
        KOKKOS_LAMBDA(const int k, const int j, const int i) {
          pack(bidx, TE::F1, Hydro::CT::Bf(), k, j, i) = is_ad ? B0 : 0.0;
        });
    // F2 (y-face, B_y): Gaussian(x) (iprob=40) or amp*sin(k x) (iprob=50).
    fib = pmb->cellbounds.GetBoundsI(IndexDomain::interior, TE::F2);
    fjb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior, TE::F2);
    fkb = pmb->cellbounds.GetBoundsK(IndexDomain::interior, TE::F2);
    pmb->par_for(
        "diffusion CT F2", fkb.s, fkb.e, fjb.s, fjb.e, fib.s, fib.e,
        KOKKOS_LAMBDA(const int k, const int j, const int i) {
          const auto &c = pack.GetCoordinates(bidx);
          const Real x = c.X<X1DIR, TE::F2>(k, j, i);
          pack(bidx, TE::F2, Hydro::CT::Bf(), k, j, i) =
              is_ad ? amp_l * std::sin(kx * x) : gnorm * std::exp(-x * x * inv4dt);
        });
    // F3 (z-face, B_z) = 0
    fib = pmb->cellbounds.GetBoundsI(IndexDomain::interior, TE::F3);
    fjb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior, TE::F3);
    fkb = pmb->cellbounds.GetBoundsK(IndexDomain::interior, TE::F3);
    pmb->par_for(
        "diffusion CT F3", fkb.s, fkb.e, fjb.s, fjb.e, fib.s, fib.e,
        KOKKOS_LAMBDA(const int k, const int j, const int i) {
          pack(bidx, TE::F3, Hydro::CT::Bf(), k, j, i) = 0.0;
        });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void UserWorkAfterLoop
//! \brief Quantitative check for the ambipolar diffusion eigenmode test (iprob=50).
//!
//! A transverse perturbation By = amp*sin(kx) on a uniform guide field B0 (density rho0=1,
//! velocity 0) evolves as a *damped Alfven wave*: the Lorentz force launches an Alfven
//! wave while ambipolar diffusion (a perpendicular resistivity eta_A = Q_A B0^2) damps it.
//! The linear dispersion relation is
//!     s^2 + eta_A k^2 s + k^2 v_A^2 = 0 ,   v_A^2 = B0^2/rho0 ,
//! with roots s1, s2. For the initial condition By(0)=amp, vy(0)=0 the exact solution is
//!     By(t) = amp/(s1 - s2) * (s1*exp(s1 t) - s2*exp(s2 t)) * sin(kx) .
//! We measure the By amplitude at the final time and compare to this analytic value. The
//! input is chosen to be (over)damped (eta_A k > 2 v_A) so the mode decays monotonically.
void UserWorkAfterLoop(Mesh *mesh, parthenon::ParameterInput *pin,
                       parthenon::SimTime &tm) {
  const auto iprob = pin->GetInteger("problem/diffusion", "iprob");
  if (iprob != 50 && iprob != 60) return;

  //--------------------------------------------------------------------------------------
  // Hall whistler / ion-cyclotron dispersion check (iprob=60).
  if (iprob == 60) {
    const auto B0 = pin->GetOrAddReal("problem/diffusion", "Bx", 0.0);
    const auto amp = pin->GetOrAddReal("problem/diffusion", "amp", 1e-6);
    const auto nwave = pin->GetOrAddReal("problem/diffusion", "nwave", 1.0);
    const auto helicity = pin->GetOrAddReal("problem/diffusion", "helicity", 1.0);
    const auto x1min = pin->GetReal("parthenon/mesh", "x1min");
    const auto x1max = pin->GetReal("parthenon/mesh", "x1max");
    const Real kpar = 2.0 * M_PI * nwave / (x1max - x1min);
    const Real QH = pin->GetReal("diffusion", "hall_coeff_code");
    const Real tfin = tm.time;
    const Real rho0 = 1.0;
    const Real vA2 = B0 * B0 / rho0;
    const Real alpha = (QH / rho0) * B0 * kpar * kpar / 2.0;
    // Combined Hall + Ohmic-floor dispersion (the floor is a perpendicular resistivity):
    //   omega^2 + (2 h alpha + i eta_O k^2) omega - k^2 v_A^2 = 0 .
    // The physical forward mode is the root with the larger real part; its real part is
    // the oscillation frequency, its imaginary part the (Ohmic) damping rate.
    const Real etaO = pin->GetOrAddReal("diffusion", "hall_ohmic_floor_code", 0.0);
    const std::complex<Real> bc(2.0 * helicity * alpha, etaO * kpar * kpar);
    const std::complex<Real> sq = std::sqrt(bc * bc + 4.0 * kpar * kpar * vA2);
    const std::complex<Real> r1 = 0.5 * (-bc + sq);
    const std::complex<Real> r2 = 0.5 * (-bc - sq);
    const std::complex<Real> root = (std::real(r1) >= std::real(r2)) ? r1 : r2;
    const Real omega_analytic = std::real(root);

    // Project By onto cos(kx) and sin(kx); for By = amp cos(kx - omega t) the phase is
    // omega*t = atan2(<By sin kx>, <By cos kx>). Also track <By^2 + Bz^2> (= amp^2 for a
    // circular mode) to verify the Hall term is non-dissipative (amplitude conserved).
    Real Sc = 0.0, Ss = 0.0, sum_bperp2 = 0.0, vol_tot = 0.0;
    for (auto &pmb : mesh->block_list) {
      IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
      IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
      IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
      auto &rc = pmb->meshblock_data.Get();
      auto u = rc->Get("cons").data.GetHostMirrorAndCopy();
      for (int k = kb.s; k <= kb.e; ++k) {
        for (int j = jb.s; j <= jb.e; ++j) {
          for (int i = ib.s; i <= ib.e; ++i) {
            const Real vol = pmb->coords.CellVolume(k, j, i);
            const Real xx = pmb->coords.Xc<1>(i);
            Sc += u(IB2, k, j, i) * std::cos(kpar * xx) * vol;
            Ss += u(IB2, k, j, i) * std::sin(kpar * xx) * vol;
            sum_bperp2 += (SQR(u(IB2, k, j, i)) + SQR(u(IB3, k, j, i))) * vol;
            vol_tot += vol;
          }
        }
      }
    }
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &Sc, 1, MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Ss, 1, MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum_bperp2, 1, MPI_PARTHENON_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &vol_tot, 1, MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    const Real phase = std::atan2(Ss, Sc); // = omega * t in (-pi, pi]
    const Real omega_meas = (tfin > 0.0) ? phase / tfin : 0.0;
    const Real amp_meas = std::sqrt(sum_bperp2 / vol_tot); // = amp for a circular mode
    const Real rel_err_omega =
        (omega_analytic != 0.0) ? std::abs(omega_meas - omega_analytic) / omega_analytic
                                : 0.0;

    if (parthenon::Globals::my_rank == 0) {
      std::cout << "\n=== Hall whistler/ion-cyclotron dispersion test (iprob=60) ===\n"
                << "  B0 = " << B0 << ", rho0 = " << rho0 << ", k = " << kpar
                << ", Q_H = " << QH << ", helicity = " << helicity << ", t = " << tfin
                << "\n"
                << "  branch = " << (helicity > 0 ? "ion-cyclotron (slow)" : "whistler (fast)")
                << "\n"
                << "  analytic omega = " << omega_analytic
                << " (ideal Alfven k*v_A = " << kpar * std::sqrt(vA2) << ")\n"
                << "  measured omega = " << omega_meas << "\n"
                << "  rel. error (omega)       = " << rel_err_omega << "\n"
                << "  amplitude in/out (Hall non-dissipative): " << amp << " / " << amp_meas
                << "  rel. change = " << std::abs(amp_meas - amp) / amp << std::endl;
    }
    return;
  }

  //--------------------------------------------------------------------------------------
  // Ambipolar diffusion damped-Alfven amplitude check (iprob=50).
  const auto B0 = pin->GetOrAddReal("problem/diffusion", "Bx", 0.0);
  const auto amp = pin->GetOrAddReal("problem/diffusion", "amp", 1e-6);
  const auto nwave = pin->GetOrAddReal("problem/diffusion", "nwave", 1.0);
  const auto x1min = pin->GetReal("parthenon/mesh", "x1min");
  const auto x1max = pin->GetReal("parthenon/mesh", "x1max");
  const Real kpar = 2.0 * M_PI * nwave / (x1max - x1min);
  const Real QA = pin->GetReal("diffusion", "ambipolar_coeff_code");
  const Real tfin = tm.time;

  const Real rho0 = 1.0; // set in the problem generator for iprob=50
  const Real eta_A = QA * B0 * B0;
  const Real vA2 = B0 * B0 / rho0;
  const Real bb = eta_A * kpar * kpar; // = eta_A k^2
  const Real cc = kpar * kpar * vA2;
  const Real disc = bb * bb - 4.0 * cc;

  // Analytic amplitude of the damped Alfven mode at t = tfin.
  Real A_pred;
  if (disc >= 0.0) { // (over)damped: two real roots
    const Real sq = std::sqrt(disc);
    const Real s1 = 0.5 * (-bb - sq); // fast
    const Real s2 = 0.5 * (-bb + sq); // slow
    A_pred = amp / (s1 - s2) * (s1 * std::exp(s1 * tfin) - s2 * std::exp(s2 * tfin));
  } else { // underdamped: complex roots, real envelope form
    const Real omega = 0.5 * std::sqrt(-disc);
    A_pred = amp * std::exp(-0.5 * bb * tfin) *
             (std::cos(omega * tfin) + 0.5 * bb / omega * std::sin(omega * tfin));
  }

  // Volume-weighted mean of By^2 over the domain. For By = A sin(kx), <By^2> = A^2/2.
  Real sum_by2 = 0.0, vol_tot = 0.0;
  for (auto &pmb : mesh->block_list) {
    IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
    IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
    IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
    auto &rc = pmb->meshblock_data.Get();
    auto u = rc->Get("cons").data.GetHostMirrorAndCopy();
    for (int k = kb.s; k <= kb.e; ++k) {
      for (int j = jb.s; j <= jb.e; ++j) {
        for (int i = ib.s; i <= ib.e; ++i) {
          const Real vol = pmb->coords.CellVolume(k, j, i);
          sum_by2 += SQR(u(IB2, k, j, i)) * vol;
          vol_tot += vol;
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &sum_by2, 1, MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &vol_tot, 1, MPI_PARTHENON_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

  const Real A_meas = std::sqrt(2.0 * sum_by2 / vol_tot);
  const Real rel_err_amp = std::abs(A_meas - std::abs(A_pred)) / std::abs(A_pred);

  if (parthenon::Globals::my_rank == 0) {
    std::cout << "\n=== Ambipolar diffusion damped-Alfven test (iprob=50) ===\n"
              << "  B0 = " << B0 << ", rho0 = " << rho0 << ", k = " << kpar
              << ", Q_A = " << QA << ", t = " << tfin << "\n"
              << "  eta_A = Q_A B0^2 = " << eta_A << ",  eta_A k^2 = " << bb
              << ",  v_A = " << std::sqrt(vA2)
              << (disc >= 0.0 ? "  (overdamped)\n" : "  (underdamped)\n")
              << "  analytic amplitude A_pred = " << std::abs(A_pred) << "\n"
              << "  measured amplitude A_meas = " << A_meas << "\n"
              << "  relative error            = " << rel_err_amp << std::endl;
  }
}
} // namespace diffusion
