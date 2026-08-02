//========================================================================================
// AthenaPK - conservation-budget history diagnostics (VALIDATION_PLAN.md WP-6).
//
// WHY THIS EXISTS. Three anomalies are visible in the production history and NONE of them
// is currently explicable from the columns that exist:
//
//   1. tot-E rises 3.2377e3 -> 6.9593e3 (+115%) over the nj4 leg. That is the EXPECTED
//      behaviour of a collapse if tot-E excludes the gravitational potential energy --
//      which it does. So energy conservation is at present not merely unverified, it is
//      UNVERIFIABLE: the conserved combination is not being written out. Adding
//      W = 1/2 int rho Phi dV makes E_tot + W either close or visibly not close.
//   2. Mass RISES by 1.25e-3 relative over the run. Outflow boundaries cannot add mass.
//      The two candidates -- inflow through a nominally outflow face, and a density floor
//      injecting mass -- are distinguished by instrumenting the boundary flux and the
//      floor SEPARATELY, which is what the cons-Mout and cons-nfloor* columns do.
//   3. Linear momentum drifts ~10% (1-mom -23.89 -> -26.37), which needs the same
//      boundary accounting, but with the FULL stress rather than just the advective part.
//
// These are BUDGETS, not constancy checks. Production uses outflow fluid boundaries, so
// mass, momentum, energy and angular momentum are all non-conserved BY CONSTRUCTION. The
// only meaningful statement is that the change in a volume integral equals the measured
// flux through the faces plus the measured source terms. A "conservation test" that simply
// asserts constancy would fail here for entirely correct reasons, which is precisely why
// the anomalies above have gone unexplained.
//
// Gated on hydro/cons_diag (default OFF). Pure read-only reductions => bit-identical off.
//
// UNITS: Heaviside-Lorentz code units (CLAUDE.md). Magnetic pressure is B^2/2 and the
// Maxwell stress is -B_i B_j, both with no 4pi. four_pi_G = 1, so the Poisson equation the
// solver actually solves is laplacian(Phi) = rho.
//
// Columns (all UserHistoryOperation::sum unless noted):
//   cons-W        = 1/2 int rho Phi dV     gravitational potential energy
//   cons-Mout     = surf rho (v.nhat) dA   mass flux OUT through the domain faces
//   cons-Poutx/y/z= surf (T.nhat)_i dA     TOTAL momentum flux out, including pressure and
//                                          Maxwell stress, not just the advective part
//   cons-nfloor   = number of cells sitting at the density floor        (see CAVEAT)
//   cons-Mfloor   = int rho dV over those cells                         (see CAVEAT)
//   cons-npfloor  = number of cells sitting at the pressure floor       (see CAVEAT)
//
// with the momentum flux tensor
//
//     T_ij = rho v_i v_j + (P + B^2/2) delta_ij - B_i B_j
//
// so the budgets to close are
//
//     d/dt(mass)  = -cons-Mout
//     d/dt(i-mom) = -cons-Pouti
//     d/dt(tot-E + cons-W) = -(energy flux, NOT instrumented -- see below)
//
// SIGN CONVENTION: all *out columns are OUTWARD-POSITIVE, so a positive value means the
// quantity is LEAVING the box, hence the minus signs above. A NEGATIVE cons-Mout is
// therefore the direct, unambiguous detection of inflow through an outflow face -- anomaly
// 2's first candidate -- with no inference required.
//
// CAVEAT ON THE FLOOR COLUMNS -- read before quoting them. These count cells whose density
// (or pressure) currently sits AT the floor value, detected as rho <= dfloor*(1+1e-9).
// That is a PROXY, not the exact mass injected. The exact quantity is the sum over the run
// of (dfloor - rho_before_flooring) * dV at each application, which cannot be recovered
// after the fact and requires an accumulator inside the EOS conserved-to-primitive kernel.
// What these columns CAN do is decide anomaly 2: if cons-nfloor is 0 for the whole run, the
// floor is not the mass source and the boundary is, and vice versa. Use them for that
// dichotomy and do not quote cons-Mfloor as "the mass the floor added".
//
// NOT INSTRUMENTED, stated so the budget is not mistaken for complete:
//   * the total-energy flux through the faces. With a TABULATED EOS the internal energy is
//     not P/(gamma-1), so the flux needs an EOS call per boundary cell rather than an
//     algebraic expression in the primitives; it is deliberately deferred rather than
//     computed with a wrong gamma-law shortcut that would look plausible and be wrong.
//   * the exact floor accumulator described above.
//
// cons-Mout-solver — THE FACTOR-2, RESOLVED (2026-07-31).
//
// RESULT: there is NO conservation bug. The apparent factor of 2 is a STAGE-SAMPLING
// mismatch between the two columns, and every budget closes once each integrator's own
// stage weights are used. An earlier note in this file (and an escalation to the Parthenon
// developers) claimed "vl2 applies the domain-boundary flux at half weight". That claim was
// WRONG and is retracted; the measurements that refute it are below.
//
// Parthenon's low-storage update is u0 <- gam0*u0 + gam1*u1 + beta*dt*L(u0), and the net
// weight applied to the flux divergence over a full step is EXACTLY 1 for rk1, vl2, rk2, rk3
// and rk34 alike (verified by hand from low_storage_integrator.cpp). There is no factor of 2
// available in the coefficients. What differs is WHICH STAGE each column samples:
//
//   cons-Mout        is built from `prim`, i.e. the END-OF-STEP state u^{n+1}.
//   cons-Mout-solver reads the flux array, i.e. the flux computed from the state entering
//                    the LAST stage, u^{(n-1)}.
//
// Measured budget closure (64^3, cfl=0.0375, t->0.3, self-gravity on, runs/wp6_{rk1,faces,rk2}):
//
//   rk1  U^{n+1}=U^n+dt*L^1              predicts -dM/dt = F_cc              measured 1.009
//   vl2  U^{n+1}=U^n+dt*L^2 (gam0=0)     predicts -dM/dt = F_solver          measured 0.991
//   rk2  U^{n+1}=U^n+dt/2(L^1+L^2)       predicts 0.5F_cc+0.5F_solver        measured 0.747
//                                                 = 0.75*F_solver
//
// PRODUCTION USES vl2, whose final stage has gam0 = 0 -- the stage-1 increment is discarded
// and U^{n+1} depends on the last stage's flux ALONE. So for production the solver flux IS
// the applied flux and the mass budget closes to 0.9% (the residual is trapezoid sampling
// between history rows, not a defect). USE cons-Mout-solver TO CLOSE THE BUDGET; cons-Mout is
// a physical cross-check of the surface integral, not the applied flux.
//
// Empirical law across five integrators (runs/wp6_{rk1,faces,rk2,rk3,rk34}), exact to 5
// digits, constant in time from the first history row, and independent of the flow:
//
//     cons-Mout / cons-Mout-solver = beta[nstages-1] / beta[nstages-2]
//        rk1 0.978(*)   vl2 2.0000   rk2 0.5000   rk3 2.6667(=8/3)   rk34 3.0001
//        (*) rk1 is single-stage, so its ratio is the one genuine physics comparison
//            (F at u^n vs F at u^{n+1}) and correctly drifts, 0.993 -> 0.974.
//
// Two competing fits -- (prod_{s<last} beta[s])/beta[last] and beta[last-1]/beta[last] --
// agreed on rk1/vl2/rk2/rk3 and were separated by rk34 (predictions 12.0 vs 3.0, measured
// 3.0001), which killed the product form. Candidate mechanisms ruled out by direct
// inspection: no in-place scaling of the flux array anywhere in Parthenon's update path, and
// SelfGravity::ApplyGravitySource only READS cons fluxes (self_gravity.cpp:544-552), never
// writes them. Gravity and the other unsplit sources are applied per stage with beta*dt and
// their net per-step weight is 1 for every integrator above (checked algebraically).
//
// HONEST LIMIT: the exact arithmetic reason the last-stage boundary flux scales as
// beta[nstages-2] is NOT isolated to a line of code and is not claimed here. It does not
// affect the production budget, because vl2's gam0=0 makes the last-stage flux the applied
// flux by construction. For an integrator with gam0 != 0 on its final stage (rk2, rk3, rk34)
// NEITHER column closes the budget on its own and a stage-weighted accumulator inside the
// update would be required; that is deliberately not built, since production is vl2.
//
// SAME CAVEAT APPLIES TO WP-4. am-FT*/am-FL* and cons-Pout* are likewise built from
// end-of-step cell-centred primitives, so under a multistage integrator they are stage-
// inconsistent with the applied fluxes in exactly the same way. This did not affect the WP-4
// acceptance result (its surface flux was <= 3.4e-11, i.e. negligible against an 8e-2 dL/dt),
// but it MUST be accounted for before any production-condition surface budget is called closed.
//
// Indexing: flux(X1DIR, IDN, k, j, i) is the flux through the LOWER face of cell i, so the
// upper face of cell i is index i+1. This is the same convention `SelfGravity::ApplyGravitySource`
// relies on when it pairs flux(...,i) with the left phi-difference and flux(...,i+1) with the
// right one. The array holds a flux DENSITY (per unit area), so the face contribution is
// flux * dA with dA = dV/dx, matching the conservative update dU = -(dt/dx)(F_{i+1} - F_i).
//
// Sign: outward-positive like every other *out column, so the -x face contributes -F_x and the
// +x face contributes +F_x.
//
// TIMING CAVEAT: this reads whatever is in the flux array when history runs, i.e. the LAST
// stage's fluxes. For VL2 that is the stage-2 flux, which is exactly the one the state update
// U^{n+1} = U^n + dt*L(U*) actually used -- so it is the right object to compare against. But
// if a later operator-split source (radiation, chemistry) clears or overwrites the fluxes, this
// column reads something stale. A zero or wildly-off value is that failure, not a physics
// result: sanity-check it against cons-Mout having the same sign and order of magnitude.
//
// CAVEAT ON cons-W. W = 1/2 int rho Phi dV is only meaningful when Phi is the physical
// potential of an isolated system. With self_gravity/*_bc = zero (Dirichlet Phi = 0 on the
// box face) the potential is clamped at the boundary and W carries a box-size-dependent
// offset; with *_bc = multipole (WS-5a) it does not. Compare W across runs only at matched
// box size and matched gravity BC.
//========================================================================================
#ifndef DIAGNOSTICS_CONS_DIAG_HPP_
#define DIAGNOSTICS_CONS_DIAG_HPP_

#include <parthenon/package.hpp>

using namespace parthenon::package::prelude;

namespace Diagnostics {

//! Which conservation-budget reduction to perform.
enum class ConsDiag {
  Wgrav,  // 1/2 int rho Phi dV                  (needs the "grav.phi" field)
  Mout,   // surf rho (v.nhat) dA
  PoutX,  // surf (T.nhat)_x dA
  PoutY,
  PoutZ,
  nfloor,  // count of cells at the density floor
  Mfloor,  // int rho dV over cells at the density floor
  npfloor, // count of cells at the pressure floor
  // The SOLVER's own mass flux through the domain faces, read from the cons flux array
  // rather than reconstructed from the cell-centred primitives. See the note below.
  MoutSolver,
  // Inner-face-only / outer-face-only splits of BOTH the cell-centred and the solver
  // versions. These exist to distinguish two explanations of the factor 2 that the
  // symmetric test problem cannot otherwise tell apart: (A) the applied face flux really is
  // half of rho_cell*v_cell on every face, versus (B) the solver's flux array is populated
  // on the inner faces only, so the outer three of the six contribute zero. In a problem
  // with six near-identical faces both produce a ratio of exactly 2.
  MoutSolverInner,
  MoutSolverOuter,
  MoutInner,
  MoutOuter
};

//! Conservation-budget reduction. Volume terms run over interior cells; the *out terms run
//! over the cells adjacent to a physical domain face. See the file header for the sign
//! convention (outward-positive) and for what is deliberately not instrumented.
Real ConsDiagReduce(MeshData<Real> *md, ConsDiag which);

} // namespace Diagnostics

#endif // DIAGNOSTICS_CONS_DIAG_HPP_
