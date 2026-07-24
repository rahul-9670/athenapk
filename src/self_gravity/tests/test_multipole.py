#!/usr/bin/env python
"""Flagship Phase 7 (gravity BCs): analytic-identity gate for the WS-5a multipole potential.

Read-only. The exterior multipole potential MultipolePhi(mm, x,y,z) in
src/self_gravity/multipole.hpp (monopole + traceless quadrupole) is the boundary lift for the
self-gravity solve. This gate reimplements the SAME formula (mirror of multipole.hpp:47-62 --
keep in sync) and checks the rigorous analytic identities a correct exterior expansion must obey,
turning the WS-5a validation (0.54% vs a full analytic solve) into a committed reproducible test.

Phi = -G (M/r + 0.5 * dQd/r^5),  G = four_pi_G/(4 pi),  dQd = d.Q.d  (Q traceless symmetric).

Gates:
  1. Monopole: Phi = -G M / r exactly, with 1/r scaling.
  2. Quadrupole: P_2 angular structure + 1/r^3 scaling (M=0, traceless Q).
  3. Superposition: Phi(M,Q) == Phi(M,0) + Phi(0,Q) (linear in the moments).
  4. Traceless-Q shell average == 0: the quadrupole term integrates to zero over a sphere
     (Int n_i n_j dOmega = 4pi/3 delta_ij, Tr Q = 0) -- the defining property that keeps the
     monopole the sole far-field mass term.

Run: /beegfs/u/bbg6470/venvs/analysis_env/bin/python src/self_gravity/tests/test_multipole.py
"""
import sys
import numpy as np

FOUR_PI_G = 1.0  # FHC normalization
G = FOUR_PI_G * (0.25 / np.pi)


def multipole_phi(M, Q, com, x, y, z):
    """Mirror of multipole.hpp::MultipolePhi. Q = [Qxx,Qyy,Qzz,Qxy,Qxz,Qyz]."""
    dx, dy, dz = x - com[0], y - com[1], z - com[2]
    r2 = dx * dx + dy * dy + dz * dz
    r = np.sqrt(r2)
    if r <= 0.0:
        return 0.0
    invr = 1.0 / r
    invr5 = invr / (r2 * r2)
    dQd = (Q[0] * dx * dx + Q[1] * dy * dy + Q[2] * dz * dz +
           2.0 * (Q[3] * dx * dy + Q[4] * dx * dz + Q[5] * dy * dz))
    return -G * (M * invr + 0.5 * dQd * invr5)


def main():
    fails = 0
    com = (0.3, -0.2, 0.1)  # arbitrary center of mass (identities are COM-independent)
    M = 6.0
    Q = np.array([-1.0, -1.0, 2.0, 0.4, -0.3, 0.2])  # traceless (Tr = -1-1+2 = 0)

    # Gate 1: monopole.
    print("--- GATE 1: monopole Phi = -G M / r, 1/r scaling ---")
    worst = 0.0
    for (px, py, pz) in [(5, 0, 0), (0, 7, 0), (3, 4, 12), (-8, 2, -3)]:
        r = np.sqrt((px - com[0])**2 + (py - com[1])**2 + (pz - com[2])**2)
        phi = multipole_phi(M, np.zeros(6), com, px, py, pz)
        worst = max(worst, abs(phi - (-G * M / r)) / abs(G * M / r))
    # radial probes from the COM so r = 10 and 20 exactly (dy=dz=0).
    r1 = multipole_phi(M, np.zeros(6), com, com[0] + 10.0, com[1], com[2])
    r2 = multipole_phi(M, np.zeros(6), com, com[0] + 20.0, com[1], com[2])
    scale_ok = abs(r1 / r2 - 2.0) < 1e-12
    g1 = worst < 1e-13 and scale_ok
    print(f"  worst |Phi - (-GM/r)|/|GM/r| = {worst:.2e} ; 1/r scaling ok = {scale_ok}  "
          f"{'PASS' if g1 else 'FAIL'}")
    fails += not g1

    # Gate 2: quadrupole angular structure + 1/r^3 scaling.
    print("--- GATE 2: quadrupole P_2 structure + 1/r^3 scaling ---")
    Qd = np.array([-1.0, -1.0, 2.0, 0.0, 0.0, 0.0])  # diagonal traceless (prolate along z)
    R = 10.0
    phz = multipole_phi(0.0, Qd, com, com[0], com[1], com[2] + R)  # along z: dQd=Qzz R^2
    phx = multipole_phi(0.0, Qd, com, com[0] + R, com[1], com[2])  # along x: dQd=Qxx R^2
    # analytic: Phi_z = -0.5 G Qzz / R^3 ; Phi_x = -0.5 G Qxx / R^3
    az = -0.5 * G * Qd[2] / R**3
    ax = -0.5 * G * Qd[0] / R**3
    ang_ok = abs(phz - az) / abs(az) < 1e-13 and abs(phx - ax) / abs(ax) < 1e-13
    # 1/r^3 scaling along z
    phz2 = multipole_phi(0.0, Qd, com, com[0], com[1], com[2] + 2 * R)
    scale3_ok = abs(phz / phz2 - 8.0) < 1e-10
    g2 = ang_ok and scale3_ok
    print(f"  Phi_z/Phi_x = {phz/phx:.6f} (=Qzz/Qxx={Qd[2]/Qd[0]:.1f}); 1/r^3 scaling ok = "
          f"{scale3_ok}  {'PASS' if g2 else 'FAIL'}")
    fails += not g2

    # Gate 3: superposition (linear in moments).
    print("--- GATE 3: superposition Phi(M,Q) == Phi(M,0) + Phi(0,Q) ---")
    worst_sup = 0.0
    for (px, py, pz) in [(6, 1, 2), (-4, 5, -7), (9, -3, 4)]:
        full = multipole_phi(M, Q, com, px, py, pz)
        parts = (multipole_phi(M, np.zeros(6), com, px, py, pz) +
                 multipole_phi(0.0, Q, com, px, py, pz))
        worst_sup = max(worst_sup, abs(full - parts) / (abs(full) + 1e-300))
    g3 = worst_sup < 1e-14
    print(f"  worst superposition residual = {worst_sup:.2e}  {'PASS' if g3 else 'FAIL'}")
    fails += not g3

    # Gate 4: traceless-Q shell average == 0.
    print("--- GATE 4: traceless quadrupole shell-average == 0 ---")
    R = 12.0
    ntheta, nphi = 200, 400
    th = (np.arange(ntheta) + 0.5) * np.pi / ntheta
    ph = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    acc = 0.0
    wsum = 0.0
    for t in th:
        w = np.sin(t)  # solid-angle weight
        for p in ph:
            x = com[0] + R * np.sin(t) * np.cos(p)
            y = com[1] + R * np.sin(t) * np.sin(p)
            zc = com[2] + R * np.cos(t)
            acc += w * multipole_phi(0.0, Q, com, x, y, zc)
            wsum += w
    avg = acc / wsum
    # normalize by the typical |Phi_quad| magnitude at R
    typ = abs(G * 0.5 * np.max(np.abs(Q)) / R**3)
    g4 = abs(avg) / typ < 1e-3
    print(f"  |shell-avg Phi_quad| / typical = {abs(avg)/typ:.2e}  (quadrature-limited)  "
          f"{'PASS' if g4 else 'FAIL'}")
    fails += not g4

    print(f"\n{'ALL GATES PASS' if fails == 0 else 'GATE FAILURES'} ({fails} failures)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
