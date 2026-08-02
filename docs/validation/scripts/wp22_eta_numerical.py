#!/usr/bin/env python
"""WP-22 — bound the NUMERICAL resistivity from the WP-14 Alfven ladder, no new compute.

An Alfven wave run with ZERO physical resistivity decays only through the scheme's own
dissipation:  B_perp(t) = B_perp(0) * exp(-eta_num k^2 t / 2).  Over one full period T the
transverse-field error is therefore

    L1(B_perp) ~= A0 * (1 - exp(-eta_num k^2 T/2)) ~= A0 * eta_num k^2 T / 2      (small dissip.)
 => eta_num ~= 2 * (L1/A0) / (k^2 T)

L1 also contains DISPERSION error, which does not damp the wave, so this over-attributes error to
dissipation: the number below is a rigorous UPPER BOUND on eta_num. That is the conservative
direction for WP-22 -- if even the bound sits below the physical eta_A in the first core, the
fossil-field argument is safe.

Geometry follows src/pgen/linear_wave_mhd.cpp exactly (angles default to -999.9 => computed).
Background state from that file: d0 = 1, bx0 = 1 => v_A = bx0/sqrt(d0) = 1 in code units
(AthenaPK is Heaviside-Lorentz: v_A = B/sqrt(rho), NO 4pi -- CLAUDE.md's #1 analysis pitfall).
"""
import os
import numpy as np

HERE = "/beegfs/u/bbg6470/athenapk/runs/wp14_order"
AMP = 1.0e-6
X1, X2, X3 = 3.0, 1.5, 1.5      # box sizes from runs/wp14_order/lw_mhd.in
COL_B2C = 11                     # B2c_L1 in linearwave-errors.dat


def wavelength():
    ang_3 = np.arctan(X1 / X2)
    ang_2 = np.arctan(0.5 * (X1 * np.cos(ang_3) + X2 * np.sin(ang_3)) / X3)
    x1 = X1 * np.cos(ang_2) * np.cos(ang_3)
    x2 = X2 * np.cos(ang_2) * np.sin(ang_3)
    x3 = X3 * np.sin(ang_2)
    return min(x1, x2, x3)


def read(path):
    with open(path) as f:
        rows = [l.split() for l in f if not l.startswith("#") and l.strip()]
    return [float(x) for x in rows[-1]]


def main():
    lam = wavelength()
    k = 2.0 * np.pi / lam
    vA = 1.0                       # bx0=1, d0=1, Heaviside-Lorentz
    T = lam / vA                   # one period (tlim = 1.0 period in the pgen's units)
    print(f"Alfven ladder: lambda = {lam:.6f}, k = {k:.6f}, v_A = {vA}, period T = {T:.6f}\n")
    print(f"{'N':>5} {'dx':>10} {'L1(B2c)':>12} {'eta_num (bound)':>17} {'eta/(v_A dx)':>14}")
    rows = []
    for n in (16, 32, 64, 128):
        p = os.path.join(HERE, f"wf1_n{n}", "linearwave-errors.dat")
        if not os.path.exists(p):
            continue
        L1 = read(p)[COL_B2C]
        dx = X2 / n                            # nx2 = nx3 = N; cell size
        eta = 2.0 * (L1 / AMP) / (k * k * T)
        rows.append((n, dx, L1, eta))
        print(f"{n:5d} {dx:10.5f} {L1:12.4e} {eta:17.4e} {eta/(vA*dx):14.4f}")
    if len(rows) >= 2:
        n = np.array([r[0] for r in rows], float)
        e = np.array([r[3] for r in rows])
        p = np.polyfit(np.log(1.0 / n), np.log(e), 1)[0]
        print(f"\n  eta_num ~ dx^{p:.2f}   (2nd-order scheme => expect ~2; 1st order => ~1)")
        C = rows[-1][3] / (vA * rows[-1][1])
        print(f"  finest rung: eta_num = {C:.4f} * v_A * dx   [dimensionless coefficient]")
        print("\n  To apply at production: eta_num(core) ~ C_eff * v * dx_core, with the SAME")
        print("  dx-scaling measured above. Compare against eta_A from the ionization model")
        print("  evaluated at first-core rho/T/B -- production uses ambipolar_coeff =")
        print("  ionization_chem, so eta_A is NOT a constant and must be recomputed offline.")


if __name__ == "__main__":
    main()
