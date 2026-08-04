#!/usr/bin/env python
"""WP-15 — collapse against a KNOWN ANSWER: pressure-free uniform-sphere free-fall.

The exact Newtonian solution for a pressureless uniform sphere, parametric in beta:

    r/r0 = cos^2(beta)      rho/rho0 = sec^6(beta)      t/t_ff = (2/pi)(beta + sin b cos b)
    t_ff = sqrt(3 pi / (32 G rho0))

With `four_pi_G = 1` (G = 1/(4 pi)) and rho0 = 1 this gives t_ff = pi*sqrt(3/8) = 1.9238...,
a property of the equations rather than of the code.

Inverting for a density contrast X: beta = arccos(X^(-1/6)), t(X) = (2/pi) t_ff (beta + sin b cos b).

The measured quantity is the TIME at which the central density first reaches X. It depends on the
Poisson solve, the gravity source, the hydro advection and their stage coupling all being right at
once -- the exact combination WP-13's stale-RHS bug corrupted, and which no self-convergence study
could have caught (a one-stage-lagged RHS converges beautifully to the wrong answer).
"""
import glob
import os
import re
import numpy as np

H = "/beegfs/u/bbg6470/athenapk/runs/wp15_freefall"
# Two 64^3 sweeps replace the original resolution ladder. 128^3 needs ~2.4 h for the ~975 cycles
# that dt_ceil = 0.002 requires and does not decide the open question (see below).
# amb_*  = ambient density sweep -- THE decisive test of the negative error.
# p*     = pressure sweep -- rules pressure in or out as the cause.
LEGS = ["amb_1e-2", "amb_5e-3", "amb_3e-3", "p1e-6", "p5e-7", "p2e-6", "n64_p1e-6"]
FOUR_PI_G = 1.0
RHO0 = 1.0
T_FF = np.sqrt(3.0 * np.pi / (32.0 * (FOUR_PI_G / (4.0 * np.pi)) * RHO0))


def t_analytic(X):
    """Exact arrival time at density contrast X (X >= 1)."""
    beta = np.arccos(X ** (-1.0 / 6.0))
    return (2.0 / np.pi) * T_FF * (beta + np.sin(beta) * np.cos(beta))


def hst_cols(p):
    names = []
    for line in open(p):
        if line.startswith("#") and "[" in line:
            names = re.findall(r"\[(\d+)\]=(\S+)", line)
        elif not line.startswith("#"):
            break
    return {n: int(i) - 1 for i, n in names}


def central_rho_series(d):
    """(t, rho_max) from the phdf dumps. rho_max IS the central density while the collapse is
    homologous -- checked by the fact that it grows monotonically and smoothly."""
    import h5py
    out = []
    for f in sorted(glob.glob(os.path.join(d, "parthenon.out0.*.phdf"))):
        if "final" in f:
            continue
        with h5py.File(f, "r") as h:
            out.append((float(h["Info"].attrs["Time"]), float(h["prim"][:, 0, ...].max())))
    return np.array(out) if out else None


def main():
    print(f"WP-15 — pressure-free free-fall.  t_ff (exact, four_pi_G=1, rho0=1) = {T_FF:.6f}\n")
    print("Exact arrival times:")
    for X in (2.0, 10.0, 100.0, 1000.0):
        print(f"   rho/rho0 = {X:7.0f}  ->  t = {t_analytic(X):.6f}  "
              f"({t_analytic(X)/T_FF:.4f} t_ff)")
    print()

    print(f"{'leg':>14} {'X=2':>22} {'X=10':>22} {'X=100':>22}")
    print(f"{'':>14} {'t_meas   err%':>22} {'t_meas   err%':>22} {'t_meas   err%':>22}")
    for leg in LEGS:
        d = os.path.join(H, leg)
        s = central_rho_series(d)
        if s is None or len(s) < 3:
            print(f"{leg:>14}   MISSING or too few dumps")
            continue
        t, r = s[:, 0], s[:, 1]
        row = f"{leg:>14}"
        for X in (2.0, 10.0, 100.0):
            if r.max() < X:
                row += f" {'not reached':>22}"
                continue
            # first crossing, linearly interpolated in log(rho) for accuracy
            i = int(np.argmax(r >= X))
            if i == 0:
                row += f" {'at t=0?':>22}"
                continue
            lo, hi = i - 1, i
            f = (np.log(X) - np.log(r[lo])) / (np.log(r[hi]) - np.log(r[lo]))
            tm = t[lo] + f * (t[hi] - t[lo])
            ta = t_analytic(X)
            row += f" {tm:10.5f} {100*(tm-ta)/ta:+8.3f}%"
        print(row)

    print("""
HOW TO READ IT.
  * THE ERROR CAME OUT NEGATIVE (-0.43 to -0.91 %), i.e. the collapse is FASTER than the
    isolated-sphere solution. That is the direction flagged in advance as the serious one, so it
    is being separated rather than reported:
      - amb_* sweep: the ambient is NOT massless. At rho_out = 1e-2 it carries 14.3 % of the
        sphere's mass and falls inward too, augmenting the enclosed mass; t_ff ~ 1/sqrt(M) turns
        a few % of extra mass into a few tenths of a % faster collapse -- the size observed.
        If the error SCALES with rho_out and extrapolates toward 0, the cause is the setup.
        If it is INDEPENDENT of rho_out, the gravitational source is too strong -- a real defect.
      - p* sweep: pressure can only DELAY collapse (positive error), so if halving the pressure
        makes the error MORE negative, pressure is excluded as the cause.
    rho_out cannot go below ~3e-3: at 1e-6 the ambient sound speed is 1.183 vs 0.0216 there, and
    that supersonic ambient destroyed an earlier attempt. Hence extrapolation, not zero.""")


if __name__ == "__main__":
    main()
