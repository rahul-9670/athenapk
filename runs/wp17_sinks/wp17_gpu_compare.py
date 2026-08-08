#!/usr/bin/env python
"""WS-1 inc6 — GPU sink accretion rate against the CPU-validated WP-17 answer.

Same physics constants and the same measurement as wp17_mdot.py (shell flux at r = 1.5, outside
r_acc and inside r_B, normalised against the INSTANTANEOUS Mdot_B(M(t)) built from the dumped
sink mass). The functions are duplicated rather than imported because wp17_mdot.py prints and
loops at import time.

PASS CRITERION. The CPU N=128 answer is Mdot/Mdot_B = 1.01 +/- 0.02 over t = 4.0-6.0. The GPU leg
passes if its plateau agrees within the CPU leg's own uncertainty band, i.e. the two ratios differ
by less than ~0.05. A GPU value that is merely "close to 1" is NOT sufficient on its own -- the
comparison is against the CPU MEASUREMENT, because that is the number the project has validated.
"""
import os, glob
import numpy as np
import h5py

G, C_INF, GAMMA = 1.0, 1.0, 1.4
RHO_INF = 1.0e-3
lam = (0.5) ** ((GAMMA + 1) / (2 * (GAMMA - 1))) * \
      ((5 - 3 * GAMMA) / 4) ** (-(5 - 3 * GAMMA) / (2 * (GAMMA - 1)))
H = "/beegfs/u/bbg6470/athenapk/runs/wp17_sinks"
R_MEAS = 1.5


def mdot_b(M):
    return 4 * np.pi * lam * (G * M) ** 2 * RHO_INF / C_INF ** 3


def read(fn, r0=R_MEAS):
    with h5py.File(fn, "r") as f:
        t = float(f["Info"].attrs["Time"])
        prim = f["prim"][:]
        xc, yc, zc = (f["VolumeLocations/x"][:], f["VolumeLocations/y"][:],
                      f["VolumeLocations/z"][:])
        dx = float(xc[0, 1] - xc[0, 0])
        rho = prim[:, 0]
        X = xc[:, None, None, :] * np.ones_like(rho)
        Y = yc[:, None, :, None] * np.ones_like(rho)
        Z = zc[:, :, None, None] * np.ones_like(rho)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        vr = (prim[:, 1] * X + prim[:, 2] * Y + prim[:, 3] * Z) / np.maximum(r, 1e-30)
        dr = max(0.25, 2 * dx)
        sh = np.abs(r - r0) < dr
        mdot = -np.sum(rho[sh] * vr[sh] * dx ** 3) / (2 * dr)
        sv = f["sinks/SwarmVars"]
        M = float(np.array(sv["mass"])[0]) if "mass" in sv else np.nan
        return t, mdot, M, dx


def leg(subdir, label):
    fs = sorted(glob.glob(os.path.join(H, subdir, "parthenon.out0.*.phdf")))
    if not fs:
        print(f"\n=== {label} ({subdir}): NO DUMPS ===")
        return None
    print(f"\n=== {label}  ({subdir}, {len(fs)} dumps) ===")
    print("     t        Mdot      M_sink    Mdot_B(M)    ratio")
    rows = []
    for fn in fs:
        t, md, M, dx = read(fn)
        mb = mdot_b(M)
        ratio = md / mb if mb > 0 else np.nan
        rows.append((t, md, M, mb, ratio))
        print(f"  {t:6.3f}  {md:9.4f}  {M:9.5f}  {mb:9.4f}   {ratio:7.4f}")
    # Plateau window t = 4.0-6.0, the same window the CPU answer was quoted over.
    pl = [r[4] for r in rows if 4.0 <= r[0] <= 6.0 and np.isfinite(r[4])]
    if not pl:
        print("  no samples in the t = 4.0-6.0 plateau window")
        return None
    m, s = float(np.mean(pl)), float(np.std(pl))
    print(f"  plateau (t=4.0-6.0, n={len(pl)}): Mdot/Mdot_B = {m:.4f} +/- {s:.4f}")
    return m, s


def main():
    cpu = leg("b128", "CPU N=128 (WP-17, job 2453375)")
    gpu = leg("g128", "GPU N=128 (WS-1 inc6)")

    print("\n================ VERDICT ================")
    if gpu is None:
        print("  INCONCLUSIVE — the GPU leg produced no usable plateau.")
    elif cpu is None:
        print(f"  GPU plateau = {gpu[0]:.4f} +/- {gpu[1]:.4f}")
        print("  CPU dumps are not on disk; compare against the RECORDED CPU answer 1.01 +/- 0.02.")
        d = abs(gpu[0] - 1.01)
        print(f"  |GPU - 1.01| = {d:.4f}  ->  {'PASS' if d < 0.05 else 'FAIL'}")
    else:
        d = abs(gpu[0] - cpu[0])
        print(f"  CPU plateau = {cpu[0]:.4f} +/- {cpu[1]:.4f}")
        print(f"  GPU plateau = {gpu[0]:.4f} +/- {gpu[1]:.4f}")
        print(f"  |GPU - CPU| = {d:.4f}  ->  {'PASS' if d < 0.05 else 'FAIL'}")
    print("=========================================")


if __name__ == "__main__":
    main()
