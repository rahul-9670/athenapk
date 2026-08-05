#!/usr/bin/env python
"""AUDIT N1 FALSIFIER: prove the 2D B_z is actually evolving, not merely agreeing.

Two frozen fields agree perfectly. If the N1 fix were inert -- if B_z in 2D were still pinned to
its initial value -- every leg of compare_2d3d.py would still pass, because the 3D reference at
these amplitudes barely moves either. This script closes that hole by measuring how far B_z
travels in the 2D run and comparing it against what the whistler eigenmode requires.

For iprob=60 the mode is circularly polarized: B_y = amp cos(kx), B_z = h amp sin(kx) at t=0,
and the pair ROTATES at the whistler frequency. Over the run the phase advances by omega*t, so
B_z must change by O(amp) -- not O(amp * 1e-12). A near-zero result here means B_z is frozen and
the fix did not take.
"""
import sys, glob, os
import numpy as np
import h5py

AMP = 1.0e-6


def snap(path):
    with h5py.File(path, "r") as f:
        prim = np.asarray(f["prim"][()])
        names = [n.decode() if isinstance(n, bytes) else n
                 for n in f["Info"].attrs["ComponentNames"]]
        t = float(f["Info"].attrs["Time"])
    return prim, names, t


def main():
    d = sys.argv[1]
    fs = sorted(glob.glob(os.path.join(d, "*.phdf")))
    if len(fs) < 2:
        print(f"  need >=2 outputs in {d}, found {len(fs)} -- cannot measure motion")
        return 1
    p0, names, t0 = snap(fs[0])
    p1, _, t1 = snap(fs[-1])
    # prim component order: 0=rho 1-3=v 4=P 5-7=B (see CLAUDE.md); find B3 by name if possible.
    try:
        ib3 = [i for i, n in enumerate(names) if n.endswith("B3") or n == "prim_B3"][0]
        iv3 = [i for i, n in enumerate(names) if n.endswith("V3") or n == "prim_v3"][0]
    except IndexError:
        ib3, iv3 = 7, 3
    dbz = np.abs(p1[:, ib3] - p0[:, ib3]).max() / AMP
    dvz = np.abs(p1[:, iv3] - p0[:, iv3]).max() / AMP
    print(f"  t: {t0:.6g} -> {t1:.6g}")
    print(f"  max|dB_z|/amp = {dbz:.4e}   max|dv_z|/amp = {dvz:.4e}")
    if dbz < 1.0e-6:
        print("      -> FAIL: B_z did not move. The 2D CT fix is INERT and every"
              " 2D-vs-3D pass above is passing for the wrong reason.")
        return 1
    print("      -> PASS: B_z is a live degree of freedom in 2D.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
