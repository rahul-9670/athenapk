#!/usr/bin/env python
"""AUDIT N1: compare a 2D CT leg against its z-invariant 3D reference.

The whistler eigenmode (iprob=60) is uniform in y and z, so the 3D solution stays z-invariant
and its k-slices are all equal. Collapsing the 3D field over k and comparing against the single
2D layer is therefore an exact comparison, not an approximation.

Tolerance. The two legs are NOT expected to be bit-identical: 2D and 3D take different code
paths (different kernels, different reduction/loop extents), and floating-point addition is not
associative. They ARE expected to agree to a few ULP accumulated over the run. The threshold
below is on the perturbation amplitude (amp = 1e-6), not on the total field, because the whole
signal lives in the perturbation -- normalising by the O(1) guide field B_x would hide a
100%-wrong B_z.
"""
import sys, glob, os
import numpy as np
import h5py

TOL = 1.0e-8  # relative to the perturbation amplitude amp = 1e-6


def load(d):
    fs = sorted(glob.glob(os.path.join(d, "*.prim.final.phdf")))
    if not fs:
        fs = sorted(glob.glob(os.path.join(d, "*.phdf")))
    if not fs:
        raise SystemExit(f"no phdf in {d}")
    with h5py.File(fs[-1], "r") as f:
        prim = np.asarray(f["prim"][()])  # [block, comp, k, j, i]
        names = [n.decode() if isinstance(n, bytes) else n
                 for n in f["Info"].attrs["ComponentNames"]]
    return prim, names, fs[-1]


def main():
    d2, d3, label = sys.argv[1], sys.argv[2], sys.argv[3]
    p2, n2, f2 = load(d2)
    p3, n3, f3 = load(d3)

    # Collapse the k axis. In 3D every k-slice must already be identical; verify that first --
    # if it is not, the reference itself is not z-invariant and the comparison is meaningless.
    kspread = np.abs(p3 - p3[:, :, :1, :, :]).max()
    amp = 1.0e-6
    if kspread / amp > 1.0e-9:
        print(f"  {label}: REFERENCE NOT z-INVARIANT (max k-spread/amp = {kspread/amp:.3e})"
              " -- comparison void")
        return 1
    a3 = p3[:, :, :1, :, :]
    a2 = p2[:, :, :1, :, :]
    if a2.shape != a3.shape:
        print(f"  {label}: SHAPE MISMATCH 2D {a2.shape} vs 3D {a3.shape}")
        return 1

    worst, worst_name = 0.0, ""
    lines = []
    ncomp = a2.shape[1]
    for c in range(ncomp):
        d = np.abs(a2[:, c] - a3[:, c]).max() / amp
        nm = n2[c] if c < len(n2) else f"comp{c}"
        lines.append(f"{nm}={d:.2e}")
        if d > worst:
            worst, worst_name = d, nm
    print(f"  {label}: max |2D-3D|/amp = {worst:.3e} (worst: {worst_name})")
    print("      " + "  ".join(lines))
    if worst > TOL:
        print(f"      -> FAIL (tol {TOL:.1e})")
        return 1
    print(f"      -> PASS (tol {TOL:.1e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
