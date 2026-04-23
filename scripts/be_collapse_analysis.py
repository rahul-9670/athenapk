#!/usr/bin/env python3
"""
Bonnor-Ebert collapse diagnostic: rho_max, phi_min, nblocks vs time.

Usage: python3 be_collapse_analysis.py <prefix>
  where <prefix> is e.g. 'parthenon.out0'
"""
import sys
import glob
import h5py
import numpy as np


def analyze(prefix):
    files = sorted(glob.glob(f"{prefix}.?????.phdf"))
    if not files:
        print(f"ERROR: no files matching {prefix}.?????.phdf found")
        return

    print(f"Found {len(files)} snapshots")
    print(f"{'t':>10} {'nblocks':>8} {'rho_max':>12} {'rho_min':>12} "
          f"{'phi_min':>12} {'phi_max':>12}")
    print("-" * 72)

    for fn in files:
        with h5py.File(fn, "r") as f:
            t = f["Info"].attrs["Time"]
            # prim shape: (nblocks, nvars, nz, ny, nx)
            prim = f["prim"][:]
            rho = prim[:, 0, ...]            # first variable is density
            phi = f["grav.phi"][:]

            rho_max = rho.max()
            rho_min = rho.min()
            phi_min = phi.min()
            phi_max = phi.max()
            nblocks = rho.shape[0]

            print(f"{t:10.4f} {nblocks:8d} {rho_max:12.4e} {rho_min:12.4e} "
                  f"{phi_min:12.4e} {phi_max:12.4e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 be_collapse_analysis.py <prefix>")
        print("  e.g.  python3 be_collapse_analysis.py parthenon.out0")
        sys.exit(1)
    analyze(sys.argv[1])
