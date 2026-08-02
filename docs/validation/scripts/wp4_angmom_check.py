#!/usr/bin/env python
"""WP-4 -- independent cross-check of the am-L* history columns.

Recomputes L = int rho (r x v) dV directly from a .phdf snapshot, in numpy, with no
shared code path with the Kokkos kernel in src/diagnostics/angmom_diag.cpp. Agreement
between the two validates the KERNEL. It says nothing about whether the SCHEME conserves
L -- that is a separate measurement (see wp4_drift below).

Moments are about the centre of the box, matching angmom_diag.cpp.

Usage:
    wp4_angmom_check.py <snapshot.phdf> [--hst <file.hst> --time <t>]
    wp4_angmom_check.py --drift <file.hst>
"""
import argparse
import sys

import h5py
import numpy as np


def angmom_from_phdf(path):
    """Return (Lx, Ly, Lz, mass) computed from a snapshot, box-centre moments."""
    with h5py.File(path, "r") as f:
        prim = f["prim"][...]  # [block, comp, k, j, i]
        # Cell-centre coordinates per block, and cell faces for the widths.
        xc = f["VolumeLocations"]["x"][...]  # [block, i]
        yc = f["VolumeLocations"]["y"][...]
        zc = f["VolumeLocations"]["z"][...]
        xf = f["Locations"]["x"][...]  # [block, i+1]
        yf = f["Locations"]["y"][...]
        zf = f["Locations"]["z"][...]
        dom = f["Info"].attrs["RootGridDomain"]  # [x0,x1,_, y0,y1,_, z0,z1,_]

    # Box centre -- the same origin angmom_diag.cpp uses.
    x0 = 0.5 * (dom[0] + dom[1])
    y0 = 0.5 * (dom[3] + dom[4])
    z0 = 0.5 * (dom[6] + dom[7])

    rho = prim[:, 0]  # [b,k,j,i]
    vx, vy, vz = prim[:, 1], prim[:, 2], prim[:, 3]

    dx = np.diff(xf, axis=1)  # [b,i]
    dy = np.diff(yf, axis=1)
    dz = np.diff(zf, axis=1)
    # dV[b,k,j,i] = dz[b,k]*dy[b,j]*dx[b,i]
    dV = dz[:, :, None, None] * dy[:, None, :, None] * dx[:, None, None, :]

    X = (xc - x0)[:, None, None, :]
    Y = (yc - y0)[:, None, :, None]
    Z = (zc - z0)[:, :, None, None]

    Lx = np.sum(rho * (Y * vz - Z * vy) * dV, dtype=np.float64)
    Ly = np.sum(rho * (Z * vx - X * vz) * dV, dtype=np.float64)
    Lz = np.sum(rho * (X * vy - Y * vx) * dV, dtype=np.float64)
    mass = np.sum(rho * dV, dtype=np.float64)
    return Lx, Ly, Lz, mass


def read_hst(path):
    """Parse a parthenon .hst into {column_name: array}. Indices come from the HEADER --
    never hardcode them; enabling any optional diagnostic shifts every later column."""
    names = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") and "=" in line and "cycle" in line:
                # e.g. "# [1]=time     [2]=dt   [3]=mass ..."
                parts = line.lstrip("#").split()
                names = []
                for p in parts:
                    if "]=" in p:
                        names.append(p.split("]=", 1)[1])
                    elif names:
                        names[-1] += "_" + p
                break
    if names is None:
        raise SystemExit(f"no column header found in {path}")
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != len(names):
        raise SystemExit(
            f"{path}: header has {len(names)} names but {data.shape[1]} columns"
        )
    return {n: data[:, i] for i, n in enumerate(names)}, names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("snapshot", nargs="?")
    ap.add_argument("--hst")
    ap.add_argument("--row", type=int, default=0, help="hst row to compare against")
    ap.add_argument("--drift", help="report L drift over a .hst and stop")
    args = ap.parse_args()

    if args.drift:
        cols, names = read_hst(args.drift)
        print(f"columns: {names}")
        for c in ("am-Lx", "am-Ly", "am-Lz"):
            if c not in cols:
                print(f"  {c}: ABSENT")
                continue
            v = cols[c]
            v0 = v[0]
            # Normalise the drift by the largest |L| in the run, not by L(0): for a run
            # whose L(0) happens to be ~0 (Lx, Ly here) a relative-to-initial number is a
            # near-zero denominator and is meaningless. Same trap as WP-5's max|rel|=5.35.
            scale = np.max(np.abs(v)) if np.max(np.abs(v)) > 0 else 1.0
            print(
                f"  {c}: L0={v0:.9e}  Lend={v[-1]:.9e}  "
                f"max|dL|/max|L| = {np.max(np.abs(v - v0))/scale:.3e}"
            )
        for c in ("am-FLx", "am-FLy", "am-FLz", "am-Tmagx", "am-Tmagy", "am-Tmagz"):
            if c in cols:
                print(f"  {c}: max|.| = {np.max(np.abs(cols[c])):.6e}")
        return

    if not args.snapshot:
        raise SystemExit("need a snapshot (or --drift)")

    Lx, Ly, Lz, mass = angmom_from_phdf(args.snapshot)
    print(f"snapshot {args.snapshot}")
    print(f"  mass  = {mass:.12e}")
    print(f"  Lx    = {Lx:.12e}")
    print(f"  Ly    = {Ly:.12e}")
    print(f"  Lz    = {Lz:.12e}")

    if args.hst:
        cols, _ = read_hst(args.hst)
        missing = [c for c in ("am-Lx", "am-Ly", "am-Lz") if c not in cols]
        if missing:
            raise SystemExit(f"hst lacks {missing} -- was hydro/angmom_diag on?")
        r = args.row
        print(f"\n  hst row {r}, t = {cols['time'][r]:.6f}")
        for name, val in (("am-Lx", Lx), ("am-Ly", Ly), ("am-Lz", Lz)):
            h = cols[name][r]
            # The snapshot is float32; the .hst prints ~6 significant figures. The
            # achievable agreement is set by whichever is coarser, NOT by float64 eps.
            den = max(abs(h), abs(val))
            rel = abs(h - val) / den if den > 0 else 0.0
            print(f"  {name}: hst={h:.6e}  numpy={val:.6e}  rel={rel:.3e}")


if __name__ == "__main__":
    main()
