#!/usr/bin/env python
"""Flux-retention table for the fossil-field paper (ms.tex "headline" pending box).

For a tier run, at matched milestones rho_max/rhocrit in {1, 1e2, 1e5, [deepest]},
measure:
  (a) Phi_core : magnetic flux through the equatorial plane (z = z_peak) within the
      FHC core radius r_core (spherical radius where the mass-weighted density
      profile about the peak first drops below rhocrit);
  (b) mu_core  : core mass-to-flux ratio in units of the critical value, reported
      for BOTH common conventions:
        NN78: (M/Phi)_crit = 1/(2 pi sqrt(G))          (Nakano & Nakamura 1978)
        MS76: (M/Phi)_crit = 0.53/(3 pi) sqrt(5/G)     (Mouschovias & Spitzer 1976)
  (c) retention = Phi_core / Phi_0(M_core), where Phi_0 is the INITIAL flux
      threading the same Lagrangian mass: from the t=0 snapshot, find r0 such
      that M(<r0) = M_core about the box center, then Phi_0 = flux through the
      midplane disc of radius r0 at t=0 (numerically; = pi r0^2 B0z since B is
      uniform at t=0).

Units: shared FHC code units (rho0=5.467e-19, v0=1.9e4, l0=2.81e16, HL B with
B_unit = sqrt(4 pi rho0 v0^2) = 4.98e-5 G). All ratios are convention-free except
mu, which is why both critical constants are quoted. Raw h5py (no yt): phdf
layout prim[block, comp, k, j, i], comps 0=rho, 5-7=B; VolumeLocations = cell
centers; Locations = faces (33 per 32-cell block edge).

Usage:
  flux_retention.py <run_dir> [--scan-stride 8] [--out results.json]
"""
import argparse
import glob
import json
import os
import sys

import h5py
import numpy as np

# shared FHC code units (identical across all tiers by construction)
RHO0 = 5.467e-19          # g/cm^3
V0 = 1.9e4                # cm/s
L0 = 2.81e16              # cm
M0 = RHO0 * L0**3         # g
B_UNIT = np.sqrt(4.0 * np.pi * RHO0 * V0 * V0)  # G (Heaviside-Lorentz B-unit)
G_CGS = 6.674e-8
MSUN = 1.989e33
RHOCRIT_CGS = 1.0e-13
RHOCRIT_CODE = RHOCRIT_CGS / RHO0
# critical mass-to-flux (cgs)
M2F_CRIT_NN78 = 1.0 / (2.0 * np.pi * np.sqrt(G_CGS))
M2F_CRIT_MS76 = (0.53 / (3.0 * np.pi)) * np.sqrt(5.0 / G_CGS)


def rho_max_code(path):
    """Global max code density of a snapshot (reads only the density component)."""
    with h5py.File(path, "r") as f:
        rho = f["prim"][:, 0, :, :, :]
        return float(rho.max()), float(f["Info"].attrs["Time"])


def load_geom(path):
    """rho, B (code), cell centers, per-block dx (code) from a phdf file."""
    with h5py.File(path, "r") as f:
        prim = f["prim"]
        rho = prim[:, 0, :, :, :].astype(np.float64)
        bx = prim[:, 5, :, :, :].astype(np.float64)
        by = prim[:, 6, :, :, :].astype(np.float64)
        bz = prim[:, 7, :, :, :].astype(np.float64)
        xc = f["VolumeLocations/x"][:]  # (nb, nx) cell centers
        yc = f["VolumeLocations/y"][:]
        zc = f["VolumeLocations/z"][:]
        t = float(f["Info"].attrs["Time"])
    # per-block cell size (blocks are uniform-dx cubes)
    dxb = xc[:, 1] - xc[:, 0]
    return rho, bx, by, bz, xc, yc, zc, dxb, t


def cell_grids(xc, yc, zc, shape):
    """Broadcast per-block 1D center arrays to full (nb,nz,ny,nx) grids."""
    nb, nz, ny, nx = shape
    X = xc[:, None, None, :] * np.ones((1, nz, ny, 1))
    Y = yc[:, None, :, None] * np.ones((1, nz, 1, nx))
    Z = zc[:, :, None, None] * np.ones((1, 1, ny, nx))
    return X, Y, Z


def measure_snapshot(path, r_core_fixed=None, center=None):
    """Core radius/mass/flux about the density peak (or a supplied center).

    r_core_fixed: skip the profile and use this radius (code units) -- used for
    the t=0 Lagrangian inversion. Returns a dict of code-unit quantities.
    """
    rho, bx, by, bz, xc, yc, zc, dxb, t = load_geom(path)
    nb, nz, ny, nx = rho.shape
    X, Y, Z = cell_grids(xc, yc, zc, rho.shape)
    dV = (dxb**3)[:, None, None, None] * np.ones_like(rho)
    dA = (dxb**2)[:, None, None, None] * np.ones_like(rho)
    dz_half = 0.5 * dxb[:, None, None, None] * np.ones_like(rho)

    if center is None:
        pk = np.unravel_index(np.argmax(rho), rho.shape)
        cxyz = (X[pk], Y[pk], Z[pk])
    else:
        cxyz = center
    r = np.sqrt((X - cxyz[0]) ** 2 + (Y - cxyz[1]) ** 2 + (Z - cxyz[2]) ** 2)

    if r_core_fixed is None:
        # mass-weighted spherical profile about the peak; r_core = first bin
        # whose mean density drops below rhocrit (linear interp in log rho).
        rmax = 6.0  # code; generous upper bound on any FHC radius (~1e4 AU)
        nbins = 240
        rb = np.logspace(np.log10(2e-4), np.log10(rmax), nbins)
        idx = np.digitize(r.ravel(), rb)
        m = (rho * dV).ravel()
        v = dV.ravel()
        msum = np.bincount(idx, weights=m, minlength=nbins + 1)
        vsum = np.bincount(idx, weights=v, minlength=nbins + 1)
        prof = np.where(vsum > 0, msum / np.maximum(vsum, 1e-300), 0.0)[1:nbins]
        rmid = 0.5 * (rb[:-1] + rb[1:])[: nbins - 1]
        below = np.where((prof > 0) & (prof < RHOCRIT_CODE))[0]
        # Bins that are BOTH non-empty and above rhocrit -- the only valid inner anchors for the
        # log-linear crossing. `below` already filters empty bins; the inner anchor must be
        # filtered the same way, and previously was not (see the bug note below).
        above = np.where(prof > RHOCRIT_CODE)[0]
        if len(below) == 0 or below[0] == 0 or len(above) == 0:
            r_core = rmid[0]
        else:
            i1 = below[0]
            # BUG FIXED 2026-08-06. This used to be `i0 = i1 - 1` with no check that the
            # preceding bin was populated. The radial bins are log-spaced from 2e-4, so at the
            # matched epoch the inner ones are far finer than the cells and many are EMPTY
            # (prof == 0). An empty i0 gives log(0) = -inf -> f = nan -> r_core = nan, and then
            # BOTH `r < r_core` and `rcyl < r_core` are all-False, so M_core and Phi_core come
            # back exactly 0.0 with no error raised.
            #
            # This was not a rare edge case, it was systematic AT THE MEASUREMENT EPOCH:
            # RHOCRIT_CODE = 1e-13/5.467e-19 = 1.829e5 is exactly the rho_max the ensemble
            # matches on, so the super-critical region is then only a few cells across and the
            # inner bins are empty by construction. Six of seven finished ensemble members were
            # silently dropped from the flux-retention distribution this way (2026-08-06).
            cand = above[above < i1]
            if len(cand) == 0:
                r_core = rmid[i1]
            else:
                i0 = cand[-1]                       # nearest populated super-critical bin
                denom = np.log(prof[i0]) - np.log(prof[i1])
                if not np.isfinite(denom) or abs(denom) < 1e-300:
                    r_core = rmid[i1]
                else:
                    f = (np.log(prof[i0]) - np.log(RHOCRIT_CODE)) / denom
                    f = min(max(float(f), 0.0), 1.0)  # interpolate, never extrapolate
                    r_core = rmid[i0] + f * (rmid[i1] - rmid[i0])
        # Final guard: a non-finite or non-positive r_core must NEVER reach the masks below,
        # where it would silently yield an empty core (M_core = Phi_core = 0) that is
        # indistinguishable downstream from a real measurement. Fail loudly instead.
        if not np.isfinite(r_core) or r_core <= 0.0:
            raise ValueError(
                f"r_core is not finite/positive ({r_core}) for {path}; refusing to return an "
                "empty core, which downstream cannot distinguish from a real measurement"
            )
    else:
        r_core = float(r_core_fixed)

    in_core = r < r_core
    M_core = float((rho * dV)[in_core].sum())

    # flux through the plane z = z_peak within cylindrical radius r_core:
    # exactly ONE cell layer per column. A symmetric |dz|<=dz/2 test double-counts
    # when the plane sits exactly on a cell face (e.g. z=0 on the t=0 uniform grid
    # -- caught by the pi R^2 B0 anchor as a 2x flux error), so nudge the plane by
    # a per-cell epsilon and use a strict window.
    rcyl = np.sqrt((X - cxyz[0]) ** 2 + (Y - cxyz[1]) ** 2)
    layer = (np.abs(Z - cxyz[2] - 2e-6 * dz_half) < dz_half) & (rcyl < r_core)
    Phi_core = float((bz * dA)[layer].sum())

    return dict(time=t, r_core=r_core, M_core=M_core, Phi_core=Phi_core,
                rho_max=float(rho.max()), center=[float(c) for c in cxyz])


def lagrangian_r0(path0, M_target_code, center):
    """r0 with M(<r0) = M_target about `center` in the t=0 snapshot."""
    rho, bx, by, bz, xc, yc, zc, dxb, t = load_geom(path0)
    X, Y, Z = cell_grids(xc, yc, zc, rho.shape)
    dV = (dxb**3)[:, None, None, None] * np.ones_like(rho)
    r = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2 + (Z - center[2]) ** 2)
    order = np.argsort(r.ravel())
    mcum = np.cumsum((rho * dV).ravel()[order])
    i = np.searchsorted(mcum, M_target_code)
    i = min(i, len(mcum) - 1)
    return float(r.ravel()[order][i])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--scan-stride", type=int, default=8)
    ap.add_argument("--out", default=None)
    ap.add_argument("--milestones", default="1,100,1e5",
                    help="rho_max/rhocrit targets (deepest snapshot always added)")
    args = ap.parse_args()

    snaps = sorted(glob.glob(os.path.join(args.run_dir, "parthenon.out1.*.phdf")))
    assert snaps, f"no phdf in {args.run_dir}"
    targets = [float(x) for x in args.milestones.split(",")]

    # --- coarse scan (+ endpoint) of rho_max(t), then refine each crossing ----
    print(f"# scanning {len(snaps)} snapshots (stride {args.scan_stride})",
          flush=True)
    scan = {}
    for i in list(range(0, len(snaps), args.scan_stride)) + [len(snaps) - 1]:
        scan[i], _ = rho_max_code(snaps[i])
        print(f"  scan {i:4d} {os.path.basename(snaps[i])} rho_max={scan[i]:.4e} "
              f"(={scan[i]/RHOCRIT_CODE:.3e} rhocrit)", flush=True)

    idxs = sorted(scan)
    chosen = {}
    for tgt in targets:
        tgt_code = tgt * RHOCRIT_CODE
        # first scanned index at/above target, then linear-refine the gap
        above = [i for i in idxs if scan[i] >= tgt_code]
        if not above:
            print(f"# milestone {tgt:g} x rhocrit NEVER reached; skipping", flush=True)
            continue
        hi = above[0]
        lo = max([i for i in idxs if i < hi], default=hi)
        while hi - lo > 1:
            mid = (hi + lo) // 2
            if mid not in scan:
                scan[mid], _ = rho_max_code(snaps[mid])
                print(f"  refine {mid:4d} rho_max={scan[mid]:.4e}", flush=True)
            if scan[mid] >= tgt_code:
                hi = mid
            else:
                lo = mid
        chosen[f"{tgt:g}"] = hi
    # deepest state always included
    chosen["deepest"] = len(snaps) - 1

    # --- measure each milestone + the Lagrangian t=0 reference ---------------
    results = dict(run_dir=args.run_dir, rhocrit_code=RHOCRIT_CODE,
                   B_unit_G=B_UNIT, milestones={})
    for name, i in chosen.items():
        m = measure_snapshot(snaps[i])
        # Lagrangian initial flux for this core mass (about the box center at t=0)
        r0 = lagrangian_r0(snaps[0], m["M_core"], center=(0.0, 0.0, 0.0))
        m0 = measure_snapshot(snaps[0], r_core_fixed=r0, center=(0.0, 0.0, 0.0))
        Phi0 = m0["Phi_core"]
        M_cgs = m["M_core"] * M0
        Phi_cgs = m["Phi_core"] * B_UNIT * L0**2
        m2f = M_cgs / Phi_cgs if Phi_cgs != 0 else np.inf
        results["milestones"][name] = dict(
            snapshot=os.path.basename(snaps[i]), index=i,
            time_code=m["time"], time_kyr=m["time"] * 1.48e12 / 3.156e7 / 1e3,
            rho_max_code=m["rho_max"],
            rho_max_over_rhocrit=m["rho_max"] / RHOCRIT_CODE,
            r_core_code=m["r_core"], r_core_au=m["r_core"] * L0 / 1.496e13,
            M_core_msun=M_cgs / MSUN,
            Phi_core_code=m["Phi_core"],
            mu_core_NN78=m2f / M2F_CRIT_NN78,
            mu_core_MS76=m2f / M2F_CRIT_MS76,
            r0_lagrangian_code=r0,
            Phi0_lagrangian_code=Phi0,
            retention=(m["Phi_core"] / Phi0 if Phi0 != 0 else np.inf),
        )
        r = results["milestones"][name]
        print(f"# {name:>8s}: snap {i} t={r['time_kyr']:.2f} kyr "
              f"rho/rhocrit={r['rho_max_over_rhocrit']:.3e} "
              f"r_core={r['r_core_au']:.1f} AU M_core={r['M_core_msun']:.4f} Msun "
              f"mu(NN78)={r['mu_core_NN78']:.2f} retention={r['retention']:.4f}",
              flush=True)

    out = args.out or os.path.join(args.run_dir, "flux_retention.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"# written {out}", flush=True)


if __name__ == "__main__":
    main()
