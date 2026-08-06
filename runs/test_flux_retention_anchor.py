#!/usr/bin/env python
"""Executable anchor test for flux_retention.measure_snapshot.

WHY THIS FILE EXISTS. `flux_retention.py` is the toolkit that produces the flagship's headline
observable (mu_core = M_core / Phi_core), and it had NO executable test -- only an inline comment
recording that a "pi R^2 B0 anchor" had once caught a 2x flux error. A validated toolkit with no
runnable check is how a bug survives: on 2026-08-06 six of seven finished ensemble members came
back with Phi_core exactly zero and nothing flagged it, because a zero core is indistinguishable
downstream from a legitimate measurement.

THE ANCHOR. At t = 0 the collapse_be IC lays down a UNIFORM vertical field B0z on a uniform grid.
The magnetic flux through a disc of radius R centred on the axis is then analytic:

        Phi(R) = pi * R^2 * B0z

so `measure_snapshot` with a fixed core radius must reproduce that to within the discretisation
error of tiling a circle with square cells (a few percent at these radii, shrinking with R/dx).
This is a genuine external check: it depends on no other part of the pipeline and on no
previously-recorded number.

It also pins the failure mode the comment describes -- if the z-plane selection picks up two cell
layers instead of one, the answer comes out 2x high, which this test would catch immediately.

USAGE
    test_flux_retention_anchor.py [snapshot.phdf]
defaults to the t=0 snapshot of ensemble member point000.
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/beegfs/u/bbg6470/athenapk/runs")
import flux_retention as fr

DEFAULT = ("/beegfs/u/bbg6470/athenapk/runs/ensemble/design01/point000/"
           "parthenon.out1.00000.phdf")
# B0z for the ensemble decks is a sampled IC parameter; read it from the member's own deck rather
# than hardcoding, so the test cannot silently drift from the run it is checking.
DECK = ("/beegfs/u/bbg6470/athenapk/runs/ensemble/design01/point000/fhc_ens.in")


def b0z_from_deck(path):
    for line in open(path):
        s = line.split("#", 1)[0].strip()
        if s.startswith("B0z"):
            return float(s.split("=", 1)[1])
    raise SystemExit(f"no B0z in {path}")


def main():
    snap = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    if not os.path.exists(snap):
        raise SystemExit(f"missing snapshot {snap}")
    b0z = b0z_from_deck(DECK)
    print(f"snapshot : {os.path.basename(snap)}")
    print(f"B0z      : {b0z:.6g} (code, Heaviside-Lorentz, from the deck)")

    rho, bx, by, bz, xc, yc, zc, dxb, t = fr.load_geom(snap)
    print(f"time     : {t:.6g}  (anchor is only valid at t=0)")
    if abs(t) > 1e-12:
        print("   ** WARNING: not a t=0 snapshot; the uniform-field assumption does not hold **")

    # Sanity: the field really is uniform and vertical at t=0.
    print(f"bz       : min={bz.min():.6g} max={bz.max():.6g}  "
          f"(uniform? {np.allclose(bz, b0z, rtol=1e-6)})")
    print(f"bx,by    : max|bx|={np.abs(bx).max():.3e}  max|by|={np.abs(by).max():.3e}")

    print(f"\n{'R (code)':>10} {'Phi measured':>15} {'pi R^2 B0z':>13} {'rel err':>10}")
    worst = 0.0
    for R in (0.5, 1.0, 2.0, 3.0):
        m = fr.measure_snapshot(snap, r_core_fixed=R)
        exact = np.pi * R * R * b0z
        rel = (m["Phi_core"] - exact) / exact
        worst = max(worst, abs(rel))
        print(f"{R:10.3f} {m['Phi_core']:15.6e} {exact:13.6e} {100*rel:9.2f}%")

    # Tiling a circle with square cells converges slowly; 5% is generous at R/dx ~ 8 and tight
    # enough to catch the 2x (100%) plane-selection error the module comment describes.
    tol = 0.05
    print(f"\nworst |rel err| = {100*worst:.2f}%   tolerance {100*tol:.0f}%")
    if worst <= tol:
        print("ANCHOR PASS — flux integration reproduces pi R^2 B0z")
        return 0
    print("ANCHOR FAIL — flux integration does NOT reproduce pi R^2 B0z.")
    print("  ~100% high => the z-plane window is selecting two cell layers, not one.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
