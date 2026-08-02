#!/usr/bin/env python
"""WP-11 (part 2) — interpolation error of the table PRODUCTION ACTUALLY USES.

Part 1 (`wp11_eos_interp.py`) measured `eos_table.bin` (180x220) against
`eos_table_hires.bin` (400x1000).  That is NOT the production answer: both production decks
set

    eos_table_file = /beegfs/u/bbg6470/athenapk/src/eos/eos_table_hires.bin

(`runs/prod_flagship_test/fhc_flagship.in:59`, `runs/root_ladder/fhc_rootladder.in:77`).
The 180x220 table is only the compiled-in default (`src/hydro/hydro.cpp:790`).  Part 1 is
therefore a bound on the DEFAULT table, and is the number that would apply if a deck ever
omitted the key.

There is no finer table on disk to use as truth for the hi-res table, so estimate its error by
SELF-CONVERGENCE: decimate the hi-res table by s = 2,3,4,5, measure each decimation's error
against the hi-res values at off-node points, fit the observed order p in err ~ h^p, and
Richardson-extrapolate back to s = 1.

Bilinear interpolation is O(h^2) for smooth data.  A measured p well below 2 is itself the
finding: it means the table is NOT smooth on its grid scale, which for a protostellar EOS is
expected exactly at the H2-dissociation and H-ionization kinks -- the regions that set first-
and second-core formation.
"""
import sys
import numpy as np

REPO = "/beegfs/u/bbg6470/athenapk"
HIRES = f"{REPO}/src/eos/eos_table_hires.bin"


def load(path):
    with open(path, "rb") as f:
        nr, ne, nT = np.frombuffer(f.read(24), dtype=np.int64)
        lr0, dlr, le0, dle, lT0, dlT = np.frombuffer(f.read(48), dtype=np.float64)
        def rd(n1, n2):
            return np.frombuffer(f.read(n1 * n2 * 8),
                                 dtype=np.float64).reshape(n1, n2).copy()
        P, cs2, logT = rd(nr, ne), rd(nr, ne), rd(nr, ne)
        espT = rd(nr, nT)
    return dict(nr=int(nr), ne=int(ne), nT=int(nT), lr0=lr0, dlr=dlr, le0=le0, dle=dle,
                lT0=lT0, dlT=dlT, P=P, cs2=cs2, logT=logT, espT=espT)


def bilin_on(a, fx, fy):
    """Bilinear sample of array `a` at fractional indices (fx, fy), clamped like the C++."""
    n1, n2 = a.shape
    i = np.clip(np.floor(fx).astype(int), 0, n1 - 2)
    j = np.clip(np.floor(fy).astype(int), 0, n2 - 2)
    tx = np.clip(fx - i, 0.0, 1.0)
    ty = np.clip(fy - j, 0.0, 1.0)
    return ((1 - tx) * (1 - ty) * a[i, j] + tx * (1 - ty) * a[i + 1, j]
            + (1 - tx) * ty * a[i, j + 1] + tx * ty * a[i + 1, j + 1])


def err_at_stride(a, s):
    """Decimate `a` by stride s, interpolate back to every fine node, report error there.

    Only off-node fine points are scored (exact on the decimated nodes by construction)."""
    coarse = a[::s, ::s]
    n1, n2 = a.shape
    ii = np.arange(n1)
    jj = np.arange(n2)
    # a fine node i maps to coarse fractional index i/s
    keep_i = ii[ii <= (coarse.shape[0] - 1) * s]
    keep_j = jj[jj <= (coarse.shape[1] - 1) * s]
    I, J = np.meshgrid(keep_i, keep_j, indexing="ij")
    off = (I % s != 0) | (J % s != 0)
    I, J = I[off], J[off]
    approx = bilin_on(coarse, I / s, J / s)
    truth = a[I, J]
    rel = np.abs(approx - truth) / np.maximum(np.abs(truth), 1e-300)
    rel = rel[np.isfinite(rel)]
    return rel


def main():
    t = load(HIRES)
    print(f"production table: {HIRES}")
    print(f"  nr={t['nr']} ne={t['ne']} nT={t['nT']}   "
          f"dlog10rho={t['dlr']:.5f} dlog10esp={t['dle']:.5f} dlog10T={t['dlT']:.5f}")
    print()

    strides = [2, 3, 4, 5]
    for name in ("P", "cs2", "logT", "espT"):
        a = t[name]
        print(f"--- {name} ---")
        maxs, rmss, p99s = [], [], []
        for s in strides:
            rel = err_at_stride(a, s)
            maxs.append(rel.max()); rmss.append(np.sqrt((rel ** 2).mean()))
            p99s.append(np.percentile(rel, 99))
            print(f"   stride {s}  h/h0={s:<2d}  max={rel.max():.3e}  "
                  f"p99={np.percentile(rel,99):.3e}  rms={np.sqrt((rel**2).mean()):.3e}")
        # observed order from the rms trend (least squares on log-log)
        ls = np.log(np.array(strides, dtype=float))
        for lbl, arr in (("rms", rmss), ("p99", p99s), ("max", maxs)):
            p = np.polyfit(ls, np.log(np.array(arr)), 1)[0]
            est = arr[0] / (2.0 ** p)          # Richardson from stride 2 back to stride 1
            print(f"   observed order p({lbl}) = {p:.2f}   "
                  f"=> production-table {lbl} error ~ {est:.3e}  ({100*est:.4f} %)")
        print()

    print("Bilinear interpolation of smooth data is O(h^2), i.e. p = 2.")
    print("p noticeably below 2 means the table is not smooth on its own grid scale.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
