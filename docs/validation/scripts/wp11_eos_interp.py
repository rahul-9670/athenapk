#!/usr/bin/env python
"""WP-11 — tabulated-EOS interpolation error.

The production run reads a tabulated protostellar EOS and interpolates it BILINEARLY in
(log10 rho, log10 esp).  WP-11 asks how much error that interpolation introduces, i.e.
whether the table is fine enough that the EOS is not itself a source of error in the
collapse.

Method: `src/eos/eos_table_hires.bin` exists on disk and is a refined version of the same
table.  Use it as truth.  For every hi-res grid node that lies inside the production table's
domain, evaluate the PRODUCTION table's bilinear interpolant there and compare against the
hi-res value at that same physical point.  Nodes of the coarse table are excluded from the
error statistics: the interpolant is exact there by construction, so including them would
dilute the error with structural zeros.

Format (from src/eos/eos_table.hpp:119-148):
    int64 nr, ne, nT
    float64 lr0, dlr, le0, dle, lT0, dlT      (log10 axis origins + spacings)
    float64 P[nr][ne], cs2[nr][ne], logT[nr][ne], espT[nr][nT]
Interpolation must mirror `bilin` in the same header, INCLUDING its edge clamping.
"""
import sys
import numpy as np

REPO = "/beegfs/u/bbg6470/athenapk"


def load(path):
    with open(path, "rb") as f:
        nr, ne, nT = np.frombuffer(f.read(24), dtype=np.int64)
        lr0, dlr, le0, dle, lT0, dlT = np.frombuffer(f.read(48), dtype=np.float64)
        def rd(n1, n2):
            return np.frombuffer(
                f.read(n1 * n2 * 8), dtype=np.float64).reshape(n1, n2).copy()
        P, cs2, logT = rd(nr, ne), rd(nr, ne), rd(nr, ne)
        espT = rd(nr, nT)
    return dict(nr=int(nr), ne=int(ne), nT=int(nT), lr0=lr0, dlr=dlr, le0=le0, dle=dle,
                lT0=lT0, dlT=dlT, P=P, cs2=cs2, logT=logT, espT=espT)


def bilin(tab, key, x, y, n1, n2, o1, d1, o2, d2):
    """Vectorised copy of EosTable::bilin, clamping at the edges exactly as the C++ does."""
    a = np.asarray(tab[key])
    fx = (x - o1) / d1
    fy = (y - o2) / d2
    i = np.clip(np.floor(fx).astype(int), 0, n1 - 2)
    j = np.clip(np.floor(fy).astype(int), 0, n2 - 2)
    tx = np.clip(fx - i, 0.0, 1.0)
    ty = np.clip(fy - j, 0.0, 1.0)
    return ((1 - tx) * (1 - ty) * a[i, j] + tx * (1 - ty) * a[i + 1, j]
            + (1 - tx) * ty * a[i, j + 1] + tx * ty * a[i + 1, j + 1])


def report(name, coarse, fine, key, axes):
    """axes = ('r','e') or ('r','T') -- which pair of log axes `key` is tabulated on."""
    a1, a2 = axes
    n1f = fine["nr"]
    n2f = fine["ne"] if a2 == "e" else fine["nT"]
    o1f, d1f = fine["lr0"], fine["dlr"]
    o2f, d2f = (fine["le0"], fine["dle"]) if a2 == "e" else (fine["lT0"], fine["dlT"])
    n1c = coarse["nr"]
    n2c = coarse["ne"] if a2 == "e" else coarse["nT"]
    o1c, d1c = coarse["lr0"], coarse["dlr"]
    o2c, d2c = (coarse["le0"], coarse["dle"]) if a2 == "e" else (coarse["lT0"], coarse["dlT"])

    x = o1f + d1f * np.arange(n1f)
    y = o2f + d2f * np.arange(n2f)
    X, Y = np.meshgrid(x, y, indexing="ij")

    # Only points strictly inside the coarse domain -- outside it `bilin` clamps, and a
    # clamping error is an extrapolation artefact, not an interpolation error.
    inside = ((X >= o1c) & (X <= o1c + d1c * (n1c - 1))
              & (Y >= o2c) & (Y <= o2c + d2c * (n2c - 1)))
    # Exclude coarse-grid nodes: exact by construction.
    on_node = (np.isclose((X - o1c) / d1c, np.round((X - o1c) / d1c), atol=1e-9)
               & np.isclose((Y - o2c) / d2c, np.round((Y - o2c) / d2c), atol=1e-9))
    m = inside & ~on_node
    if m.sum() == 0:
        print(f"  {name:<10} no off-node interior points -- grids coincide; skipped")
        return None

    approx = bilin(coarse, key, X[m], Y[m], n1c, n2c, o1c, d1c, o2c, d2c)
    truth = np.asarray(fine[key])[m]

    denom = np.maximum(np.abs(truth), 1e-300)
    rel = np.abs(approx - truth) / denom
    finite = np.isfinite(rel)
    rel = rel[finite]
    print(f"  {name:<10} N={rel.size:>9d}   max={rel.max():.3e}   "
          f"p99={np.percentile(rel,99):.3e}   median={np.median(rel):.3e}   "
          f"rms={np.sqrt((rel**2).mean()):.3e}")
    return rel


def main():
    cp = f"{REPO}/src/eos/eos_table.bin"
    fp = f"{REPO}/src/eos/eos_table_hires.bin"
    c, f = load(cp), load(fp)

    for tag, t, p in (("production", c, cp), ("hi-res", f, fp)):
        span_r = t["lr0"] + t["dlr"] * (t["nr"] - 1)
        span_e = t["le0"] + t["dle"] * (t["ne"] - 1)
        span_T = t["lT0"] + t["dlT"] * (t["nT"] - 1)
        print(f"{tag:<11} nr={t['nr']:<5} ne={t['ne']:<5} nT={t['nT']:<5}  "
              f"log10 rho [{t['lr0']:.3f},{span_r:.3f}] d={t['dlr']:.5f}  "
              f"log10 esp [{t['le0']:.3f},{span_e:.3f}] d={t['dle']:.5f}  "
              f"log10 T [{t['lT0']:.3f},{span_T:.3f}] d={t['dlT']:.5f}")
    print()

    ratio_r = f["dlr"] / c["dlr"] if c["dlr"] else float("nan")
    ratio_e = f["dle"] / c["dle"] if c["dle"] else float("nan")
    print(f"hi-res spacing / production spacing:  rho {ratio_r:.4f}   esp {ratio_e:.4f}")
    if ratio_r >= 1.0 or ratio_e >= 1.0:
        print("*** hi-res table is NOT finer on every axis -- it cannot serve as truth. "
              "WP-11 needs a separately generated reference. ***")
        return 2
    print()

    print("Relative interpolation error of the PRODUCTION table, hi-res table as truth")
    print("(coarse-grid nodes excluded -- exact there by construction):")
    out = {}
    for name, key, axes in (("P", "P", ("r", "e")), ("cs2", "cs2", ("r", "e")),
                            ("logT", "logT", ("r", "e")), ("espT", "espT", ("r", "T"))):
        out[name] = report(name, c, f, key, axes)

    print()
    print("ACCEPTANCE (VALIDATION_PLAN WP-11): the EOS must not be a leading error term.")
    print("Compare against WP-18's realization scatter, sigma ~ 16% on MEtor/MEpol.")
    worst = max((r.max() for r in out.values() if r is not None), default=float("nan"))
    print(f"worst max relative interpolation error over all four tables: {worst:.3e}"
          f"  ({100*worst:.4f} %)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
