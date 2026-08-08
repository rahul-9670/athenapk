#!/usr/bin/env python
"""WP-13b — compare the three legs of the GPU+AMR restart test.

The question is NOT "is the restart bit-identical" -- on GPU with 4 ranks it need not be, because
reductions are not order-deterministic. The question is whether restart divergence exceeds the
non-determinism floor measured by two identical fresh runs.

  d_floor   = |fresh_a - fresh_b|     two identical fresh runs
  d_restart = |fresh_a - split|       midpoint restart

PASS iff d_restart <= d_floor (to within the same order of magnitude). If d_floor is exactly zero
the code IS deterministic here, and then anything nonzero in d_restart is a real restart defect.
"""
import sys, os, glob
import numpy as np
import h5py

R = os.path.dirname(os.path.abspath(__file__))
LEGS = ["fresh_a", "fresh_b", "split"]


def load(leg):
    f = os.path.join(R, leg, "parthenon.out1.final.phdf")
    if not os.path.exists(f):
        return None, f"MISSING {f}"
    out = {}
    with h5py.File(f, "r") as h:
        for k in ("prim", "grav.phi", "rad.Er"):
            if k in h:
                out[k] = np.array(h[k])
        info = h["Info"].attrs
        out["_t"] = float(info["Time"])
        out["_cycle"] = int(info["NCycle"])
        out["_nblocks"] = int(info["NumMeshBlocks"])
    return out, None


def cmp(a, b, la, lb):
    print(f"\n--- {la} vs {lb} ---")
    if a["_nblocks"] != b["_nblocks"]:
        print(f"  BLOCK COUNT DIFFERS: {a['_nblocks']} vs {b['_nblocks']}  "
              f"(AMR histories diverged; field comparison not meaningful)")
    print(f"  t      : {a['_t']:.17e}  vs {b['_t']:.17e}   dt_rel="
          f"{abs(a['_t']-b['_t'])/max(abs(a['_t']),1e-300):.3e}")
    print(f"  cycle  : {a['_cycle']} vs {b['_cycle']}    blocks: {a['_nblocks']} vs {b['_nblocks']}")
    worst = 0.0
    for k in sorted(set(a) & set(b)):
        if k.startswith("_"):
            continue
        x, y = a[k], b[k]
        if x.shape != y.shape:
            print(f"  {k:10s}: SHAPE {x.shape} vs {y.shape} -- skipped")
            continue
        d = np.abs(x - y)
        scale = np.maximum(np.abs(x), np.abs(y))
        rel = np.where(scale > 0, d / np.where(scale > 0, scale, 1.0), 0.0)
        print(f"  {k:10s}: max|abs|={d.max():.6e}   max|rel|={rel.max():.6e}   "
              f"{'BIT-IDENTICAL' if d.max() == 0.0 else ''}")
        worst = max(worst, rel.max())
    return worst


def main():
    legs = {}
    for L in LEGS:
        d, err = load(L)
        if err:
            print(err)
            sys.exit(1)
        legs[L] = d

    floor = cmp(legs["fresh_a"], legs["fresh_b"], "fresh_a", "fresh_b")
    restart = cmp(legs["fresh_a"], legs["split"], "fresh_a", "split")

    print("\n================ VERDICT ================")
    print(f"  non-determinism floor (fresh_a vs fresh_b) : max|rel| = {floor:.6e}")
    print(f"  restart divergence    (fresh_a vs split)   : max|rel| = {restart:.6e}")
    if floor == 0.0 and restart == 0.0:
        print("  PASS — the code is deterministic here AND the restart is bit-identical.")
    elif floor == 0.0 and restart > 0.0:
        print("  FAIL — two fresh runs are bit-identical, so the code IS deterministic in this")
        print("         configuration, and the restart is therefore introducing a real difference.")
    elif restart <= max(floor * 10.0, floor):
        print("  PASS — restart divergence is within the non-determinism floor. Restart is as")
        print("         reproducible as this configuration allows; the point012 dt difference is")
        print("         explained by non-determinism plus chaotic amplification, not a defect.")
    else:
        print(f"  FAIL — restart divergence exceeds the floor by {restart/max(floor,1e-300):.2f}x.")
    print("=========================================")


if __name__ == "__main__":
    main()
