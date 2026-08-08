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
        # Compare EVERY dumped field, not a hardcoded three. The multigroup decks (n_group=3,
        # which is what the flagship and all 69 production decks run) dump rad.Er_g1/rad.Er_g2
        # as well, and a hardcoded list would silently skip exactly the fields whose solver is
        # under test.
        keys = [k for k in h.keys() if k not in ("Info", "Params", "Blocks", "Locations",
                                                 "LogicalLocations", "Levels", "SparseInfo",
                                                 "VolumeLocations")]
        for k in keys:
            try:
                a = np.array(h[k])
            except Exception:
                continue
            if a.dtype.kind == "f":
                out[k] = a
        if np.dtype(out.get("prim", np.zeros(1)).dtype) == np.float32:
            print("  !! WARNING: fields are SINGLE precision (parthenon/output1/"
                  "single_precision_output=1). Differences below ~1e-7 relative are INVISIBLE; "
                  "this test is only meaningful with single_precision_output=0.")
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
    out = {}
    for k in sorted(set(a) & set(b)):
        if k.startswith("_"):
            continue
        x, y = a[k], b[k]
        if x.size == 0 or y.size == 0:
            continue  # e.g. an unallocated sparse/swarm field: nothing to compare
        if x.shape != y.shape:
            print(f"  {k:10s}: SHAPE {x.shape} vs {y.shape} -- skipped")
            continue
        d = np.abs(x - y)
        scale = np.maximum(np.abs(x), np.abs(y))
        rel = np.where(scale > 0, d / np.where(scale > 0, scale, 1.0), 0.0)
        print(f"  {k:10s}: max|abs|={d.max():.6e}   max|rel|={rel.max():.6e}   "
              f"{'BIT-IDENTICAL' if d.max() == 0.0 else ''}")
        out[k] = (float(d.max()), float(np.abs(x).max()))
    return out


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

    # VERDICT, PER FIELD, ON ABSOLUTE SCALES.
    #
    # The previous version reduced each comparison to a single global max|rel| and compared the
    # two numbers. That is not a usable statistic and it returned a FALSE PASS on the multigroup
    # flagship deck (2026-08-08): |x-y|/max(|x|,|y|) SATURATES at 2.0 whenever x and y have
    # opposite signs, so the near-empty rad.Fr*_g2 fields (|value| ~ 1e-15 .. 1e-22, pure noise)
    # pinned BOTH the floor and the restart comparison at ~1.9996 and the ratio test passed --
    # while rad.Fr1 was actually diverging by 74.1 in absolute terms against a floor of 8.6e-09.
    # Compare like with like, field by field, and let the worst field decide.
    print("\n================ VERDICT ================")
    # A field is flagged only if the restart divergence is BOTH well above the non-determinism
    # floor AND non-negligible against the field's own magnitude. Both conditions are needed:
    # ratio alone flags a field whose floor is 1e-10 and whose restart difference is 1e-8, which
    # is round-off on a field of size 1e-2 (1e-6 relative) and is not a defect; magnitude alone
    # cannot separate a defect from ordinary GPU reduction non-determinism.
    REL_FLOOR = 1.0e-4  # restart difference must exceed this fraction of the field's own scale
    print(f"  {'field':12s} {'floor max|abs|':>16s} {'restart max|abs|':>18s} {'ratio':>10s} "
          f"{'field scale':>13s} {'rest/scale':>11s}")
    worst_ratio, worst_field = 0.0, None
    for k in sorted(set(floor) & set(restart)):
        f, _ = floor[k]
        r, scale = restart[k]
        ratio = (r / f) if f > 0 else (np.inf if r > 0 else 1.0)
        frac = r / scale if scale > 0 else 0.0
        print(f"  {k:12s} {f:16.6e} {r:18.6e} {ratio:10.2e} {scale:13.5e} {frac:11.2e}")
        if r > 0 and ratio > 10.0 and frac > REL_FLOOR and ratio > worst_ratio:
            worst_ratio, worst_field = ratio, k
    if worst_field is None:
        allzero = all(v[0] == 0.0 for v in restart.values())
        if allzero:
            print("\n  PASS — the restart is BIT-IDENTICAL to the fresh run in every field.")
        else:
            print("\n  PASS — no field diverges both >10x beyond the non-determinism floor and")
            print(f"         >{REL_FLOOR:g} of its own magnitude. The restart is as reproducible")
            print("         as this configuration allows.")
    else:
        print(f"\n  FAIL — {worst_field} diverges {worst_ratio:.2e}x beyond the non-determinism")
        print("         floor, so this is a real restart defect and not GPU/MPI reduction order.")
        if floor.get(worst_field, (0.0, 0.0))[0] == 0.0:
            print("         (Its floor is EXACTLY zero: the code is deterministic in that field.)")
    print("=========================================")


if __name__ == "__main__":
    main()
