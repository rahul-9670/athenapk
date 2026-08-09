#!/usr/bin/env python
"""WP-13b long-baseline A/B — does the DiodeBC radiation-ghost fix move the science?

THE ONLY QUESTION. old_a vs old_b measures the GPU/MPI non-determinism FLOOR at this depth.
old_a vs new measures the effect of the fix. The fix "does not move the science" ONLY if the
second sits inside the first. Without the floor leg the comparison is uninterpretable, because
4-rank GPU reductions are not order-deterministic and ~500 cycles of a collapsing turbulent flow
amplify a last-bit difference into a visible one.

WHY THE .hst AND NOT THE .phdf. The history file is full double precision and written every
dt=0.001, so it gives a dense time series of the conserved integrals -- exactly the quantities a
boundary-flux error would accumulate in. The phdf is used only to report the density reached.

WHY MATCHED TIME AND NOT MATCHED CYCLE. The three legs take different timestep sequences (that is
what non-determinism does), so cycle N is not the same physical state. Rows are matched on TIME.

WHY NEAREST-ROW AND NOT INTERPOLATION -- this cost a wrong number once, so it is written down.
The first version of this script interpolated all three legs onto a common linspace time grid.
That FABRICATED a signal: `mass` came out at 3.9e-06 against a floor of exactly 0.0, i.e. an
infinite ratio, on a quantity the code conserves to the last bit. The cause is that interpolation
mixes two neighbouring rows with weights set by where the grid point falls, and the legs' rows sit
at DIFFERENT times, so each leg gets a different blend of its own values -- a difference built by
the analysis, not by the physics. Nearest-row matching removed it: mass became exactly 0.0 in all
three legs and no other ratio moved by more than two digits. The .hst cadence (dt=0.001) is far
finer than any real drift, so nearest-row costs nothing in resolution. The worst |dt| of the match
is printed below; if it is ever comparable to the cadence, the match is not safe and the number
should not be believed.

    analyze_ab.py [run_dir]      default: this script's own directory
"""
import os, sys, glob
import numpy as np

R = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
R = os.path.abspath(R)
LEGS = ["old_a", "old_b", "new"]


def read_hst(leg):
    fs = glob.glob(os.path.join(R, leg, "*.hst"))
    if not fs:
        return None, None
    path = fs[0]
    cols = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") and "[" in line:
                import re
                cols = re.findall(r"\[\d+\]=(\S+)", line)
                break
    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d[None, :]
    return cols, d


def main():
    data = {}
    for L in LEGS:
        cols, d = read_hst(L)
        if d is None:
            print(f"MISSING hst for {L}")
            sys.exit(1)
        data[L] = (cols, d)
        print(f"  {L:6s}: {d.shape[0]} rows, t = {d[0,0]:.6f} .. {d[-1,0]:.6f}, {len(cols)} columns")

    cols = data["old_a"][0]
    # Common time window: every leg must cover it. old_a's own rows ARE the sample points -- no
    # synthetic grid, so old_a is never resampled and cannot acquire an interpolation error.
    tmax = min(data[L][1][-1, 0] for L in LEGS)
    ta = data["old_a"][1][:, 0]
    keep = ta <= tmax
    tgrid = ta[keep]
    idx = {}
    worst_dt = 0.0
    for L in LEGS:
        tl = data[L][1][:, 0]
        j = np.searchsorted(tl, tgrid)
        j = np.clip(j, 1, len(tl) - 1)
        # pick whichever of the two bracketing rows is closer
        j = np.where(np.abs(tl[j] - tgrid) < np.abs(tl[j - 1] - tgrid), j, j - 1)
        idx[L] = j
        worst_dt = max(worst_dt, float(np.abs(tl[j] - tgrid).max()))
    # Match quality must be judged against the LOCAL row spacing, not a single median cadence.
    # The .hst spacing here spans 0 to 6.5e-3 -- sparse early, then one row per cycle once the
    # collapse drives dt down -- so the median (2.0e-5, set by the dense late rows) is not a
    # meaningful denominator for an offset that occurs early. Measured against the median, the
    # deep run's largest offset (1.3e-4, at t=1.0991 where the local spacing is 1.35e-4) reads as
    # "650 % of one row" against the median and trips a NOT SAFE flag, when locally it is 0.96 of
    # one row. The first version of this guard used the median and cried wolf on a sound result,
    # which is its own kind of failure.
    # What is reported below is the worst LOCAL ratio over all rows, which need not occur at the
    # largest |dt| -- it lands on whichever row has the tightest spacing relative to its offset
    # (2.00 for the deep run). That is the honest worst case, so the flag can still fire; the
    # tight-subset cross-check underneath it is what decides whether the flag matters.
    loc = np.gradient(tgrid) if len(tgrid) > 2 else np.array([np.nan])
    mfrac = np.zeros(len(tgrid))
    for L in LEGS:
        tl = data[L][1][:, 0]
        with np.errstate(divide="ignore", invalid="ignore"):
            mfrac = np.maximum(mfrac, np.where(loc > 0, np.abs(tl[idx[L]] - tgrid) / loc, 0.0))
    worst_frac = float(np.nanmax(mfrac))
    print(f"\n  common time window: {tgrid[0]:.6f} .. {tmax:.6f}  ({len(tgrid)} matched rows)")
    print(f"  nearest-row match: worst |dt| = {worst_dt:.3e} = {worst_frac:.2f} x the LOCAL row "
          f"spacing where it occurs" + ("" if worst_frac < 1.5 else "   ** MATCH NOT SAFE"))
    # Independent falsifier: redo the comparison using ONLY rows matched to within 10 % of local
    # spacing. If the verdict depends on the loosely-matched rows, it is a matching artefact.
    tight = np.ones(len(tgrid), bool)
    for L in LEGS:
        tl = data[L][1][:, 0]
        tight &= (np.abs(tl[idx[L]] - tgrid) < 0.1 * loc)
    print(f"  tight subset (all legs within 10 % of local spacing): {tight.sum()} of {len(tgrid)}"
          f" rows -- ratios recomputed on it below as a cross-check")

    def interp(L, ci):
        return data[L][1][idx[L], ci]

    # The conserved/integral quantities a boundary flux error accumulates in.
    want = ["mass", "KE", "tot-E", "ME", "1-mom", "2-mom", "3-mom"]
    print(f"\n{'column':10s} {'scale':>12s} {'FLOOR a-b':>12s} {'SIGNAL a-new':>13s} "
          f"{'ratio':>8s} {'tight':>7s} {'sig/scale':>11s}  verdict")
    print("-" * 96)
    worst = 0.0
    worst_col = None
    for name in want:
        if name not in cols:
            continue
        ci = cols.index(name)
        a, b, n = interp("old_a", ci), interp("old_b", ci), interp("new", ci)
        scale = np.abs(a).max()
        floor = np.abs(a - b).max()
        sig = np.abs(a - n).max()
        ratio = (sig / floor) if floor > 0 else (np.inf if sig > 0 else 1.0)
        frac = sig / scale if scale > 0 else 0.0
        # Same ratio on the tightly-matched rows only. If the two disagree, the verdict is being
        # carried by rows where the time match is loose, i.e. it is an artefact of the matching.
        fl_t = np.abs(a[tight] - b[tight]).max() if tight.any() else 0.0
        sg_t = np.abs(a[tight] - n[tight]).max() if tight.any() else 0.0
        rat_t = (sg_t / fl_t) if fl_t > 0 else (np.inf if sg_t > 0 else 1.0)
        # A column is only flagged when the fix moves it well beyond the floor AND the move is
        # non-negligible against the column's own size. Ratio alone flags round-off.
        flag = (ratio > 10.0) and (frac > 1.0e-6)
        if flag and ratio > worst:
            worst, worst_col = ratio, name
        print(f"{name:10s} {scale:12.5e} {floor:12.5e} {sig:13.5e} {ratio:8.2f} {rat_t:7.2f} "
              f"{frac:11.2e}  {'** ABOVE FLOOR' if flag else 'within floor'}")

    print("\n" + "=" * 88)
    if worst_col is None:
        print("  RESULT: no integral quantity moves more than 10x the non-determinism floor")
        print("          AND more than 1e-6 of its own magnitude. Over this baseline the")
        print("          DiodeBC fix is indistinguishable from run-to-run non-determinism.")
    else:
        print(f"  RESULT: {worst_col} moves {worst:.1f}x the non-determinism floor.")
        print("          The fix CHANGES the science at this depth; the ensemble numbers")
        print("          produced with 84a6d248 must be re-examined, not just re-labelled.")
    print("=" * 88)

    # Depth actually reached -- the honest statement of how far this bounds anything.
    print("\n  depth reached (rho_max from the last readable snapshot):")
    try:
        import h5py
        for L in LEGS:
            fs = sorted(glob.glob(os.path.join(R, L, "parthenon.out1.*.phdf")),
                        key=os.path.getmtime)
            for f in reversed(fs):
                try:
                    with h5py.File(f, "r") as h:
                        r = float(np.array(h["prim"][:, 0, ...]).max())
                        t = float(h["Info"].attrs["Time"])
                        c = int(h["Info"].attrs["NCycle"])
                    print(f"    {L:6s} {os.path.basename(f):26s} cycle={c:5d} t={t:.5f} "
                          f"rho_max={r*5.467e-19:.4e} g/cm3 = {r*5.467e-19/1e-13:6.2f}x rhocrit")
                    break
                except Exception:
                    continue
    except ImportError:
        print("    (h5py unavailable)")


if __name__ == "__main__":
    main()
