#!/usr/bin/env python
"""WP-8 — does the DENSITY SPLIT restore convergence for the DISSIPATION budget?

THE QUESTION, and why this is the one that matters. `mag_diag.hpp` diagnosed the global
`mag-dissO`/`mag-dissA` as ill-conditioned: 90 % of int |J|^2 dV comes from a volume fraction of
~1e-6, so these "volume integrals" are point samples of the innermost core, taken in exactly the
region whose resolution changes between ladder rungs. Its prescribed remedy was to split the
integrals on DENSITY -- regions defined by physics rather than by the grid -- so the core and
envelope budgets converge (or fail to) independently and visibly.

That remedy was never measured. Two OTHER splits were tried on `Jsq` and both failed:
  * density split on Jsq (2026-08-06): the low-density bin holds 97-99 % at f_eff ~ 1e-7 --
    the pathology under a new name (-51.6 %/-39.9 % vs the global -52.3 %/-39.3 %).
  * current-sheet split on Jsq (2026-08-08): no bin converges; smooth@0.1 is not even monotone.
Neither touched dissO/dissA, because eta is in no output file and cannot be reconstructed offline
(dust_coupling=true makes it depend on the evolved grain state, which nothing writes). The fix was
to RESTART each rung with the split columns enabled -- eta is {Derived}, so a restart regenerates
it -- using a binary that is the ladder's own physics plus the backported read-only diagnostic.

WHAT COUNTS AS SUCCESS, stated before looking. A bin is a useful convergence metric only if
  (i) it converges: |change| shrinks with refinement and is monotone, AND
  (ii) it is not a point sample: f_eff = (int q dV)^2 / (V_box * int q^2 dV) is not ~1e-7.
A bin holding 99 % of the integral at f_eff ~ 1e-7 is the original pathology relabelled. Both are
reported for every bin so the second cannot be quietly skipped.

EPOCH MATCHING. Every rung is read at rho_max = 1e-12 g/cm^3, the ladder's matched epoch. nj4's
restart sits ON it, so its first history row is the epoch. nj8/nj16 pass THROUGH it, so their
epoch time is interpolated from the (time, rho_max) pairs of their own snapshots and the nearest
history row is taken. The offset actually used is printed -- if it is ever large compared with the
local row spacing, the number should not be believed.

    wp8_diss_convergence.py [rho_target_cgs]        default 1e-12
"""
import sys, os, glob, re
import numpy as np

RHO0 = 5.467e-19          # code density unit [g/cm^3]
VBOX = 16.0 ** 3          # ladder domain [-8,8]^3 in code units; AMR does not change it
W8 = "/beegfs/u/bbg6470/athenapk/runs/wp8_dissplit"
# All three rungs now come from one identical code path (submit_rung.sh with NJ=4/8/16) on the
# round-3 binary, so every leg carries the sheet-split dissipation and per-bin sq columns. The
# earlier nj4 point lived in the v2_on gate run, which predates those columns.
LEGS = [("nj4", f"{W8}/nj4"), ("nj8", f"{W8}/nj8"), ("nj16", f"{W8}/nj16")]


def read_hst(d):
    fs = glob.glob(os.path.join(d, "*.hst"))
    if not fs:
        return None, None
    cols = []
    with open(fs[0]) as fh:
        for line in fh:
            if line.startswith("#") and "[" in line:
                cols = re.findall(r"\[\d+\]=(\S+)", line)
                break
    a = np.loadtxt(fs[0], comments="#")
    return cols, (a[None, :] if a.ndim == 1 else a)


def ladder_epoch_time(leg, rho_t):
    """Matched-epoch TIME from the ladder's own epoch_scan.txt — the authoritative source.

    Preferred over this run's snapshots for two reasons. (1) It is exactly what
    wp8_split_convergence.py and wp8_sheet_convergence.py used, so these numbers are directly
    comparable to the published density- and sheet-split results instead of being anchored to a
    slightly different epoch. (2) The restarted legs write few snapshots, so their own rho(t)
    sampling is coarse: nj4's gate run wrote only `final`, which would have pinned it at
    1.16e-12 (0.07 dex) when its restart sits at 9.80e-13 (0.009 dex). The restarts follow the
    ladder's own trajectory, so the ladder's epoch time transfers directly.
    """
    p = "/beegfs/u/bbg6470/athenapk/runs/convergence_ladder/epoch_scan.txt"
    best = None
    try:
        for line in open(p):
            f = line.split()
            if len(f) < 6 or f[0] != leg:
                continue
            d = abs(np.log10(float(f[5]) / rho_t))
            if best is None or d < best[0]:
                best = (d, float(f[2]), f[1])
    except OSError:
        return None, None, None
    return (None, None, None) if best is None else (best[1], best[0], best[2])


def epoch_time(d, rho_t):
    """(time of rho_max = rho_t, how it was obtained). None if the leg never reached it."""
    pts = []
    for f in sorted(glob.glob(os.path.join(d, "parthenon.out1.*.phdf"))):
        try:
            import h5py
            with h5py.File(f, "r") as h:
                pts.append((float(h["Info"].attrs["Time"]),
                            float(np.array(h["prim"][:, 0, ...]).max()) * RHO0))
        except Exception:
            continue
    if not pts:
        return None, "no snapshots"
    pts.sort()
    t = np.array([p[0] for p in pts]); r = np.array([p[1] for p in pts])
    if r.max() < rho_t:
        return None, f"never reached target (max {r.max():.3e})"
    if r.min() > rho_t:                       # already past at the first snapshot
        return t[0], f"restart already at/past target ({r[0]:.3e})"
    # interpolate in log rho, which is near-linear in t during collapse
    return float(np.interp(np.log10(rho_t), np.log10(r), t)), "interpolated between snapshots"


def main():
    rho_t = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-12
    print(f"WP-8 DISSIPATION split on the njeans ladder — matched epoch "
          f"rho_max = {rho_t:.3e} g/cm^3\n")

    res = {}
    for leg, d in LEGS:
        cols, a = read_hst(d)
        if a is None:
            print(f"  {leg}: NO .hst in {d} — leg missing, cannot complete the ladder")
            return 1
        # Prefer the ladder's own matched epoch, so these numbers sit on exactly the epoch the
        # published density- and sheet-split results used. Fall back to this run's snapshots only
        # if epoch_scan.txt cannot supply it.
        te, dex, snap = ladder_epoch_time(leg, rho_t)
        how = f"ladder epoch_scan {snap}, {dex:.3f} dex from target" if te is not None else None
        if te is None:
            te, how = epoch_time(d, rho_t)
        if te is None:
            print(f"  {leg}: {how} — cannot place the matched epoch")
            return 1
        # The epoch must lie inside the restarted leg's own time span, or we would be reading a
        # row the leg never reached and silently reporting the nearest endpoint instead.
        if not (a[0, 0] - 1e-9 <= te <= a[-1, 0] + 1e-9):
            print(f"  {leg}: matched epoch t={te:.6f} is OUTSIDE this leg's span "
                  f"[{a[0,0]:.6f}, {a[-1,0]:.6f}] — the leg does not cover the epoch")
            return 1
        # EXCLUDE THE RESTART'S FIRST HISTORY ROW. The eta cache is filled during the first stage
        # of a step, so on the row written at restart time it is still zero and every eta-weighted
        # column (mag-dissO, mag-dissA and all four split bins) reads EXACTLY 0. Jsq is unaffected
        # because it needs no eta. This bites hardest on nj4, whose restart sits ON the matched
        # epoch: the nearest row to the epoch IS row 0, so an unguarded argmin returns 0.0 for the
        # whole dissipation budget and it looks like a physical result, not a startup artefact.
        # Verified to be row 0 only, in both legs (nj4 and nj8): the very next row is fully
        # populated (nj4 3.50e-03 at t=1.093700; nj8 1.57e-05 at t=1.081480).
        valid = np.ones(len(a), bool)
        oi = cols.index("mag-dissO") if "mag-dissO" in cols else None
        if oi is not None:
            valid = a[:, oi] != 0.0
            if not valid.any():
                print(f"  {leg}: mag-dissO is zero in EVERY row — the eta cache never filled")
                return 1
        i = int(np.where(valid)[0][np.argmin(np.abs(a[valid, 0] - te))])
        if i != int(np.argmin(np.abs(a[:, 0] - te))):
            print(f"        (skipped the restart row, where the eta cache is still empty)")
        spacing = np.median(np.diff(a[:, 0])) if len(a) > 2 else float("nan")
        off = abs(a[i, 0] - te)
        flag = "" if not (spacing > 0) or off < 2 * spacing else "  ** OFFSET > 2 rows"
        print(f"  {leg:5s} rows={len(a):4d}  epoch t={te:.6f} ({how})")
        print(f"        using row {i} at t={a[i,0]:.6f}, offset {off:.2e} "
              f"(median row spacing {spacing:.2e}){flag}")
        res[leg] = dict(zip(cols, a[i]))
    print()

    def g(leg, name):
        return res[leg].get(name, float("nan"))

    bins = [
        ("dissO  global", "mag-dissO",    None),
        ("dissO  core",   "mag-dissO-hi", "mag-dissO"),
        ("dissO  envel",  "mag-dissO-lo", "mag-dissO"),
        ("dissA  global", "mag-dissA",    None),
        ("dissA  core",   "mag-dissA-hi", "mag-dissA"),
        ("dissA  envel",  "mag-dissA-lo", "mag-dissA"),
        ("dissO  sheet",  "mag-dissO-sheet",  "mag-dissO"),
        ("dissO  smooth", "mag-dissO-smooth", "mag-dissO"),
        ("dissA  sheet",  "mag-dissA-sheet",  "mag-dissA"),
        ("dissA  smooth", "mag-dissA-smooth", "mag-dissA"),
        ("Jsq    global", "mag-Jsq",       None),
        ("Jsq    sheet",  "mag-Jsq-sheet", "mag-Jsq"),
        ("Jsq    smooth", "mag-Jsq-smooth","mag-Jsq"),
    ]
    names = [l for l, _ in LEGS]
    print(f"{'bin':15s} " + " ".join(f"{n:>22s}" for n in names))
    print(f"{'':15s} " + " ".join(f"{'value    share':>22s}" for _ in names))
    print("-" * (15 + 23 * len(names)))
    for lbl, key, tot in bins:
        row = f"{lbl:15s} "
        for n in names:
            v = g(n, key)
            sh = (v / g(n, tot) * 100) if tot and g(n, tot) else float("nan")
            row += f"{v:13.4e}{sh:8.2f}% " if tot else f"{v:13.4e}{'':9s} "
        print(row)

    # f_eff -- the falsifier. Only the GLOBAL integrals have int q^2 dV companions, so this is a
    # global concentration measure; a bin cannot be rescued by a split if the global is a point
    # sample AND the bin holds nearly all of it.
    print(f"\n{'f_eff':15s} " + " ".join(f"{n:>22s}" for n in names))
    for base, sq in (("dissO", "mag-dissOsq"), ("dissA", "mag-dissAsq"), ("Jsq", "mag-Jsqsq")):
        row = f"{base:15s} "
        for n in names:
            I = g(n, f"mag-{base}"); Isq = g(n, sq)
            fe = I * I / (VBOX * Isq) if Isq and Isq == Isq and Isq > 0 else float("nan")
            row += f"{fe:13.3e}{'  PT' if fe == fe and fe < 1e-3 else '    '}{'':5s} "
        print(row)
    # PER-BIN f_eff, now formable: each density bin has its own sq column and carrying volume.
    # This is what decides whether a CONVERGING bin is genuinely resolved or just a smaller
    # point sample -- the question global f_eff could not answer.
    print(f"\n{'f_eff per bin':15s} " + " ".join(f"{n:>22s}" for n in names))
    for base in ("dissO", "dissA"):
        for tag, vol in (("hi", "mag-Vhi"), ("lo", "mag-Vlo")):
            row = f"{base+'-'+tag:15s} "
            for n in names:
                I = g(n, f"mag-{base}-{tag}"); Isq = g(n, f"mag-{base}{tag}sq"); V = g(n, vol)
                fe = I * I / (V * Isq) if (Isq == Isq and Isq > 0 and V == V and V > 0) else float("nan")
                row += f"{fe:13.3e}{'  PT' if fe == fe and fe < 1e-3 else '    '}{'':5s} "
            print(row)
    row = f"{'V_hi/V_box':15s} "
    for n in names:
        row += f"{g(n,'mag-Vhi')/VBOX:13.3e}{'':9s} "
    print(row)

    print("\n=== CONVERGENCE (change per refinement rung) ===")
    print(f"{'bin':15s} {'nj4->nj8':>12s} {'nj8->nj16':>12s}  verdict")
    for lbl, key, _ in bins:
        a, b, c = (g(n, key) for n in names)
        if not (a == a and b == b and c == c) or a == 0 or b == 0:
            continue
        r1 = 100.0 * (b - a) / a
        r2 = 100.0 * (c - b) / b
        mono = (r1 < 0) == (r2 < 0)
        conv = abs(r2) < abs(r1) / 2.0
        v = ("CONVERGING" if (mono and conv and abs(r2) < 15) else
             "not monotone" if not mono else "monotone, not converged")
        print(f"{lbl:15s} {r1:11.1f}% {r2:11.1f}%  {v}")

    print("\nA bin is only a usable metric if it CONVERGES *and* is not flagged PT (point sample).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
