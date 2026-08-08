#!/usr/bin/env python
"""WP-8 — does the CURRENT-SHEET split restore convergence for Jsq, where density did not?

WHAT THIS CLOSES. `wp8_split_convergence.py` (2026-08-06) applied the DENSITY split to the njeans
ladder and found it FAILS for Jsq: the low-density bin carries 97-99 % of the integral and keeps
f_eff ~ 1e-7, so `Jsq-lo` is the original pathology under a new name (-51.6 %/-39.9 % vs the
global -52.3 %/-39.3 %). Its conclusion was explicit: *"Jsq needs a split on a current-sheet
indicator, not on density"*. The indicator was implemented in the code on 2026-08-08
(`hydro/mag_diag_sheet_thresh`, columns mag-Jsq-sheet/smooth, mag-Vsheet, mag-Jsqsq) but, exactly
like the density split before it, was never MEASURED on the ladder. This script measures it.

THE INDICATOR, and why it is the right variable.
    s = |J| * dx / |B|
= the fraction of the local field that reverses across one cell. s -> 1 means B flips over a
single zone: the current is at the grid scale and is a resolution artefact as much as a physical
structure. s << 1 means a current spread over many cells -- resolved, physical. It is measured IN
UNITS OF THE GRID, which is precisely the thing that changes between ladder rungs; density is
blind to it, which is why a density threshold cannot separate grid-scale sheets that live in the
diffuse envelope.

THE FALSIFIER, carried over deliberately. A split can look convergent for a trivial reason: if
one bin holds essentially the whole integral, it just reproduces the global number and the other
is noise. So every bin reports BOTH its share of the global integral AND its f_eff. A "converged"
bin holding 99 % at f_eff ~ 1e-7 is the pathology relabelled, not a fix. Both splits are printed
side by side on the SAME snapshots so the comparison is like-for-like.

METHOD is identical to wp8_split_convergence.py (same stencil, same interior-only differencing,
same matched-density epoch selection from epoch_scan.txt) so the two results are directly
comparable. Reused verbatim where possible rather than reimplemented.

    wp8_sheet_convergence.py [rho_target_cgs] [thresh ...]      default 1e-12, 0.1 0.3 0.5
"""
import sys, os
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wp8_split_convergence import (RHO0, LEGS, LADDER, RHO_SPLIT_CODE,
                                   curl_sq_interior, f_eff, scan)


def measure_both(path, threshes):
    """Density split AND current-sheet split from one read of one snapshot."""
    with h5py.File(path, "r") as h:
        names = [str(s) for s in h["Info"].attrs["ComponentNames"]]
        off = names.index("prim_density")
        ip = lambda n: names.index(n) - off
        p = h["prim"]
        rho = p[:, ip("prim_density"), ...]
        b1 = p[:, ip("prim_magnetic_field_1"), ...]
        b2 = p[:, ip("prim_magnetic_field_2"), ...]
        b3 = p[:, ip("prim_magnetic_field_3"), ...]
        xf = h["Locations/x"][:]
    dxb = xf[:, 1] - xf[:, 0]
    jsq = curl_sq_interior(b1, b2, b3, dxb)
    sl = (slice(None), slice(1, -1), slice(1, -1), slice(1, -1))
    rho_i = rho[sl]
    bmag = np.sqrt(b1[sl] ** 2 + b2[sl] ** 2 + b3[sl] ** 2)
    dV = np.broadcast_to((dxb ** 3)[:, None, None, None], jsq.shape)
    Vtot = float(dV.sum())

    out = {"V": Vtot}
    fe, I = f_eff(jsq, dV, Vtot)
    out["glob"] = (I, fe)

    # --- density split (the one already shown to fail) ---
    hi = rho_i > RHO_SPLIT_CODE
    for lab, m in (("rho_hi", hi), ("rho_lo", ~hi)):
        if m.sum() == 0:
            out[lab] = (0.0, float("nan"), 0.0); continue
        fe_b, I_b = f_eff(jsq[m], dV[m], float(dV[m].sum()))
        out[lab] = (I_b, fe_b, float(dV[m].sum()) / Vtot)

    # --- current-sheet split ---
    # |B| == 0 with J != 0 is a pure grid artefact => counts as sheet (s = inf). Matches the
    # kernel in src/diagnostics/mag_diag.cpp exactly.
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(bmag > 0.0, np.sqrt(jsq) * dV ** (1.0 / 3.0) / bmag, np.inf)
    for th in threshes:
        sheet = s > th
        for lab, m in ((f"sheet@{th}", sheet), (f"smooth@{th}", ~sheet)):
            if m.sum() == 0:
                out[lab] = (0.0, float("nan"), 0.0); continue
            fe_b, I_b = f_eff(jsq[m], dV[m], float(dV[m].sum()))
            out[lab] = (I_b, fe_b, float(dV[m].sum()) / Vtot)
    return out


def main():
    rho_t = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-12
    threshes = [float(x) for x in sys.argv[2:]] or [0.1, 0.3, 0.5]
    print(f"WP-8 CURRENT-SHEET split on the njeans ladder — matched epoch "
          f"rho_max = {rho_t:.3e} g/cm^3\n")
    res = {}
    for leg in LEGS:
        got = scan(leg, rho_t)
        if got is None:
            print(f"  {leg}: no snapshot found"); return
        d, path, rmax_code, t = got
        print(f"  {leg:5s} {os.path.basename(path):28s} t={t:.5f} "
              f"rho_max={rmax_code*RHO0:.3e} g/cm3 ({d:.3f} dex from target)")
        res[leg] = measure_both(path, threshes)
    print()

    labels = ["glob", "rho_hi", "rho_lo"] + \
             [f"{k}@{th}" for th in threshes for k in ("sheet", "smooth")]
    print(f"{'bin':14s} " + " ".join(f"{L:>26s}" for L in LEGS))
    print(f"{'':14s} " + " ".join(f"{'I    share   f_eff':>26s}" for _ in LEGS))
    print("-" * (14 + 27 * len(LEGS)))
    for lab in labels:
        row = f"{lab:14s} "
        for L in LEGS:
            v = res[L].get(lab)
            if v is None: row += f"{'--':>26s} "; continue
            I = v[0]; fe = v[1]
            share = I / res[L]["glob"][0] if res[L]["glob"][0] else float("nan")
            row += f"{I:9.3e}{share*100:7.2f}%{fe:10.2e} "
        print(row)

    print("\n=== CONVERGENCE (change per refinement rung) ===")
    print(f"{'bin':14s} {'nj4->nj8':>12s} {'nj8->nj16':>12s}  verdict")
    for lab in labels:
        a, b, c = (res[L].get(lab, (0.0,))[0] for L in LEGS)
        if not a or not b:
            continue
        r1 = 100.0 * (b - a) / a
        r2 = 100.0 * (c - b) / b if b else float("nan")
        mono = (r1 < 0) == (r2 < 0)
        conv = abs(r2) < abs(r1) / 2.0
        v = ("CONVERGING" if (mono and conv and abs(r2) < 15) else
             "not monotone" if not mono else "monotone, not converged")
        # a bin is only meaningful if it is not itself a point sample
        fe = res["nj16"].get(lab, (0, float("nan")))[1]
        if fe == fe and fe < 1e-3:
            v += f"  [f_eff={fe:.1e} — POINT SAMPLE, treat with suspicion]"
        print(f"{lab:14s} {r1:11.1f}% {r2:11.1f}%  {v}")


if __name__ == "__main__":
    main()
