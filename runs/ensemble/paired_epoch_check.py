#!/usr/bin/env python3
"""Isolate the EPOCH effect from the SAMPLE effect in the two-epoch mu_core comparison.

The 1e-13 pass keeps 10 members and the 1e-12 pass keeps 5, overlapping in only point000 and
point001. So comparing the two medians (48.96 vs 97.13) conflates "mu_core changes between
rho_crit and 10x rho_crit" with "these are different members". The only clean statement available
from the existing data is the WITHIN-MEMBER change on the two members measured at both epochs.
"""
import sys
sys.path.insert(0, "/beegfs/u/bbg6470/athenapk/runs")
import flux_retention as fr

pairs = {
    "point000": ("parthenon.out1.00007.phdf", "parthenon.out1.00013.phdf"),
    "point001": ("parthenon.out1.00008.phdf", "parthenon.out1.00013.phdf"),
}
print(f"{'member':>9} {'epoch':>6} {'rho_max[code]':>14} {'M_core':>12} {'Phi_core':>12} {'mu_core':>10}")
res = {}
for m, (lo, hi) in pairs.items():
    for tag, f in (("1e-13", lo), ("1e-12", hi)):
        p = f"/beegfs/u/bbg6470/athenapk/runs/ensemble/design01/{m}/{f}"
        d = fr.measure_snapshot(p)
        mu = d["M_core"] / d["Phi_core"]
        res[(m, tag)] = mu
        print(f"{m:>9} {tag:>6} {d['rho_max']:14.4e} {d['M_core']:12.5e} {d['Phi_core']:12.5e} {mu:10.4f}")
print()
for m in pairs:
    a, b = res[(m, "1e-13")], res[(m, "1e-12")]
    print(f"  {m}: mu_core {a:.4f} -> {b:.4f}   ratio {b/a:.3f}  ({100*(b/a-1):+.1f}% from 1x to 10x rho_crit)")
