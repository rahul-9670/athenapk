#!/usr/bin/env python3
"""Flagship Phase 7 ensemble — UQ analysis harness.

Turns the ensemble run outputs into the flagship deliverable: flux retention as a PREDICTIVE
DISTRIBUTION (not one curve), plus a first-order SENSITIVITY ranking (which IC parameter drives
the spread) and a separation of the IC-variance uncertainty class from the fiducial.

Method (consistent with convergence.py): each run's primary observable mu_core = M_core/Phi_core
is measured at MATCHED PHYSICAL STATE (a target rho_max, not matched time), reusing the validated
flux_retention.measure_snapshot. Across the ensemble at that state we report:
  * the distribution of mu_core (median, 16/84 percentiles, min/max) -> the predictive spread;
  * IC-variance sub-class: the spread ACROSS SEEDS at fixed LHS point (turbulence realization
    noise), separated from the parameter-driven spread;
  * first-order sensitivity: Spearman rank correlation of mu_core with each sampled IC parameter
    -> which knob dominates the flux-retention spread.

Usage:  uq.py <design_dir> [--target-rho 1.829e5] [--stride 1] [--max-log-dist 0.3]
  <design_dir> holds run_matrix.json + pointNNN/ run directories.
"""
import os, sys, json, glob, argparse
import numpy as np
sys.path.insert(0, "/beegfs/u/bbg6470/athenapk/runs")
import flux_retention as fr


def mu_at_rho(run_dir, target_rho, stride=1, max_log_dist=0.3):
    """mu_core = M_core/Phi_core at the snapshot whose rho_max is closest to target_rho (code).

    TWO BUGS FIXED 2026-08-06, both of which produced a WRONG NUMBER SILENTLY rather than an error.

    (1) `stride` defaulted to 8, which subsamples the snapshot list. That is harmless for a run
        with hundreds of snapshots and destructive for these: ensemble members that stop at the
        matched epoch carry only 7-17 snapshots, so stride=8 kept 1-3 of them. point006 (7
        snapshots) was reduced to a SINGLE candidate -- `parthenon.out1.00000.phdf`, the t=0
        initial condition, five decades below target -- and, being the only candidate, it was
        selected. Default is now 1; stride remains available for genuinely long runs.

    (2) There was NO tolerance check. `best` was whatever minimised |log10 rho - log10 target|
        over the candidates, however far that was, and `return None` fired only when nothing was
        readable at all. So a member nowhere near the target still contributed a mu_core to the
        predictive distribution. The printed message even claimed "no snapshot near rho" while the
        code never tested nearness. A member is now REJECTED unless its best snapshot lands within
        `max_log_dist` decades of the target, and the distance is reported either way so the match
        quality is visible rather than assumed.
    """
    snaps = sorted(glob.glob(os.path.join(run_dir, "parthenon.out1.*.phdf")))[::stride]
    best = None
    for s in snaps:
        try:
            rm, t = fr.rho_max_code(s)
        except Exception:
            continue                      # unreadable/truncated -> skip, never treat as a match
        dist = abs(np.log10(rm) - np.log10(target_rho))
        if best is None or dist < best[0]:
            best = (dist, s, rm)
    if best is None:
        return {"reject": "no readable snapshot", "dist": None}
    if best[0] > max_log_dist:
        return {"reject": f"closest snapshot is {best[0]:.2f} decades from target "
                          f"(rho={best[2]:.3e} code, {os.path.basename(best[1])}); "
                          f"tolerance is {max_log_dist:.2f}",
                "dist": best[0]}
    d = fr.measure_snapshot(best[1])
    if abs(d["Phi_core"]) < 1e-300:
        return {"reject": f"Phi_core is zero at {os.path.basename(best[1])} "
                          f"(r_core={d.get('r_core')})", "dist": best[0]}
    if d.get("r_core_degenerate"):
        # r_core pinned to the innermost profile bin: the rhocrit crossing was never bracketed,
        # so r_core reflects the binning rather than the core and mu_core is not a first-core
        # measurement. Excluding is the honest choice -- including it would widen the predictive
        # distribution with a discretisation artifact and misreport it as physical spread.
        return {"reject": f"degenerate core at {os.path.basename(best[1])} "
                          f"(r_core={d['r_core']:.4g} is the innermost profile bin; the rhocrit "
                          f"crossing was not bracketed)", "dist": best[0]}
    return {"mu": d["M_core"] / d["Phi_core"], "M_core": d["M_core"], "Phi_core": d["Phi_core"],
            "rho_max": best[2], "snap": os.path.basename(best[1]), "dist": best[0]}


def spearman(x, y):
    """Spearman rank correlation (no scipy dependency needed for a rank corr)."""
    x = np.asarray(x); y = np.asarray(y)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan
    rx = np.argsort(np.argsort(x[m])); ry = np.argsort(np.argsort(y[m]))
    rx = rx - rx.mean(); ry = ry - ry.mean()
    denom = np.sqrt((rx**2).sum() * (ry**2).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else np.nan


def analyze(design_dir, target_rho, stride, max_log_dist):
    meta = json.load(open(os.path.join(design_dir, "run_matrix.json")))
    params = meta["params"]
    rows = []
    for r in meta["runs"]:
        rd = os.path.join(design_dir, f"point{r['point']:03d}")
        res = mu_at_rho(rd, target_rho, stride, max_log_dist)
        if "reject" in res:
            print(f"  point{r['point']:03d}: EXCLUDED — {res['reject']}")
            continue
        print(f"  point{r['point']:03d}: matched at {res['snap']} "
              f"(rho={res['rho_max']:.3e} code, {res['dist']:.3f} decades from target)")
        rows.append({**r, **res})
    if not rows:
        print("\nNo ensemble points have reached the target state yet -- run the ensemble first.")
        print("(This harness is validated by unit test; see uq_selftest.)")
        return
    mu = np.array([x["mu"] for x in rows])
    print(f"\n=== flux-retention predictive distribution (mu_core at rho_max~{target_rho:.1e} code) ===")
    print(f"  N={len(rows)} ensemble members")
    print(f"  median={np.median(mu):.4g}  16-84pct=[{np.percentile(mu,16):.4g},{np.percentile(mu,84):.4g}]"
          f"  min/max=[{mu.min():.4g},{mu.max():.4g}]  CoV={np.std(mu)/abs(np.mean(mu)):.2%}")
    # IC-variance sub-class: spread across seeds at fixed LHS index
    byidx = {}
    for x in rows:
        byidx.setdefault(x["lhs_index"], []).append(x["mu"])
    seed_spreads = [np.std(v) / abs(np.mean(v)) for v in byidx.values() if len(v) >= 2]
    if seed_spreads:
        print(f"  IC-variance (turbulence-seed) sub-class: median across-seed CoV = "
              f"{np.median(seed_spreads):.2%}  (vs total CoV {np.std(mu)/abs(np.mean(mu)):.2%})")
    # first-order sensitivity: Spearman of mu vs each sampled parameter
    print(f"\n=== first-order sensitivity (Spearman rank corr of mu_core vs IC parameter) ===")
    sens = []
    for pn in params:
        rho_s = spearman([x[pn] for x in rows], mu)
        sens.append((pn, rho_s))
    for pn, rs in sorted(sens, key=lambda z: -abs(z[1]) if np.isfinite(z[1]) else 0):
        bar = "#" * int(abs(rs) * 30) if np.isfinite(rs) else ""
        print(f"  {pn:>16}: rho_S = {rs:+.3f}  {bar}")
    print("\n  => the top |rho_S| parameter dominates the flux-retention spread. Report mu_core")
    print("     as the distribution above, with IC-variance separated from parameter-driven spread.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("design_dir")
    ap.add_argument("--target-rho", type=float, default=1.829e5, help="target rho_max [code] (1e-13 cgs)")
    ap.add_argument("--stride", type=int, default=1,
                help="subsample snapshots; 1 = use all. >1 only for runs with "
                     "hundreds of snapshots -- it silently discards candidates.")
    ap.add_argument("--max-log-dist", type=float, default=0.3,
                help="reject a member whose closest snapshot is farther than this "
                     "many decades in rho from the target (default 0.3)")
    args = ap.parse_args()
    analyze(args.design_dir, args.target_rho, args.stride, args.max_log_dist)
