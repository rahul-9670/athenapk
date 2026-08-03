#!/usr/bin/env python
"""WP-3 — temporal (CFL) convergence of the production configuration.

Three legs at 256^3 on the new production physics (binary f181c0a1, diode fluid BCs, relative
gravity residual), differing ONLY in `parthenon/time/cfl`:

    0.3    production  -- runs/root_ladder/r256_sw   (reused, not re-run)
    0.15   half        -- runs/wp3_cfl/cfl0.15
    0.075  quarter     -- runs/wp3_cfl/cfl0.075

READ AT t = 0.90, NOT AT THE t = 1.0 ENDPOINT. WP-7 established that the last 0.01 t0 of these
uniform-grid runs sits on the collapse singularity -- dt falls 1-2 decades and each leg stalls at
a different point, so a comparison taken there measures the stall, not the parameter under test.

The observed order in TIME is p = log2(|v(0.3) - v(0.15)| / |v(0.15) - v(0.075)|). It is bounded
above by the scheme's formal 2 (WP-14 measured 2.06-2.13 for the base MHD update in space+time),
but the full production system carries operator-split physics packages that are 1st order in time
by construction, so a value between 1 and 2 is the expected outcome -- and WP-14 is what makes a
shortfall attributable rather than mysterious.

The acceptance question is NOT the order, though: it is whether production's cfl = 0.3 is already
converged, i.e. whether |v(0.3) - v(0.075)| is small against WP-18's sigma = 16 %.
"""
import os
import numpy as np

LEGS = [
    ("cfl=0.3  (production)", "/beegfs/u/bbg6470/athenapk/runs/root_ladder/r256_sw"),
    ("cfl=0.15 (half)      ", "/beegfs/u/bbg6470/athenapk/runs/wp3_cfl/cfl0.15"),
    ("cfl=0.075 (quarter)  ", "/beegfs/u/bbg6470/athenapk/runs/wp3_cfl/cfl0.075"),
]
# hst columns (verified from the header of these runs, not assumed)
C = {"time": 0, "dt": 1, "cycle": 2, "mass": 4, "KE": 8, "ME": 10, "Jsq": 22,
     "tor": 24, "pol": 25}
T_READ = 0.90
SIGMA = 16.0  # WP-18 seed-to-seed spread on MEtor/MEpol, in per cent


def load(d):
    p = os.path.join(d, "parthenon.out0.hst")
    return np.loadtxt(p) if os.path.exists(p) else None


def at(a, col, t=T_READ):
    return np.interp(t, a[:, 0], a[:, col])


def main():
    data = []
    for label, d in LEGS:
        a = load(d)
        if a is None:
            print(f"  {label}: MISSING ({d})")
            continue
        data.append((label, d, a))
    if len(data) < 2:
        print("\nNot enough legs yet.")
        return

    print(f"WP-3 — temporal convergence, matched state at t = {T_READ}\n")
    print(f"{'leg':>22} {'cycles':>8} {'t_end':>9} {'MEtor/MEpol':>13} {'ME':>12} "
          f"{'KE':>12} {'wall/cycle':>11}")
    vals = []
    for label, _d, a in data:
        r = at(a, C["tor"]) / at(a, C["pol"])
        vals.append(r)
        print(f"{label:>22} {a[-1, 2]:8.0f} {a[-1, 0]:9.5f} {r:13.6e} "
              f"{at(a, C['ME']):12.5e} {at(a, C['KE']):12.5e} {'-':>11}")

    ref = vals[0]
    print(f"\nvs production (cfl = 0.3), on MEtor/MEpol:")
    for (label, _d, _a), v in zip(data[1:], vals[1:]):
        d = (v - ref) / ref * 100.0
        print(f"  {label}: {d:+8.4f} %   ({abs(d)/SIGMA*100:.3f} % of sigma = {SIGMA} %)")

    if len(vals) >= 3:
        d1, d2 = vals[1] - vals[0], vals[2] - vals[1]
        if abs(d2) > 0:
            p = np.log2(abs(d1 / d2))
            rich = d2 / (2.0 ** p - 1.0) if 2.0 ** p > 1.0001 else float("nan")
            print(f"\nobserved order in time  p = log2(|d(0.3->0.15)| / |d(0.15->0.075)|) "
                  f"= {p:.2f}")
            print(f"Richardson residual at cfl = 0.075: {rich / vals[2] * 100:+.4f} %")
            print("  (p < 2 is EXPECTED: several physics packages are operator-split and 1st")
            print("   order in time by construction. WP-14 fixes the base MHD update at 2.06-2.13,")
            print("   which is what makes any shortfall attributable to a named package.)")

    print("\nOther observables vs production, at t = 0.90:")
    for (label, _d, a), _ in zip(data[1:], vals[1:]):
        row = "  " + label + ":"
        for k in ("ME", "KE", "Jsq", "mass"):
            x, y = at(data[0][2], C[k]), at(a, C[k])
            row += f"  d{k}={((y - x) / x * 100):+.3f} %"
        print(row)

    print("\nCost (dt at t = 0.90, total cycles, and wall time from the run log):")
    for label, d, a in data:
        wall = ""
        try:
            with open(os.path.join(d, "run.log")) as f:
                for line in f:
                    if "walltime used" in line:
                        wall = line.split("=")[-1].strip()
        except OSError:
            pass
        print(f"  {label}: dt = {at(a, C['dt']):.4e}, {a[-1, 2]:.0f} cycles to "
              f"t = {a[-1, 0]:.3f}, wall = {wall or 'n/a'} s")
    print("\n  NOTE the wall times are NOT directly comparable unless the legs ran on the same")
    print("  rank count AND without contention -- check the job records before quoting a ratio.")


if __name__ == "__main__":
    main()
