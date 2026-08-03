#!/usr/bin/env python
"""WP-10 follow-up — the GPU cost and physics of raising `chemistry/nsub_max`.

WHY THIS EXISTS. runs/wp10_chem/ measured nsub_max 400 / 3200 / 100000 at identical wall time on a
32^3 CPU deck and concluded "free". That does not transfer to production, for two reasons the CPU
test cannot see:

  1. COST MODE. The sub-cycler is a SERIAL per-cell loop inside a Kokkos kernel. On a GPU the
     lanes of a warp advance in lockstep, so a warp costs its SLOWEST cell: one cell wanting 1e4
     sub-steps makes its 31 neighbours pay 1e4 too. One cell per thread on CPU, no lockstep, no
     penalty.
  2. REGIME. Production reaches rho/rho_0 ~ 1e10 with x_e collapsing through the C->CO ionization
     minimum; the sub-step count the accuracy criterion demands is a function of local state.

PHYSICS IS NOT AT RISK, and that is established by inspection rather than measurement: the update
is y_new = (y + dt*P)/(1 + dt*L), linearized-implicit, unconditionally stable and
positivity-preserving at ANY step size. A smaller dt_floor can only improve accuracy. So the only
question is cost -- plus the size of the physics SHIFT, which tells us how wrong nsub_max = 400 was.

Read at t = 0.90, never the t = 1.0 endpoint (WP-7: that endpoint is a collapse singularity).
"""
import os
import re
import numpy as np

HERE = "/beegfs/u/bbg6470/athenapk/runs/wp10_chem"
LEGS = [400, 4000, 40000]
C = {"time": 0, "dt": 1, "cycle": 2, "mass": 4, "KE": 8, "ME": 10, "Jsq": 22,
     "tor": 24, "pol": 25}
T_READ = 0.90
SIGMA = 16.0


def at(a, col, t=T_READ):
    return np.interp(t, a[:, 0], a[:, col])


def wsec_step(log):
    """Median wsec_step over the run -- robust to the startup and output-dump outliers that
    make a single cycle's value or a raw total misleading."""
    v = []
    with open(log) as f:
        for line in f:
            # NOT r"wsec_step=..." -- each cycle line carries BOTH
            #   zone-cycles/wsec_step=7.14e+05   (throughput)
            #   wsec_step=2.94e+00              (the per-step wall time we want)
            # and re.search returns the FIRST match, i.e. the throughput. That silently
            # reported 1.02e6 "seconds per step". Require the token not to be preceded by
            # '/' or a word character.
            m = re.search(r"(?<![/\w])wsec_step=([0-9.eE+-]+)", line)
            if m:
                v.append(float(m.group(1)))
    return (np.median(v), len(v)) if v else (float("nan"), 0)


def trunc_cells(log):
    """Cells reported by the B10 instrument as nsub_max-limited (first occurrence; the warning
    prints once per run per rank)."""
    tot = None
    with open(log) as f:
        for line in f:
            m = re.search(r"WARNING Chemistry: (\d+) cell", line)
            if m:
                tot = (tot or 0) + int(m.group(1))
    return tot


def main():
    data = []
    for n in LEGS:
        # `seq_ns*` = the three legs run BACK TO BACK IN ONE JOB on the same 2 GPUs, which is
        # what makes wsec_step comparable between them. `gpu_ns*` were the earlier concurrent
        # 3-GPU legs (only ns400 completed) and are NOT comparable across nodes -- see WP-3's
        # 1.02e6-vs-1.76e6 throughput split for why.
        d = os.path.join(HERE, f"seq_ns{n}")
        h, lg = os.path.join(d, "parthenon.out0.hst"), os.path.join(d, "run.log")
        if not os.path.exists(h):
            print(f"  nsub_max={n}: MISSING")
            continue
        data.append((n, np.loadtxt(h), lg))
    if not data:
        print("no legs yet")
        return

    print(f"WP-10 follow-up — nsub_max on the PRODUCTION config (128^3, 3 GPUs), t = {T_READ}\n")
    print(f"{'nsub_max':>9} {'cycles':>7} {'wsec_step':>11} {'vs 400':>8} "
          f"{'floored cells':>14} {'MEtor/MEpol':>13}")
    base_w, base_v = None, None
    for n, a, lg in data:
        w, _ = wsec_step(lg)
        v = at(a, C["tor"]) / at(a, C["pol"])
        tc = trunc_cells(lg)
        if base_w is None:
            base_w, base_v = w, v
        print(f"{n:9d} {a[-1,2]:7.0f} {w:11.4f} {w/base_w:7.3f}x "
              f"{(str(tc) if tc is not None else '0'):>14} {v:13.6e}")

    print(f"\nPhysics shift vs nsub_max = 400 (production), on MEtor/MEpol:")
    for n, a, _ in data[1:]:
        v = at(a, C["tor"]) / at(a, C["pol"])
        d = (v - base_v) / base_v * 100.0
        print(f"  nsub_max={n:6d}: {d:+8.4f} %   ({abs(d)/SIGMA*100:.3f} % of sigma = {SIGMA} %)")

    print("\nOther observables vs nsub_max = 400, at t = 0.90:")
    for n, a, _ in data[1:]:
        row = f"  nsub_max={n:6d}:"
        for k in ("ME", "KE", "Jsq", "mass"):
            x, y = at(data[0][1], C[k]), at(a, C[k])
            row += f"  d{k}={((y - x) / x * 100):+.3f} %"
        print(row)

    print("""
HOW TO READ IT.
  * cost ratio < ~1.1  -> raise nsub_max; cfl_cool then becomes the binding tolerance it was
    always meant to be, and the chemistry is converged for free.
  * cost ratio large   -> the answer is an intermediate nsub_max, not 400 and not 40000. Pick the
    largest value whose cost is acceptable and report the residual error honestly.
  * physics shift large -> a finding in its own right: production's chemistry has never been
    converged, and every flux-retention number carries that.
  * floored cells: 100 % of the domain at nsub_max = 400 means cfl_cool is INERT there (B10).""")


if __name__ == "__main__":
    main()
