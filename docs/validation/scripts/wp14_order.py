#!/usr/bin/env python
"""WP-14 — observed order of accuracy from the linear-wave ladder.

Reads runs/wp14_order/wf<W>_n<N>/linearwave-errors.dat and reports, per wave family, the
observed order p between successive rungs:  p = log2( L1(h) / L1(h/2) ).

PLM + vl2 is formally 2nd order, so the acceptance target is p -> 2 as h -> 0. A ladder that
sits at p ~ 1 would mean the production configuration is effectively first order, which is the
outcome WP-14 exists to rule in or out.

Note on the error file: column 5 (0-indexed 4) is the RMS L1 over all conserved variables;
the per-variable L1 columns follow. We report the RMS as the headline and density separately,
because a scheme can be 2nd order in the smooth variables while a single component (e.g. a
cell-centred B under GLM) drags the RMS down -- and that distinction matters for the paper.
"""
import os, sys
import numpy as np

HERE = "/beegfs/u/bbg6470/athenapk/runs/wp14_order"
NAMES = {0: "fast-", 1: "Alfven-", 2: "slow-", 3: "entropy", 4: "slow+", 5: "Alfven+", 6: "fast+"}
COLS = {"RMS": 4, "rho": 5, "M1": 6, "E": 9, "B2c": 11}


def read(path):
    with open(path) as f:
        rows = [l.split() for l in f if not l.startswith("#") and l.strip()]
    return [float(x) for x in rows[-1]]


def main():
    waves = [int(w) for w in (sys.argv[1:] or [0, 1, 2, 3])]
    ns = [16, 32, 64, 128]
    any_found = False
    for wf in waves:
        data = {}
        for n in ns:
            p = os.path.join(HERE, f"wf{wf}_n{n}", "linearwave-errors.dat")
            if os.path.exists(p):
                data[n] = read(p)
        if len(data) < 2:
            continue
        any_found = True
        print(f"\n=== wave_flag={wf} ({NAMES.get(wf,'?')}) ===")
        hdr = "  N     " + "".join(f"{k:>13}" for k in COLS) + "   |   " + \
              "".join(f"p({k}){'':>3}" for k in COLS)
        print(hdr)
        got = sorted(data)
        for i, n in enumerate(got):
            row = "".join(f"{data[n][c]:13.4e}" for c in COLS.values())
            if i == 0:
                ords = "".join(f"{'--':>8}" for _ in COLS)
            else:
                prev = data[got[i - 1]]
                ords = "".join(
                    f"{np.log2(prev[c]/data[n][c]):8.2f}" if data[n][c] > 0 else f"{'--':>8}"
                    for c in COLS.values())
            print(f"{n:5d} {row}   |   {ords}")
        # asymptotic order = between the two finest rungs available
        if len(got) >= 2:
            a, b = got[-2], got[-1]
            p_rms = np.log2(data[a][4] / data[b][4])
            verdict = "PASS (>=1.8)" if p_rms >= 1.8 else (
                "MARGINAL (1.5-1.8)" if p_rms >= 1.5 else "FAIL (<1.5)")
            print(f"  asymptotic p(RMS) over {a}->{b}: {p_rms:.3f}   {verdict}")
    if not any_found:
        print("no completed ladders found under", HERE)


if __name__ == "__main__":
    main()
