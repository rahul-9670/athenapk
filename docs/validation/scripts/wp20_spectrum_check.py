#!/usr/bin/env python
"""WP-20 ON-state verification of the `turb_ksample` wavenumber sampler.

WHY NOT AN FFT. The obvious check -- FFT the initial velocity field and fit E(k) --
DOES NOT WORK for this field construction, and the first attempt at it produced
nonsense (measured slopes +0.28 and +0.89 where -3.667 and -1.667 were expected).
Two reasons, both structural:
  1. `collapse_be.cpp` draws |k| from a CONTINUOUS distribution and the direction
     isotropically, so the plane waves are not commensurate with the box. A
     non-periodic field leaks power across every FFT bin, flattening (indeed
     inverting) any measured slope.
  2. `turb_region = sphere` windows the field, convolving the spectrum further.
Neither is a defect -- the IC is not meant to be a periodic turbulent box -- but
both make an FFT-based slope fit invalid here.

WHAT THIS CHECKS INSTEAD. The pgen prints `vrms_analytic = sqrt(0.5 * sum_m k_m^-alpha)`
computed from the ACTUAL drawn modes. Its distribution under each sampler is known in
closed form, so comparing the printed value against a Monte-Carlo of that distribution
tests the sampler directly and needs no periodicity:

    flat : k ~ U(kmin,kmax)                      -> dN/dk = const
    k2   : k = cbrt(kmin^3 + u(kmax^3-kmin^3))   -> dN/dk ~ k^2

Result 2026-07-31 (turb_alpha=3.667, N=128, k in [1,8], seed 42):
    flat printed 1.86246 -> 55.0th percentile of its own distribution   PASS
    k2   printed 0.54643 -> 32.8th percentile of its own distribution   PASS

It also reports the realization scatter, which is the load-bearing side effect: the
k2 sampler DOUBLES the scatter at fixed mode count (23.3% vs 11.7% at N=128), because
correct k^2 weighting puts only ~1.8 modes below k=2 while those modes carry ~49% of
the energy. N ~= 512-1024 is needed for the k2 sampler to match the status-quo scatter.
"""
import sys
import numpy as np
from scipy.integrate import quad

ALPHA, KMIN, KMAX = 3.667, 1.0, 8.0


def mc(sampler, N, alpha=ALPHA, kmin=KMIN, kmax=KMAX, trials=200000, seed=0):
    u = np.random.default_rng(seed).random((trials, N))
    k = (kmin + u * (kmax - kmin) if sampler == "flat"
         else np.cbrt(kmin**3 + u * (kmax**3 - kmin**3)))
    return np.sqrt(0.5 * (k**-alpha).sum(axis=1))


def expected(sampler, N, alpha=ALPHA, kmin=KMIN, kmax=KMAX):
    if sampler == "flat":
        f = lambda k: k**-alpha / (kmax - kmin)
    else:
        f = lambda k: k**-alpha * 3 * k**2 / (kmax**3 - kmin**3)
    return np.sqrt(0.5 * N * quad(f, kmin, kmax)[0])


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 128
    obs = {"flat": float(sys.argv[2]), "k2": float(sys.argv[3])} if len(sys.argv) > 3 else \
          {"flat": 1.86246, "k2": 0.546426}
    print(f"turb_nmodes = {N}, turb_alpha = {ALPHA}, k in [{KMIN:g}, {KMAX:g}]\n")
    print(f"{'sampler':>8} {'printed':>10} {'analytic':>10} {'MC mean':>9} {'sd':>8} "
          f"{'pctile':>8} {'verdict':>8}")
    print("-" * 68)
    for s in ("flat", "k2"):
        v = mc(s, N)
        pct = (v < obs[s]).mean() * 100
        ok = "PASS" if 2.0 < pct < 98.0 else "FAIL"
        print(f"{s:>8} {obs[s]:10.5f} {expected(s, N):10.4f} {v.mean():9.4f} "
              f"{v.std():8.4f} {pct:7.1f}% {ok:>8}")
    print("\nrealization scatter (sd/mean) vs mode count:")
    print(f"{'N':>6} {'flat':>9} {'k2':>9} {'k2 modes k<=2':>15}")
    for n in (128, 256, 512, 730, 1024, 2048):
        a, b = mc("flat", n, trials=40000, seed=1), mc("k2", n, trials=40000, seed=1)
        nlow = n * (2**3 - KMIN**3) / (KMAX**3 - KMIN**3)
        print(f"{n:6d} {a.std()/a.mean()*100:8.1f}% {b.std()/b.mean()*100:8.1f}% {nlow:15.1f}")


if __name__ == "__main__":
    main()
