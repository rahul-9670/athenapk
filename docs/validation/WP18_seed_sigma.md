# WP-18 — seed-variance σ, the acceptance threshold for the campaign

**Status: COMPLETE. 2026-07-31.** 12 realizations, all COMPLETED, all passing the k2 IC guard.
Binary `bffdf8cd`, corrected IC (`turb_ksample=k2`, `turb_nmodes=2048`, `turb_zeta=0.6667`),
256³ uniform, t → 1.0. Raw: `runs/wp18_seed_ensemble/WP18_sigma.txt`.

## The headline

**Realization scatter at MATCHED STATE (the number to use):**

| quantity | rel σ |
|---|---|
| `MEtor/MEpol` (magnetic braking observable) | **≈ 16%** |
| `mag-Jsq` | ≈ 13–15% |
| `KE` | ≈ 12–20% |
| `mag-Hc` (current helicity) | **sign-indefinite** — spans −3.99 to +8.29 across seeds |
| `mass` | 0.00% (sanity check: the ensemble differs only in the turbulent field) |

## Matched TIME vs matched STATE — the load-bearing methodological point

At matched *time* the scatter appears to explode near collapse onset:

| t | ME | MEtor/MEpol | Jsq |
|---|---|---|---|
| 0.25 | 0.08% | 16.34% | 15.30% |
| 0.50 | 0.34% | 16.63% | 12.68% |
| 0.75 | 0.94% | 18.06% | 14.67% |
| **1.00** | **7.44%** | **65.99%** | **131.14%** |

That blow-up is **epoch mismatch, not intrinsic variance.** Seeds reach a given collapse
stage at different times: ME first reaches 50 at t = 0.898 ± 0.051 (full spread 0.156).
Comparing at matched *state* instead removes it:

| comparison basis | MEtor/MEpol rel σ |
|---|---|
| matched time t = 1.0 | 65.99% |
| matched state ME = 50 | **15.91%** |
| matched state ME = 52 | **15.79%** |

15.9% at matched state equals the 16.3% seen at t = 0.25, i.e. the intrinsic scatter is
**flat at ~16% for the whole braking phase**.

> **Always compare at matched state (density or ME), never at matched time.** WP-8 reached
> the same conclusion independently from the njeans ladder. This is now two separate lines of
> evidence for the same rule.

## What this licenses and what it kills

**The critical distinction is PAIRED vs UNPAIRED.**

*Paired* comparisons — same seed, one knob varied (resolution, CFL, STS ratio, CT vs GLM) —
are controlled experiments. The seed scatter does **not** apply to them, because the
realization is held fixed. Those comparisons remain valid at their own precision.

*Unpaired* statements — quoting a physical value, or comparing runs with different ICs — carry
the full ~16% error bar.

Consequences:

1. **The CT-vs-GLM 1.6%** was a paired comparison, so it is a valid *paired difference*. But
   it is ~10× below the realization scatter, so it **cannot be reported as a physical result**.
   The honest statement is: "at fixed IC, dropping CT changes flux retention by 1.6%, an order
   of magnitude below the ~16% realization scatter."
2. **WP-8's njeans conclusion** (nj8 ↔ nj16 agree to 1.2% on MEtor/MEpol) is likewise paired —
   all three legs used seed 42 — so "njeans = 8 is converged" stands *at that seed*. Whether
   it generalizes is untested; that would need the ladder repeated at ≥ 3 seeds.
3. **Any quoted physical value** — e.g. MEtor/MEpol ≈ 1.4e-2 — must carry ± ~16%.
4. **`mag-Hc` is sign-indefinite across realizations.** Current helicity cannot be reported as
   a single-realization result at all; only its ensemble distribution is meaningful.

## Why the scatter is this large

Traceable to the IC construction (WP-20). With the corrected k² sampler, the energy-bearing
low-k shells are sparsely populated: at `turb_nmodes = 2048` only ~28 modes lie at k ≤ 2, and
those carry roughly half the turbulent energy. The largest scales — which set the collapse
geometry — are therefore determined by a few tens of random draws. (At the old `nmodes = 128`
it would have been ~2 modes; the MC estimate of scatter in `vrms_analytic` alone was 23.3%,
versus 6.0% at 2048.) Raising `turb_nmodes` further would reduce it at one-off IC cost.

## Practical acceptance rule for the remaining WPs

For WP-1 (RKL2 ratio), WP-2 (creduc), WP-3 (CFL), WP-7 (root grid), WP-9 (n_group): these are
all **paired** sweeps at fixed seed, so judge them against their own run-to-run reproducibility,
not against 16%. But **any conclusion phrased as a statement about the physical system** must
be checked against the 16% band, and ideally repeated at ≥ 3 seeds before it goes in the paper.
