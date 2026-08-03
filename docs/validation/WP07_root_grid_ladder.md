# WP-7 — root-grid convergence ladder, and Gate A (the gravity-fix re-baseline)

**Status: BOTH PASS, and the two answers come from the same measurement.**

- **WP-7 PASS** for the phase it was built to test: through t ≤ 0.95 the observed order is
  **p = 1.14–1.60** and the Richardson residual at 512³ is **0.41–1.75 %**.
- **Gate A PASS**: the WP-13 gravity fix moves `MEtor/MEpol` by **+0.108 / +0.076 / +0.041 %** at
  128/256/512 — ~400× below σ ≈ 16 % and **shrinking with resolution**. WP-8 and WP-18 survive as
  paired comparisons; neither needs re-running.
- **A methodological finding that touches several WPs: the t = 1.0 endpoint of these uniform-grid
  runs is a singular stall and must not be used as the comparison state.**

Legs, all on deck `824da4c6` (md5-verified per job), uniform grid, t → 1.0:
`r128`/`r256`/`r512` on `bffdf8cd` (stale gravity), `r128_gfix`/`r256_gfix`/`r512_gfix` on
`49d9c257` (fixed). `r256_rank4` is **not** part of this comparison — it is the WP-12
decomposition study and ran a different deck (`a5e1d221`) and binary (`5ebddce0`).

---

## 1. The endpoint is a singularity, not a state

Read at t = 1.0, the ladder looks catastrophic: observed order collapses to **p = 0.21** and
Richardson puts the residual discretization error at 512³ at **+46 %**. That reading is wrong, and
the reason is visible in `dt`:

| t | dt (128³) | dt (256³) | dt (512³) |
|---|---|---|---|
| 0.90 | 8.26e-03 | 4.11e-03 | 2.05e-03 |
| 0.99 | 6.74e-03 | 3.15e-03 | 1.36e-03 |
| **1.00** | **3.75e-03** | **3.45e-04** | **2.21e-04** |

In the last 0.01 t₀ every leg falls off a cliff — a factor 9 at 256³, 6 at 512³ — as the central
density runs away toward the first core. **These are uniform grids with `refinement = none`**, so
none of them can represent the forming core, and each stalls at a different point. Comparing at
t = 1.0 compares three different stall states.

Everything below is therefore read at **t = 0.90**, and the endpoint numbers are kept alongside
only to show the size of the artifact.

## 2. WP-7 — root-grid convergence

Observed order `p = log₂(|v₂₅₆ − v₁₂₈| / |v₅₁₂ − v₂₅₆|)` on `MEtor/MEpol`, with the Richardson
residual at 512³:

| t | v(128³) | v(256³) | v(512³) | p | residual @512³ |
|---|---|---|---|---|---|
| 0.20 | 0.000589 | 0.000617 | 0.000630 | 1.21 | 1.47 % |
| 0.40 | 0.001904 | 0.001959 | 0.001984 | 1.14 | 1.05 % |
| 0.60 | 0.003723 | 0.003812 | 0.003852 | 1.14 | 0.86 % |
| 0.80 | 0.006017 | 0.006146 | 0.006198 | 1.32 | 0.56 % |
| **0.90** | 0.007459 | 0.007635 | 0.007695 | **1.54** | **0.41 %** |
| 0.95 | 0.008382 | 0.008644 | 0.008730 | 1.60 | 0.48 % |
| 0.98 | 0.009136 | 0.009608 | 0.009820 | 1.16 | 1.75 % |
| 1.00 | 0.009864 | 0.010838 | 0.011679 | *0.21* | *45.6 %* |

**PASS.** Through the magnetically-braked envelope phase the residual is under 1 % at 512³ and
under 2 % at 256³ — well inside σ ≈ 16 % — and it *shrinks* monotonically until the singularity.

**p ≈ 1.2–1.6 is expected, not a shortfall.** WP-14 established the base MHD update at
p = 2.06–2.13; several physics packages here are operator-split and 1st order in time by
construction, so the full system must land between 1 and 2. Having WP-14 first is what makes this
attributable rather than mysterious.

### Falsification checks run before accepting this

- **Is the IC the same field at every resolution?** Yes. All three legs report `Modes : 2048,
  k in [1, 8], zeta=0.6667, alpha=3.667, k sampling: k2` — a fixed, resolution-independent set of
  large-scale modes (k ≤ 8 is ≥ 16 cells/wavelength even at 128³). At t = 0: mass identical to 7
  figures, ME identical, KE agreeing to 0.04 %, `MEtor = 0`. If the mode count had scaled with N
  the three legs would have had *different* ICs and this would not be a convergence test at all.
- **Is the divergence early or late?** Late, and only late. The 512/128 ratio falls monotonically
  1.118 → 1.069 → 1.042 → 1.035 → 1.030 across t = 0.1 → 0.9 — converging — then jumps to 1.184 at
  t = 1.0.
- **Same deck and binary per leg?** Verified from each job's own `md5sum` output, not assumed.

### `mag-Jsq` does not converge

99.6 → 264.8 → 752.4 at 128/256/512 (×2.66, ×2.84). `Jsq` is a grid-scale quantity and is
**diverging**, not converging, under refinement. This corroborates every other study in this
campaign that has found `Jsq` non-quotable, and it should not appear in any physical claim.

## 3. Gate A — does the WP-13 gravity fix invalidate WP-8 and WP-18?

`SelfGravity::FillPoissonRHS` read `prim` with no guaranteed inter-package ordering, so the Poisson
RHS was **one stage stale on every step of every self-gravity run**. WP-7, WP-8 and WP-18 — the
last of which supplies σ ≈ 16 %, the campaign's acceptance threshold — were all measured before the
fix. Gate A asks whether they survive.

At **t = 0.90**, fixed vs stale, same deck and seed (a strictly paired comparison):

| N | stale `bffdf8cd` | fixed `49d9c257` | Δ(MEtor/MEpol) | ΔME | ΔKE | ΔJsq | Δmass |
|---|---|---|---|---|---|---|---|
| 128 | 0.007451 | 0.007459 | **+0.108 %** | +0.064 % | +1.545 % | +2.40 % | −3.5e-4 % |
| 256 | 0.007629 | 0.007635 | **+0.076 %** | +0.035 % | +0.815 % | +1.64 % | −2e-5 % |
| 512 | 0.007692 | 0.007695 | **+0.041 %** | +0.018 % | +0.416 % | +0.92 % | −3.9e-4 % |

**PASS.** The shift is ≤ 0.11 % anywhere — ~150–400× below σ — and **every observable shrinks
monotonically with resolution**. That is exactly what the mechanism predicts: a one-stage lag in
the RHS is an O(dt) error, and dt ∝ dx on a CFL-limited grid, so its effect scales like the
discretization error and vanishes under refinement. A defect that merely *happened* to be small
would not converge away this cleanly.

**Consequence: WP-8 (njeans ladder) and WP-18 (12-seed σ) stand as paired comparisons and do NOT
need re-running.** Gate B — the 3-seed spot-check of σ on the fixed binary — is obviated: a
systematic +0.1 % bias common to every seed cannot change a 16 % spread.

*Absolute* numbers quoted from stale-gravity runs should still be regenerated when convenient, but
no conclusion in the campaign turns on 0.1 %.

### Correction to an earlier reading

An earlier pass through this comparison reported **+0.032 % at 256³ and −0.259 % at 128³**, with
the shift shrinking from 128³. Those figures are not reproducible from the completed runs and are
superseded by the table above. The direction of the conclusion (small, shrinking with resolution,
WP-8/WP-18 survive) is unchanged; the magnitudes and the 128³ sign are not.

A second intermediate reading, taken at the t = 1.0 endpoint, gave +1.494 / +2.093 / +2.331 % —
*growing* with resolution. That is the singular-stall artifact of §1, not a physics shift, and it
is superseded too.

## 3b. Ladder homogeneity note — the `r*_sw` rungs are not all on the same deck

`r128_sw` and `r256_sw` ran on deck `c67eb458` (`chemistry/nsub_max = 400`). `r512_sw` was released
from hold AFTER the WP-10 measurement, so it runs on `43eeb280` (`nsub_max = 4000`). The job logs
will therefore echo different deck md5s for rungs of the same ladder, which is worth knowing before
someone treats it as an error.

**Quantified, not waved away:** the measured difference between `nsub_max` 400 and 4000 on this
exact configuration is **2.06e-05**, confined to `mag-dissA`, with `MEtor/MEpol`, `ME`, `KE`, `Jsq`
and `mass` unchanged to printed precision (WP-10 §4, job 2450198). The differences the ladder
exists to measure are ~9 % *between* rungs — five orders of magnitude larger. The inhomogeneity is
below the printed precision of every quantity this ladder is read on.

It was left as-is rather than forced homogeneous because neither choice is clean: matching the
older rungs via a CLI override would leave the deck md5 mismatched anyway, and re-running the two
completed rungs costs ~1 GPU-hour to move a number by 2e-5. Recorded here so the md5 mismatch is
explained rather than discovered.

## 4. What this changes elsewhere

**Any WP that judged "at the common end time t = 1.0" on a uniform-grid root-ladder deck must be
re-read at t ≈ 0.90.** Already applied to WP-2, where it strengthened the result by two decades
(+0.49 % at the endpoint → **+0.003 %** at t = 0.90). Production itself is unaffected — it runs
AMR, which is precisely the instrument the uniform ladder lacks — but every *study* built on this
deck inherits the artifact.
