# WP-13b — restart reproducibility in the PRODUCTION configuration (GPU, 4 ranks, AMR, RKL2)

**Status: CLOSED 2026-08-08 with a FIX. Two independent defects, in two different code paths.
The one that reaches the FLAGSHIP is `DiodeBC`: it filled the domain-boundary ghost zones of
`cons` only, so every M1 radiation moment had NO boundary condition at all. Fixed in
`src/bvals/boundary_conditions_apk.hpp`; after the fix the flagship restart is reproducible to
the non-determinism floor (worst field 1.48x the floor, was 1.35e+10x). The second defect is in
the GRAY matter coupling, which wrote an unconverged Newton iterate; it is off the flagship path
(all 69 production decks set `n_group = 3`) and is fixed separately in
`src/radiation/radiation_moments.cpp`.**

## Why this exists

WP-13 closed restart reproducibility on `build_cpu`, one rank, `diffusion/integrator = unsplit`,
12 cycles. Production is GPU + 4 MPI ranks + AMR + RKL2. The trigger was ensemble member point012
(2026-08-08): a fresh leg and a leg restarted from the t=0 restart file reached cycle 200 with
`t = 1.0262248` vs `0.9922842`. GPU reductions need not be order-deterministic, so the test below
always measures the non-determinism floor explicitly.

## Design

Three legs, same deck, same binary, same rank count:

| leg | what it runs |
|---|---|
| `fresh_a` | `-i <deck>`, nlim=N |
| `fresh_b` | `-i <deck>`, nlim=N — IDENTICAL to `fresh_a` |
| `split` | `-i <deck>` nlim=N/2, then `-r <rst>` nlim=N |

`fresh_a` vs `fresh_b` is the **non-determinism floor**; `fresh_a` vs `split` is the restart
divergence. Without leg B the test cannot tell a restart defect from ordinary non-determinism.
Harness: `runs/wp13b_gpu_restart/submit_wp13b.sh` (parameterised by
`WP13B_DIR/DECK/BIN/DIFF/NLIM/SPLIT/FREEZE/EXTRA`) and `compare_wp13b.py`.

### Traps this test walked into, in order

* **The first run was VACUOUS.** The deck's `tlim = 1.5` is reached at cycle 22, before the
  nlim=30 split, so the "restart" began already at `tlim` and advanced ZERO cycles in 0.108 s.
  All three legs compared bit-identical for a trivial reason. Fixed with `tlim=1000` plus positive
  checks that every leg ends at `nlim` and the restart leg STARTS at the split cycle.
* `rkl2_freeze_eta=true` hard-fails unless `integrator=rkl2`, so an `unsplit` leg must also pass
  `WP13B_FREEZE=false`. `physics/radiation=false` is not usable at all: `eos=hydrogen requires
  <physics> radiation=true` (`hydro.cpp:1001`).
* **The deck under test was not the flagship's radiation configuration.** `runs/wp13_restart/
  straight/fhc.in` sets no `n_group`, so it runs the GRAY `MatterCoupling`. The flagship and all
  24 ensemble members set `n_group = 3` and run `MatterCouplingMultigroup`, a *different solver*.
  Everything measured on the gray deck says nothing about the flagship, and vice versa. The
  harness now takes `WP13B_DECK` so the production deck can be tested directly.
* **`self_gravity/residual_tolerance` is silently ignored.** The real key is
  `self_gravity/solver_params/residual_tolerance`; the short form is accepted and then reported
  as "set but unused", which made one experiment a no-op that looked like a null result.
* **The ensemble deck dumps `out1` in single precision.** `single_precision_output = 1` hides
  every difference below ~1e-7 relative; the flagship legs must override it to 0.
* **`compare_wp13b.py`'s own verdict was wrong.** It reduced each comparison to one global
  `max|rel|`. That statistic SATURATES at 2.0 whenever two values have opposite signs, so the
  near-empty `rad.Fr*_g2` fields (magnitude ~1e-15) pinned BOTH the floor and the restart
  comparison at ~1.9996 and the ratio test returned **PASS** — while `rad.Fr1` was diverging by
  74.1 against a floor of 8.6e-09. The verdict is now per field, on absolute scales, and requires
  a flagged field to exceed both the floor (>10x) and 1e-4 of its own magnitude.

---

## Defect 1 (FLAGSHIP) — `DiodeBC` never filled the radiation ghost zones

`DiodeBC` packed only `"cons"`:

```cpp
auto cons = mbd->PackVariables(std::vector<std::string>{"cons"}, coarse);
```

Parthenon's stock `outflow` uses `GenericBC<..., variable_names::any>` and copies **every**
variable. So with `ix1_bc = diode` the domain-boundary ghosts of every other `FillGhost` field —
all of `rad.Er`, `rad.Fr1..3` and their per-group copies — were written by no boundary condition
at all. A FRESH run got away with it because the problem generator initialises the whole array,
ghosts included, so ghost ~ interior. A RESTART does not: restart files store INTERIOR cells only,
so the ghosts came up zero-initialised and the first `CalculateRadFluxes` after the restart saw a
full-amplitude jump across the domain face.

**42 decks use `diode`, and all 42 also enable radiation** — including
`runs/root_ladder/fhc_rootladder.in` and all 24 `runs/ensemble/design01/point*/fhc_ens.in`.

### Measurement — flagship deck (ensemble point000), 4x H100, 256^3, 31 cycles, split at 30

| field | floor max&#124;abs&#124; | restart max&#124;abs&#124; | ratio |
|---|---|---|---|
| grav.phi | 2.583e-09 | 2.279e-09 | 0.88 |
| prim | 3.231e-07 | **4.458e-02** | **1.4e+05** |
| rad.Er | 1.902e-10 | **9.772e-02** | **5.1e+08** |
| rad.Fr1 | 8.618e-09 | **7.414e+01** | **8.6e+09** |
| rad.Fr2 | 3.291e-08 | **7.414e+01** | **2.3e+09** |
| rad.Fr3 | 5.510e-09 | **7.414e+01** | **1.4e+10** |

The divergent cells are **exactly the outermost layer of meshblocks** (L-infinity radius 6.031 to
7.969, i.e. 8.0 minus one 32-cell block) and nothing inside it. The magnitude is predicted by the
free-streaming flux a zeroed ghost produces: `chat*Er/2 = 1578*0.104/2 ~ 79` vs 74.1 measured.
Block ordering is identical (`LogicalLocations` and `Levels` match, 512 blocks, all level 0), so
the comparison is cell-for-cell.

### Isolation

| experiment | job | result |
|---|---|---|
| radiation package OFF (`RAD_DISABLE_TRANSPORT=1`) | 2491808 | **PASS** — worst field 1.11x the floor |
| `diode` -> `outflow`, unfixed binary | 2491862 | rad.Fr1 7.414e+01 -> **2.66e-09** (2.8e10x better) |
| **`DiodeBC` FIXED**, `diode` kept | 2491892 | **PASS** — every field 0.42–1.48x the floor |

### The fix

`DiodeBC` now calls Parthenon's own `GenericBC<DIR, SIDE, BCType::Outflow,
variable_names::any>` first — bit-for-bit the stock `outflow` copy, for every field — and then
applies the existing normal-momentum clamp to `cons`. The `cons` path is unchanged (same copy
from the same reference cell, then the same clamp), so a radiation-free run is bit-identical.
This is the behaviour `outflow` decks already had and which is already validated.

### Impact on results already taken

The fix is **result-changing**, because the old ghosts were stale initial values in fresh runs
too, not just after restarts. Fresh flagship run, same deck, 31 cycles, old binary `84a6d248` vs
fixed `967fced6`:

| quantity | relative change |
|---|---|
| total mass | -7.0e-10 |
| magnetic energy | -1.6e-08 |
| rho_max | -5.0e-07 |
| kinetic energy | -1.5e-06 |
| time reached | -3.4e-06 |
| **total radiation energy** | **+5.5e-05** |
| `rad.Fr` (pointwise, at the boundary) | 7–15 % |

So at this epoch the correction to the bulk state is 1e-6 or smaller and the correction to the
total radiation energy is 5.5e-5. **This is a 31-cycle measurement and must not be extrapolated**:
the boundary error is a flux, so it accumulates, and nothing here bounds it at 1e4–1e5 cycles.

---

## Defect 2 (NOT on the flagship path) — the gray coupling wrote an unconverged iterate

On the gray deck (`n_group` unset), `MatterCoupling`'s Newton loop ends with `if (!conv)
lnfail += 1;` and then uses the last iterate. That iterate is not a function of the cell state —
it is wherever the stalled iteration happened to be — so a last-bit input change moves it by a
FINITE amount.

Evidence, all at one step past the restart on `runs/wp13_restart/straight/fhc.in`
(64 blocks, floor EXACTLY zero in every field, so any nonzero divergence is real):

| exp | job | change | grav.phi rel | prim rel | rad.Er rel |
|---|---|---|---|---|---|
| B | 2490927 | production (rkl2) | 6.06e-09 | 4.69e-05 | 8.54e-02 |
| A' | 2490929 | `unsplit`, freeze off | 1.17e-04 | 1.87e-02 | 2.21e-03 |
| E | 2490934 | `chemistry=false` | 6.13e-09 | 5.17e-04 | 1.69e-03 |
| **D** | 2490933 | **`self_gravity=false`** | — | **0** | **0** |
| **F** | 2490935 | restart advances **0** cycles | **0** | **0** | **0** |
| G | 2490936 | gravity tolerance 1e-6 -> 1e-12 | 2.02e-15 | 2.67e-10 | 4.04e-02 |
| **Q** | 2491803 | **`inner_iteration_tol` 1e-8 -> 1e-5** (every cell converges) | 1.97e-15 | 2.42e-10 | **4.45e-12** |

**Q is the proof.** Holding the gravity perturbation at 1e-15 and merely letting the coupling
converge drops the rad.Er restart divergence from 1.7743e-02 to 2.5975e-12 — a factor of 6.8e9 —
and the run emits no non-convergence warnings.

Two supporting measurements:

* At tolerance 1e-12 exactly **one** cell of 262144 moves (by 4 %); the other 236 differing cells
  are at 5.6e-17, i.e. round-off. `max/median = 3.2e14`. That cell is the densest in its block
  (rho = 60.673, the block maximum) and its rho and P agree between the legs to 1e-15. The log
  reports exactly one non-convergent cell.
* The jump AMPLITUDE saturates while the COUNT collapses — 1e-6 / 1e-9 / 1e-12 gravity tolerance
  gives rad.Er max|abs| = 4.99e-02 / 1.7735e-02 / 1.7743e-02 over 261705 / 65567 / 237 cells.
  Fixed-size jumps in fewer cells, not smooth amplification.

### The fix

The gray kernel now falls back to **pure bisection on the monotone gray residual**
(`GrayResidual`, which is `MGResidual` specialised to one group) whenever the Newton reports
non-convergence. 80 fixed iterations, no early exit (warp-divergence-free), guaranteed to
converge, and its answer depends only on `(rho, e0, Er0, dt)`. Cells that already converged are
untouched, so this is bit-identical wherever `conv` is true. This mirrors what
`MatterCouplingMultigroup` has done via rtsafe since the multigroup work.

---

## Falsified along the way — do not re-derive

* **"The gravity solve is the flagship's problem."** It is the dominant term on the *gray* deck
  (exp G: tightening the BiCGSTAB tolerance 1e6x drops `grav.phi` divergence 3e6x, tracking it —
  the tolerance is ABSOLUTE, `relative_residual = false`). On the FLAGSHIP deck `grav.phi`
  restart divergence is 0.88x its own floor, i.e. gravity is not implicated there at all.
* **"Radiation converts 1e-15 into 1e-2 in one step, independent of input size."** The
  observation was right; the earlier reading of it as a bracket flip or table-index flip was
  wrong. It is the unconverged-iterate write (exp Q).
* **Newton iteration count.** `inner_iteration_max` 100 -> 5000 -> 101 -> 103 are all
  BIT-IDENTICAL. The iteration is stalled at a fixed point, and it is not a parity effect.
* **Tabulated-EOS `Cv` discontinuity.** `Eint` is bilinear, so a forward difference of it is
  continuous; the measured jump across table node lines at the flipped cell's density is ~1e-5
  fractional, not O(1). (`runs/wp13b_gpu_restart/diag_cv_discontinuity.py`.)
* **Bell & Lin regime skip.** The default `bell_lin_fix_regime_skip = false` walk really does have
  a 4.88-decade kappa cliff at T = 4123 K at this cell's density (rho = 3.32e-17 g/cm^3, regime 6
  skipped, kappa 4.56e-06 -> 0.348) — but turning the fix on is BIT-IDENTICAL in both tolerance
  legs, so no iterate crosses it. Real, latent, not this bug. No deck in the tree sets it.
* **Radiation substep count.** `nsub = ceil(dt/dt_rad)` is an integer (~384, drifting 383->387
  over 31 cycles) and was a good candidate, but `RAD_PRINT_NSUB=1` shows the sequence is
  IDENTICAL across all three legs at every step (job 2491807).
* **RKL2 STS state, block ordering, chemistry, and the restart FILE.** `unsplit` fails too (it
  only changes which field carries the error); `LogicalLocations`/`Levels` are identical; exp E
  and exp F.

## Impact

* Every ensemble member restarted several times, so each is a valid realisation but not THE
  realisation an uninterrupted run would have produced. This is the mechanism behind the point012
  fresh-vs-restart divergence.
* The matched-epoch mu_core measurement is taken at matched physical state (rho_max), not matched
  time, and the ensemble's spread already carries a measured turbulence-seed noise term (57.6 %
  CoV) of the same character, so realisation-level differences do not move the statistics.
* What is NOT bounded by this work: how the +5.5e-5 radiation-energy correction accumulates over
  the 1e4–1e5 cycles a production member actually runs. That needs a long-baseline A/B, not this
  31-cycle harness.

## Binaries

* `athenaPK_PRESERVED_84a6d248` — the binary the ensemble ran (= source `0d3a559`). Unchanged.
* `athenaPK_PRESERVED_967fced6` — candidate carrying both fixes, built from the working tree on
  top of `5761651`. This is the binary all "FIXED" rows above were measured with.

---

# 2026-08-08 — LONG-BASELINE A/B: does the fix move the science?

**Status: ANSWERED at ~500 cycles / 4.4x rho_crit. The DiodeBC fix is INDISTINGUISHABLE from
run-to-run non-determinism on every conserved integral. This does NOT bound the full 1e4-1e5
cycle production baseline.**

## Why this was needed

The fix's fresh-run effect had only ever been measured at **31 cycles** (Er_tot +5.5e-05, mass
-7.0e-10, KE -1.5e-06). Those numbers are tiny and they are also **useless as a bound**: the
defect is a BOUNDARY FLUX error, so it accumulates every step. Nothing said whether 5.5e-05 at
cycle 31 becomes 5e-05 or 5e-01 at cycle 1e4.

## Design — three legs, and the third is the whole point

`runs/wp13b_ab/submit_ab.sh`, jobs 2492236 / 2492237 / 2492238, all COMPLETED, 2h15m each,
4x H100, ensemble point000 deck, **fresh from t=0 and restart-free** (a restart in the old-binary
leg would confound the fresh-run change with the zeroed-ghost corruption, which is exactly what
this test separates).

| leg | binary | end cycle | rho_max reached |
|---|---|---|---|
| `old_a` | `84a6d248` (what the 24-member ensemble ran) | 514 | 4.44e-13 g/cm3 = **4.44x** rho_crit |
| `old_b` | `84a6d248`, IDENTICAL to old_a | 476 | 3.81e-13 = 3.81x |
| `new` | `967fced6` (the fix) | 489 | 4.33e-13 = 4.33x |

`old_a` vs `old_b` is the **non-determinism floor**; the legs already diverge in cycle count
(514/476/489), which is what makes the floor leg indispensable.

## Result

Full-precision `.hst`, all 515 rows matched between legs to within |dt| < 2e-4 with **no
interpolation**:

| column | scale | FLOOR (a vs b) | SIGNAL (a vs new) | ratio | signal/scale |
|---|---|---|---|---|---|
| mass | 2.053e+03 | 0 | **0 (bit-identical)** | — | 0 |
| KE | 3.354e+03 | 1.544e+01 | 2.305e+01 | **1.49** | 6.9e-03 |
| tot-E | 6.574e+03 | 9.818e+01 | 4.520e+01 | **0.46** | 6.9e-03 |
| ME | 7.167e+01 | 1.875e-01 | 3.420e-01 | **1.82** | 4.8e-03 |
| 1/2/3-mom | — | — | — | 0.91–1.10 | ≤2.1e-04 |

**Every ratio is <= 1.82, and tot-E moves LESS than the floor.** Mass is bit-identical in all
three legs. At this depth two *identical* runs already differ by 0.46 % in KE; the fix adds
nothing beyond that.

For scale: the ensemble's turbulent-seed CoV on mu_core is **57.6 %**. A <=0.7 % shift is ~80x
smaller and cannot move the published distribution.

### The interpolation trap, checked and cleared

The first pass interpolated all three legs onto a common time grid and reported a mass signal of
3.9e-06 against a floor of exactly 0. That was an **interpolation artifact**: recomputed by
nearest-row matching with no interpolation, the mass difference is exactly 0. Every other ratio is
unchanged to two digits (KE 1.48 -> 1.49, tot-E 0.46 -> 0.46, ME 1.82 -> 1.82), so the conclusion
does not rest on the interpolation scheme -- but the one number that looked anomalous was entirely
produced by it.

## What this does and does not license

**Licensed:** the 24-member ensemble's published `mu_core` statistics do not need to be
re-measured on account of the DiodeBC defect. The effect is below the floor at 16x the previously
tested baseline, and ~80x below the seed scatter that already dominates those numbers.

**NOT licensed:** this reaches cycle ~500 and 4.4x rho_crit. A production member runs 1e4-1e5
cycles and to 1e5 x rho_crit. The accumulation question is *bounded much better than before* and
is *not closed*. The legs carry no restart files by construction, so this baseline cannot be
extended -- a deeper answer needs a new, longer pair.

---

# Long-baseline A/B: the DiodeBC fix does not move the science to 55-65x rhocrit (2026-08-09)

## The question this closes

The shallow A/B (jobs 2492236/7/8) found the fix indistinguishable from the non-determinism floor
at 476-514 cycles / 4.4x rhocrit. That was explicitly **not** a bound: the defect is a BOUNDARY
FLUX error, so it accumulates, and a production member runs 1e4-1e5 cycles. Nothing said whether
5.5e-05 at cycle 31 becomes 5e-05 or 5e-01 later.

## Result — jobs 2492767 (old_a), 2492768 (old_b), 2493189 (new)

All three ran their full internal `-t 05:10:00` and ended COMPLETED 0:0, reaching **1031 / 1023 /
1026 cycles** and **54.66x / 55.41x / 64.91x rhocrit** — a factor ~12 deeper in density than the
shallow test and past every epoch the ensemble measures at.

| column | scale | FLOOR (a-b) | SIGNAL (a-new) | ratio | tight-subset ratio | verdict |
|---|---|---|---|---|---|---|
| mass | 2.05267e+03 | 0.00000e+00 | 0.00000e+00 | 1.00 | 1.00 | within floor |
| KE | 3.46851e+03 | 5.21000e+00 | 8.99000e+00 | 1.73 | 1.73 | within floor |
| tot-E | 6.63224e+03 | 2.66400e+01 | 5.72400e+01 | 2.15 | 3.10 | within floor |
| ME | 7.30853e+01 | 7.72000e-02 | 1.69600e-01 | 2.20 | 2.08 | within floor |
| 1-mom | 3.92386e+01 | 6.40000e-03 | 8.40000e-03 | 1.31 | 1.30 | within floor |
| 2-mom | 1.79991e+01 | 5.00000e-03 | 1.03000e-02 | 2.06 | 1.66 | within floor |
| 3-mom | 1.09241e+02 | 3.00000e-03 | 2.00000e-03 | 0.67 | 0.67 | within floor |

**No quantity moves as much as 4x the run-to-run non-determinism floor**, against a flag threshold
of 10x, and `mass` is bit-identical in all three legs. The largest shift is 0.86 % of tot-E
against a 0.40 % floor — ~67x below the ensemble's 57.6 % seed CoV. **The published ensemble
mu_core statistics do not need re-measuring**, and the caveat that they were produced with the
defective binary `84a6d248` is now quantified rather than merely acknowledged.

The old_a-vs-old_b floor leg is what makes this readable: 4-rank GPU reductions are not
order-deterministic, and ~1000 cycles of a collapsing turbulent flow amplify a last-bit
difference. "The fix does not move the science" only means "old-vs-new sits inside old-vs-old".

## Two analysis corrections made while producing this

1. **Interpolation fabricates signal.** The original script interpolated all three legs onto a
   synthetic time grid, which produced a `mass` signal of 3.9e-06 against a floor of exactly 0 —
   on a quantity the code conserves bit-for-bit. Each leg gets a different blend of its own
   neighbouring rows, so the difference is built by the analysis. Nearest-row matching removed it
   (mass exactly 0.0) and moved no other ratio by more than two digits. This is now in
   `analyze_ab.py` rather than being done ad hoc alongside it.

2. **The match-quality guard cried wolf.** Judging the match offset against the *median* `.hst`
   cadence flagged this run `MATCH NOT SAFE` at "650 % of one row". The spacing here spans 0 to
   6.5e-3 — sparse early, one row per cycle once the collapse drives dt down — so the median
   (2.0e-5, set by the dense late rows) is not a valid denominator for an offset occurring early.
   Judged against *local* spacing the worst case is 2.00x, and the independent falsifier settles
   it: restricting to the 661 of 999 rows matched within 10 % of local spacing leaves every ratio
   essentially unchanged (column above). A false alarm on a sound result is its own failure mode,
   so the guard now uses local spacing and always prints the tight-subset cross-check.

## What this does NOT claim

The legs stop at ~1030 cycles. Production runs 1e4-1e5. This bounds accumulation over the epoch
where the ensemble is *measured* (1e-12 to 2e-12 g/cm3, passed here by a factor of ~30), not over
a full second-core run.
