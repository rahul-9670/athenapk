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
