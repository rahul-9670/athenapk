# WP-13b — restart reproducibility in the PRODUCTION configuration (GPU, 4 ranks, AMR, RKL2)

**Status: FAIL, root cause localized to TWO independent mechanisms. The restart FILE is perfect;
the defect is entirely in the first post-restart step. (1) The self-gravity BiCGSTAB solve is
determined only to its residual tolerance, so a restart lands on a different iterate inside the
tolerance ball — this scales exactly with the tolerance and is a POLICY issue, not a bug.
(2) The radiation solve converts that perturbation into an O(1e-2) relative change in `rad.Er` in
a SINGLE step, and does so essentially independently of the perturbation's size — the signature of
a discrete branch selection, not smooth amplification. (2) is the real defect.**

## Why this exists

WP-13 is closed and its fix is real, but it was validated on `build_cpu`, one rank,
`diffusion/integrator = unsplit`, 12 cycles. Production is GPU + 4 MPI ranks + RKL2 STS. The
trigger was ensemble member point012 (2026-08-08): a fresh leg and a leg restarted from the t=0
restart file reached cycle 200 with `t = 1.0262248` vs `0.9922842` and `dt = 1.940e-4` vs
`9.970e-5`. That is not by itself evidence of a bug — GPU reductions need not be
order-deterministic — so the test below measures the non-determinism floor explicitly.

## Design

Three legs, same deck (`runs/wp13_restart/straight/fhc.in`), same binary
(`athenaPK_PRESERVED_84a6d248`, the one the ensemble ran), same rank count:

| leg | what it runs |
|---|---|
| `fresh_a` | `-i fhc.in`, nlim=N |
| `fresh_b` | `-i fhc.in`, nlim=N — IDENTICAL to `fresh_a` |
| `split` | `-i fhc.in` nlim=N/2, then `-r <rst>` nlim=N |

`fresh_a` vs `fresh_b` is the **non-determinism floor**; `fresh_a` vs `split` is the restart
divergence. Without leg B the test cannot tell a restart defect from ordinary non-determinism.
Harness: `runs/wp13b_gpu_restart/submit_wp13b.sh` (parameterized by `WP13B_DIR/DIFF/NLIM/SPLIT/
FREEZE/EXTRA`) and `compare_wp13b.py`.

**The floor is EXACTLY ZERO in every experiment below.** The code is fully deterministic in this
configuration, so every nonzero restart divergence is real.

### Two traps this test walked into first

* **The first run was VACUOUS.** The deck's `tlim = 1.5` is reached at cycle 22, before the
  nlim=30 split, so the "restart" began already at `tlim` and advanced ZERO cycles in 0.108 s.
  All three legs compared bit-identical for a trivial reason. Fixed with `tlim=1000` plus positive
  checks that every leg ends at `nlim` and the restart leg STARTS at the split cycle.
* `rkl2_freeze_eta=true` hard-fails unless `integrator=rkl2`, so the `unsplit` leg must also pass
  `WP13B_FREEZE=false` (job 2490926 aborted 134 on all three legs before this was understood).
* `physics/radiation=false` is not a usable single-variable experiment: `eos=hydrogen requires
  <physics> radiation=true` (`hydro.cpp:1001`).

## Results — one variable per experiment, all at one step past the restart

| exp | job | change | `grav.phi` rel | `prim` rel | `rad.Er` rel | verdict |
|---|---|---|---|---|---|---|
| B | 2490927 | production (rkl2) | 6.06e-09 | 4.69e-05 | **8.54e-02** | FAIL |
| A' | 2490929 | `unsplit`, freeze off | 1.17e-04 | 1.87e-02 | 2.21e-03 | FAIL |
| E | 2490934 | `chemistry=false` | 6.13e-09 | 5.17e-04 | 1.69e-03 | FAIL |
| **D** | 2490933 | **`self_gravity=false`** | — | **0** | **0** | **BIT-IDENTICAL** |
| **F** | 2490935 | restart advances **0** cycles | **0** | **0** | **0** | **BIT-IDENTICAL** |
| G | 2490936 | `residual_tolerance` 1e-6 → 1e-12 | 2.02e-15 | 2.67e-10 | 4.04e-02 | FAIL |

At nlim=60 (job 2490880) the production configuration reaches `prim` max|abs| = 2.704266,
`t = 2.06492729943642361` vs `2.06493272048885412`.

### What each row establishes

* **D — self-gravity is necessary for the failure.** With it off the restart is bit-identical.
* **F — the restart file is NOT the problem.** Load it, write it straight back, and every field
  round-trips bit-exactly on GPU, reproducing WP-13's CPU result. The defect is in the first step.
* **A' — RKL2 is exonerated.** `unsplit` fails too; the integrator changes which field carries the
  error, not whether there is one.
* **E — chemistry is exonerated.**
* **G — the gravity mechanism is CONFIRMED and quantified.** Tightening the tolerance by 1e6 drops
  `grav.phi` divergence by ~3e6 (6.06e-09 → 2.02e-15) and `prim` by ~2e4. The Poisson solution is
  determined only to the solver's residual tolerance, and the restart lands on a different iterate
  inside that ball. `solver = BiCGSTAB`, `residual_tolerance = 1.0e-6`, `max_iterations = 200`,
  and CLAUDE.md notes the tolerance is ABSOLUTE (`relative_residual = false`).

### The second mechanism: radiation saturates

Across the 1e6 tolerance change, `rad.Er` moved only 8.54e-02 → 4.04e-02 — a factor of **2.1**,
while the perturbation feeding it shrank by a factor of **3e6**:

| | tol 1e-6 | tol 1e-12 | ratio |
|---|---|---|---|
| `grav.phi` rel | 6.06e-09 | 2.02e-15 | 3.0e+06 |
| `rad.Er` rel | 8.54e-02 | 4.04e-02 | 2.1 |

An O(1e-15) relative perturbation in `phi` produces an O(1e-2) relative change in `rad.Er` **in a
single step**. One step cannot amplify chaotically, and the response is nearly INDEPENDENT of the
input size, which is the signature of a **discrete branch selection** — a bracket flip, a table
cell index change, or an iteration-count change — rather than smooth ill-conditioning. Exp D is
consistent: with gravity off there is no perturbation at all, and `rad.Er` is bit-identical.

Prime suspects, in the order I would instrument them: the safeguarded Newton+bisection in the
multigroup matter coupling (CLAUDE.md records plain Newton overshooting to NaN on the tabulated-EOS
H2 kinks, which is exactly the kind of surface where an infinitesimal input change selects a
different bracket), the opacity-table lookup index, and the EOS table cell selection.

## Impact

* **Not a physics error.** Each leg is a valid solution of the same IC; nothing here says either
  answer is wrong.
* **It does invalidate the claim "restart is reproducible"** for the production configuration.
  Every ensemble member restarted several times, so each is a valid realization but not THE
  realization a single uninterrupted run would have produced. This is the mechanism behind the
  point012 fresh-vs-restart divergence.
* The matched-epoch μ_core measurement is unaffected: it is taken at matched physical state
  (ρ_max), not matched time, and the ensemble's spread already carries a measured
  turbulence-seed noise term of the same character.

## Not yet done

Instrumenting the radiation matter-coupling to identify which discrete decision flips. That is a
code change, not another run of this harness.
