# WP-5 — Gravity-solver convergence reporting

**Status: instrument IMPLEMENTED + gated; invariance sweep PASSED at smoke scale.
Deep-leg measurement still open (GPU-blocked). 2026-07-31.**
Source pair: `docs/provenance/binary_5ebddce0/` + the WP-5 diff (uncommitted, see below).

## The problem, confirmed at source

`VALIDATION_PLAN.md` flagged that across 1975 cycles of `runs/convergence_ladder/nj8/run.log`
there are zero occurrences of "residual", "converge" or an iteration count. Confirmed:
`grep -ic` returns **0, 0, 0** over 7080 lines.

The reason is structural, in `external/parthenon/src/solvers/bicgstab_solver.hpp` (~line 397):

```cpp
if ((rms_res < rel_tol) || (rms_res < *abs_res_tol) ||
    (solver->iter_counter >= max_iter)) {
  solver->final_residual  = rms_res;
  solver->final_iteration = solver->iter_counter;
  return TaskStatus::complete;          // <-- IDENTICAL return in all three cases
}
```

**Converging and exhausting `max_iterations` produce the same `TaskStatus::complete`, and
neither is reported.** A solve that silently bails at the ceiling is indistinguishable from a
converged one to every caller. `SolverBase` does expose `GetFinalResidual()` /
`GetFinalIterations()`, but nothing in AthenaPK was reading them.

### Second finding — the production tolerance is ABSOLUTE, not relative

From `bicgstab_solver.hpp:49-66`: the deck does **not** set `relative_residual`, so it defaults
to `false` and the else-branch assigns

```
absolute_residual_tolerance = residual_tolerance = 1.0e-6
relative_residual_tolerance = 0.0
```

So the production criterion is an **absolute** tolerance on
`rms_res = sqrt( Σr² / pmesh->GetTotalCells() )` for ∇²φ = 4πG(ρ − ρ̄), with `four_pi_G = 1`.

This matters: the RHS grows by many decades as the core collapses while the tolerance stays
pinned at 1e-6, so the *relative* accuracy demanded of the solve tightens monotonically through
the run. That is precisely the regime where a fixed iteration ceiling starts to bite — and it is
a testable prediction, not a worry (see the measurement below).

## What was implemented

New files `src/diagnostics/grav_diag.{hpp,cpp}` (house style of `mag_diag.*`), wired into
`src/self_gravity/self_gravity.cpp` (registration) and
`src/self_gravity/self_gravity_driver.cpp` (a `once_per_region` post-solve task), plus one line
in `src/CMakeLists.txt`.

Gate key: **`<self_gravity> solver_diag`, default `false`.** When off, no Params are added, no
history columns are registered and no task is inserted.

Three new history columns (`UserHistoryOperation::max` — the scalars come from a global
all-reduce inside the solver and are identical on every partition, so `sum` would multiply by
the partition count):

| column | meaning |
|---|---|
| `grav-iters` | BiCGSTAB iterations used by the most recent solve |
| `grav-res` | final rms residual of the most recent solve |
| `grav-nonconv` | 1.0 if that solve stopped because it hit `max_iterations`, else 0.0 |

Because the solver cannot distinguish the two exit paths from outside, `grav-nonconv` is
reconstructed from the only observable that differs: `iters >= max_iterations`. Non-convergence
additionally emits a rank-0 warning, rate-limited geometrically (occurrences 1–4, 10, 100, 1000,
then every 10000) so a persistently non-converging run neither hides nor floods the log.

## Gate — OFF-state bit-identical: **PASS**

Method per `VALIDATION_PLAN.md` §2 landmine 9: an **explicitly built control binary**, matched
thread count (`OMP_NUM_THREADS=1`), never a number from an old log.

| binary | md5 | build |
|---|---|---|
| control (WP-5 changes reverted) | `ff6888eef09874ce23a88c44e6d98452` | `build_cpu`, GCC 13.3, Kokkos OpenMP |
| WP-5 | `ffb68c0ccebf47938d3ef606764daf79` | same tree + WP-5 diff |

Incidental: rebuilding the WP-5 binary after the control round-trip reproduced md5
`ffb68c0c…` exactly — **the CPU build is deterministic**, which is a useful data point for the
outstanding WP-0 rebuild-reproduction item.

Case: `runs/eos_smoke/fhc.in` — 12 cycles, 32³ root, self-gravity + EOS table + chemistry + RT
+ non-ideal MHD. Runs live in `runs/wp5_gate/{ctrl,off,on,tol1e8,maxit400}/`.

| comparison | `.hst` | snapshot datasets |
|---|---|---|
| ctrl vs **off** (`solver_diag` absent ⇒ default false) | **byte-identical** | **15/15 bit-identical** in both `00000` and `final` |
| ctrl vs **on** (`solver_diag=true`), physics columns 1–13 | **byte-identical** | — |

So the instrument is bit-identical both when off *and* when on — expected, since the added task
only reads solver state and writes package Params, never field data.

Dataset comparison used `docs/validation/scripts/h5_bitcmp.py`, which compares raw bits via
`.tobytes()` rather than `==` (so ±0.0 and NaN payloads are not glossed over). Byte-comparing
`.phdf` files directly is meaningless — parthenon embeds wall-clock metadata, so two
bit-identical physics runs always differ as bytes.

## First measurement — the solver is healthy at this scale

`runs/wp5_gate/on/parthenon.out0.hst`, `max_iterations = 200`, `residual_tolerance = 1e-6`:

| t | grav-iters | grav-res | grav-nonconv |
|---|---|---|---|
| 0.00000 | −1 | −1 | 0 |
| 0.13282 | 6 | 3.921e−08 | 0 |
| 0.25469 | 5 | 9.824e−08 | 0 |
| 0.37998 | 4 | 6.813e−07 | 0 |
| 0.49426 | 5 | 1.544e−07 | 0 |
| 0.59664 | 5 | 1.748e−07 | 0 |
| 0.68977 | 5 | 1.882e−07 | 0 |
| 0.77540 | 5 | 2.187e−07 | 0 |
| 0.85474 | 5 | 2.538e−07 | 0 |
| 0.92873 | 5 | 5.253e−07 | 0 |
| 0.99751 | 5 | 8.954e−07 | 0 |
| 1.06148 | 5 | 9.508e−07 | 0 |
| 1.12114 | 6 | 1.438e−07 | 0 |

(The t = 0 row is the `-1` "no solve yet" sentinel — history is written before the first solve.)

**The solve converges every cycle, using 4–6 of 200 available iterations.** No non-convergence
warning fired. At this scale there is no problem.

**But the predicted trend is already visible.** `grav-res` climbs monotonically
3.92e−08 → 9.51e−07 between t = 0.13 and t = 1.06 — a **24× rise, reaching 95 % of the fixed
1e−6 absolute tolerance** — while the iteration count stays flat at 5. At t = 1.12 the solver
spends a 6th iteration and drops back to 1.44e−07. That is the absolute-tolerance mechanism
above, behaving exactly as predicted and self-correcting by spending more iterations.

This is a 12-cycle smoke test in which the density has barely grown. In the production collapse
ρ_max rises through ten-plus decades, so the iteration count needed to hold a *fixed absolute*
residual must keep climbing. **Whether it climbs to 200 is now measurable, and was not before.**

## Invariance sweep — **PASS** at this scale

### `max_iterations` 200 → 400: bit-identical

Predicted before running (the ceiling is never reached, so raising it can change nothing) and
then checked, per working rule 3 — construct the observation that would hold if the diagnosis
were wrong. It held:

| comparison | `.hst` physics cols 1–13 | final snapshot |
|---|---|---|
| `on` (200) vs `maxit400` (400) | **byte-identical** | **15/15 datasets bit-identical** |

Iteration counts are identical row-for-row (6, 5, 4, …), residuals identical to all printed
digits. **`max_iterations = 200` is not a binding constraint at this scale**, and the instrument
correctly reports so.

### `residual_tolerance` 1e-6 → 1e-8: differences at output precision

| t | iters @1e-6 | res @1e-6 | iters @1e-8 | res @1e-8 |
|---|---|---|---|---|
| 0.13282 | 6 | 3.921e−08 | 7 | 4.748e−09 |
| 0.49426 | 5 | 1.544e−07 | 7 | 8.756e−09 |
| 0.99751 | 5 | 8.954e−07 | 8 | 4.518e−09 |
| 1.06148 | 5 | 9.508e−07 | 9 | 2.329e−10 |
| 1.12114 | 6 | 1.438e−07 | 8 | 2.305e−09 |

A 100× tighter tolerance costs **+2 to +4 iterations** (5–6 → 7–9), still an order of magnitude
below the 200 ceiling, and drives the residual to ~1e−9.

Resulting state difference after 12 cycles, on the **final** snapshot:

| dataset | max abs. diff | max abs. / RMS(field) | in float32 ulps |
|---|---|---|---|
| `grav.phi` | 1.526e−05 | 5.94e−07 | **2.00** |
| `prim` | 1.907e−05 | 8.53e−05 | **5.00** |
| `rad.Er` | 5.960e−08 | 1.55e−07 | **2.00** |

Snapshot output is **float32**, and the differences are exactly 2, 5 and 2 units in the last
place of it — at, or barely above, the resolution of the output itself.

**Two measurement traps avoided here, both worth recording:**

- The `.hst` file prints only ~6 significant figures (`-8.16822e+00`). A naive column diff
  reported "differences" of exactly 1.0e−5, 1.0e−4, 1.0e−6 — suspiciously round powers of ten,
  because they were 1 ulp *of the printed text*, not of the state. Never conclude a physics
  difference from `.hst` text at this level; go to the snapshots.
- Comparing the `00000` snapshot proves nothing — it is the initial condition, written before
  the first solve, and is trivially identical. Only the `final` snapshot is a test.
- A pointwise *relative* error on `prim` reported `max|rel| = 5.35`, which looks catastrophic
  but is a near-zero denominator where a velocity or field component crosses zero. Normalize by
  the field RMS (or count ulps), not by the local value.

**Verdict: at eos_smoke scale, `residual_tolerance = 1e-6` and `max_iterations = 200` are both
adequate**, and the production settings are not limiting the answer. This is far below any
plausible WP-18 seed scatter.

**Caveat, stated plainly:** 12 cycles at 32³ is a very short lever arm and the density has barely
grown. This result does **not** transfer to the deep collapse — which is exactly what item 2
below is for.

## Remaining work

1. **Blocked on GPU nodes:** the measurement that matters is `grav-nonconv` on a *deep* leg
   carried to ρ = 1e-13 and beyond, where the RHS is decades larger. That needs a GPU rebuild,
   which must NOT go into `build_gpu` while the ladder jobs re-exec it at slot boundaries
   (same constraint as the WP-0 rebuild check — do both in one `build_gpu_repro` pass).
2. Re-run the tolerance sweep on that deep leg. The smoke-scale PASS above does not transfer.
3. Consider whether production should switch to `relative_residual = true`. That is a
   **result-changing** deck change, so it is the user's call, not a numerical fix — but the
   trend above is the argument for it.

## Incidental: running AthenaPK on the front-end

Not documented anywhere and cost several failed attempts. A serial front-end run needs **all**
of:

```bash
export OMPI_MCA_pml=ob1        # else MPI_Init dies in psm2_ep_open (Bus error).
                               # OMPI_MCA_mtl=^psm2 from athenapk_env.sh is NOT enough:
                               # the *ofi* MTL loads the psm2 provider via libfabric anyway.
export OMPI_MCA_io=romio341    # else every MPI_File_open fails ("could not be opened"),
                               # including reading the input deck. This is in every submit
                               # script but NOT in ~/athenapk_env.sh.
export FI_PROVIDER=tcp
env -C <rundir> <binary> -i <rundir>/deck.in    # `cd` in a compound shell command does not
                                                # reliably reach the process in this harness
```

## Confidence

*Verified* — the source mechanism (identical return for converged vs. bailed), the tolerance
semantics (absolute, not relative), the OFF-state bit-identical gate against an explicitly
built control, and the measured residual trend are all from artifacts produced this session.
*Inferred* — that the iteration count will climb toward the ceiling in the deep collapse. That
is the prediction item 2 exists to test; it is **not** yet observed.
