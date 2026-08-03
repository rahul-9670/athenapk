# B3 — is the self-gravity solver's ABSOLUTE residual tolerance safe at production depth?

**Status: RESOLVED. The feared failure mode does not occur; the recommendation is still to switch
to the relative criterion, for a different and better reason.**

- **The absolute criterion never fails to converge**, even with the Poisson source scaled over
  **10 decades**. 1–2 BiCGSTAB iterations at every scale, `grav-nonconv = 0` throughout.
- **But its accuracy is scale-dependent by accident**: fractional residual `res/rms(rhs)` drifts
  from 4.4e-11 at ρ = 1 to 3.1e-19 at ρ = 1e10 — it demands whatever the scale happens to imply.
- **The relative criterion holds a constant fractional accuracy of 4.4e-14 across all 10 decades**,
  to five significant figures. That is the physically meaningful demand, since φ and ∇φ scale with
  the source.
- Cost of switching: 3 iterations instead of 1–2 here, 7.75 instead of 4.92 on the production smoke
  deck. Result shift: **1.8e-5** relative — negligible against σ ≈ 16 %.

---

## 1. The original concern, and my first analysis of it — which was wrong

B3 was raised as: the deck omits `relative_residual`, so `1e-6` is an **absolute** tolerance. As the
collapse deepens, `rms(rhs) = rms(4πGρ)` grows by decades while the ceiling does not, so the
effective *relative* demand tightens monotonically — and the solve should eventually be unable to
meet it, exhaust `max_iterations = 200`, and exit at an uncontrolled accuracy with (before the B2
fix) no signal at all.

My first pass argued the **opposite** — that `relative_residual = true` should be *rejected*,
because `rel_tol = 1e-6 × rms(rhs)` gets looser in absolute terms as the source grows. That
reasoning weighed the absolute error without weighing the growing physical scale it is measured
against. A fixed *fractional* accuracy in φ is the meaningful demand; a fixed *absolute* one is
scale-dependent, and which way it errs depends entirely on the units the problem happens to be in.
Both framings were half-right and neither settled anything, which is why this was measured.

## 2. What the code actually does

`external/parthenon/src/solvers/bicgstab_solver.hpp`:

```cpp
// :56-66  — the two criteria are MUTUALLY EXCLUSIVE, not additive
if (relative_residual) {
  *relative_residual_tolerance = GetOrAddReal(block, "relative_residual_tolerance", *residual_tolerance);
  *absolute_residual_tolerance = GetOrAddReal(block, "absolute_residual_tolerance", 0.0);
} else {
  *relative_residual_tolerance = GetOrAddReal(block, "relative_residual_tolerance", 0.0);
  *absolute_residual_tolerance = GetOrAddReal(block, "absolute_residual_tolerance", *residual_tolerance);
}
// :405-411
Real rel_tol = *rel_res_tol * std::sqrt(solver->rhs2 / pmesh->GetTotalCells());
const bool tol_met = (rms_res < rel_tol) || (rms_res < *abs_res_tol);
```

Whichever criterion is not selected is driven to **0** so it can never fire. There is **no way to
demand that both hold** — the accept test is an OR, so if both are set the *looser* one governs.
Note also that `rms(rhs)` is a plain cell-count average over `GetTotalCells()`, so under Jeans AMR
it is dominated by core cells.

## 3. The measurement

Four attempts to reach the deep regime by actually collapsing gas all failed, for four different
reasons, and are recorded because each one cost a job:

| attempt | configuration | outcome |
|---|---|---|
| 1 | full physics, numlevel=5 | 290 s/step; cycle 12 in 47 min — could not reach depth in budget |
| 2 | adiabatic γ=5/3 | pressure-supported core; `nblk` froze at 288 by t = 0.98 |
| 3 | tabulated EOS, radiation off | refused at init — `eos=hydrogen` **requires** `radiation=true` |
| 4 | adiabatic γ=1.1 | `nblk` **also** froze at 288; `grav-iters` plateaued at 8/200 |

The mistake was insisting on a *physically driven* source. **The Poisson operator is linear**:
ρ → λρ scales φ and the residual by exactly λ. So the question is answerable with a *static* source
of controlled amplitude, using `src/pgen/poisson_test.cpp` — a uniform sphere, multipole BCs (as in
production), `nlim = 1`, `rho_in` as a free knob. Seconds per leg.
(`runs/switch_probe/submit_scale.sh`, job 2448617.)

### One invalid run, discarded

Job 2448614 passed `relative_residual_tolerance` on the command line in **every** leg. Since
`BiCGSTABParams` reads that key with `GetOrAddReal`, supplying it kept the relative test live even
in the legs labelled "absolute", and the OR meant it won. Two tells: `abs6` and `rel6` came out
identical to every printed digit, and `abs6` "converged" at `rho_in = 1e8` with
`grav-res = 4.43e-6`, four times **above** its own ceiling. The deck (`sphere_multipole.in`) is a
second source of the same contamination — it sets all three keys — so both tolerances must be
pinned explicitly, with the unselected one driven to 0.

### Result

| ρ_in | **absolute 1e-6** iters / res / res·ρ⁻¹ | **absolute 1e-8** | **relative 1e-6** |
|---|---|---|---|
| 1 | 1 / 4.36e-11 / 4.36e-11 | 1 / 4.36e-11 / 4.36e-11 | 3 / 4.43e-14 / **4.4257e-14** |
| 1e2 | 1 / 4.36e-09 / 4.36e-11 | 1 / 7.33e-12 / 7.33e-14 | 3 / 4.43e-12 / **4.4257e-14** |
| 1e4 | 1 / 7.33e-10 / 7.33e-14 | 1 / 5.52e-12 / 5.52e-16 | 3 / 4.43e-10 / **4.4257e-14** |
| 1e6 | 1 / 5.54e-10 / 5.54e-16 | 1 / 1.45e-10 / 1.45e-16 | 3 / 4.43e-08 / **4.4257e-14** |
| 1e8 | 1 / 1.44e-08 / 1.44e-16 | 2 / 3.39e-11 / 3.39e-19 | 3 / 4.43e-06 / **4.4257e-14** |
| 1e10 | 2 / 3.13e-09 / 3.13e-19 | 2 / 3.13e-09 / 3.13e-19 | 3 / 4.43e-04 / **4.4257e-14** |

`grav-nonconv = 0` in all 18 legs.

**Two things this shows.**

1. **The feared failure does not happen.** The absolute criterion converges in 1–2 iterations over
   10 decades of source. The reason is that multigrid-preconditioned BiCGSTAB converges
   superlinearly on this problem and *overshoots* whatever threshold it is given — at ρ_in = 1e10 it
   asked for 1e-6 and delivered 3.13e-9. The exact stopping rule is barely load-bearing here.
2. **The relative criterion is exactly scale-invariant**, as the linearity of the operator predicts:
   `res/ρ_in = 4.4257e-14` to five significant figures across ten decades. The absolute criterion's
   fractional accuracy, by contrast, wanders over eight decades (4.4e-11 → 3.1e-19) — not because
   anything is wrong, but because a fixed absolute number means different things at different
   scales.

## 4. Scope — what this does not prove

This is a **smooth uniform sphere on a uniform grid**, where multigrid is near-ideal (1–3
iterations). Production is a deeply refined AMR hierarchy with a near-singular profile, where the
preconditioner is much less effective: the depth probes reached t = 1.17–1.30 with `grav-iters`
plateauing at **8 of 200**, `grav-nonconv = 0`, and `grav-res` peaking at **9.98e-7 — 99.8 % of the
1e-6 absolute ceiling**. So the criterion is scale-safe, but the *margin* under the absolute rule at
production depth is thin, and it is thin for exactly the reason B3 identified.

Production runs (`prod_v9`, `prod_t4_full`) reached ρ = 3.6e5 ρ_crit with no reported
non-convergence, but those binaries predate the B2 unconditional-warning fix, so absence of a
warning there is weak evidence and is not being counted.

## 5. Recommendation

**Set `relative_residual = true`** (leaving `relative_residual_tolerance` at its default, which
equals `residual_tolerance = 1e-6`).

Not because the absolute rule is failing — it is not — but because:

1. It holds a **constant fractional accuracy in φ** at every depth (measured: 4.4e-14 over 10
   decades). φ's gradient drives the collapse, and a fractional demand is the one that means the
   same thing at every density.
2. The absolute rule's accuracy is an accident of units, and its measured margin at production depth
   is 99.8 % of tolerance — no headroom for a deeper run than has already been done.
3. The cost is small and bounded: 3 iterations vs 1–2 here, 7.75 vs 4.92 on the production smoke
   deck (+52 %), against `grav-iters` of only 4–8 out of `max_iterations = 200`. The solver is not
   the bottleneck (B8's cost hypothesis was independently falsified).

**It is result-changing**, at the 1.8e-5 relative level measured on the production smoke deck —
four decades below σ ≈ 16 %, and folded into the single authorised re-baseline together with the
diode BC.

> *Confidence:* the scaling table, the iteration counts and `grav-nonconv = 0` — **measured** this
> session on `build_cpu`. The mutual exclusivity of the two criteria — **verified** from
> `bicgstab_solver.hpp:56-66` and `:405-411`. The production-depth margin (99.8 % of tolerance) —
> **measured**, but on a 32³ AMR probe, not at production resolution.
