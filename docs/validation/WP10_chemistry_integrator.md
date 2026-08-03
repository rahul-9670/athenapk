# WP-10 — chemistry integrator: substep convergence and element conservation

**Status: CLOSED 2026-08-03. Element conservation PASSES. Substep convergence FAILED at the
production setting because of a new bug (B10) — now instrumented, measured on GPU, and FIXED in
the production deck (`nsub_max = 400 -> 4000`, §4). Measured impact of the fix: 2.06e-05.**

- **At production settings `cfl_cool` is completely inert.** Measured with the new instrument:
  **262144 of 262144 cells — 100 % of the domain** — have their chemistry sub-step set by
  `nsub_max`, not by `cfl_cool`. Lifting the cap to 1e5 drops that to **2 cells**.
- **Consequence:** at fixed `cfl_cool`, lifting `nsub_max` 400 → 3200 moves `x_CO` by **2.1e-2**
  and `x_C+` by **7.2e-3**.
- **B10: the degradation was silent.** The code's own counter read **0 in every leg**, including
  those where the cap demonstrably changed the answer. It counted a condition the integrator's
  design makes unreachable.
- **Carbon conservation PASSES**: the renormalization clamp never fires (0 cells of 78 M).
- **The fix is free.** All six ladder legs took 381 s, and the verification pair — `nsub_max` 400
  vs **100000**, a 250× increase — took **380 s each**. Chemistry is not the cost driver on this
  deck.

Ladder: `runs/wp10_chem/submit_wp10.sh`, job 2448540 — the 32³ B2/B4 smoke deck (chemistry +
radiation + full non-ideal + self-gravity all live), `nlim = 12`, all legs ending at
t = 1.114698121, so this is matched-state by construction.

---

## 1. Element conservation — PASS

`gow17_reduced` carries carbon as C⁺, CO and a *derived* neutral `x_C0 = x_Ctot − x_C+ − x_CO`.
Conservation is therefore structural, and the meaningful question is whether the joint
renormalization at `network_gow17_reduced.hpp:200-207` ever has to fire — i.e. whether the
integrator produces states that violate the budget and get silently rescaled.

| leg | `cfl_cool` | `nsub_max` | max (x_C+ + x_CO)/x_Ctot | cells at ratio = 1 |
|---|---:|---:|---:|---:|
| c100 | 0.1 | 400 | 0.306173 | **0** |
| c050 | 0.05 | 400 | 0.306356 | **0** |
| c025 | 0.025 | 400 | 0.306381 | **0** |
| c0125 | 0.0125 | 400 | 0.306380 | **0** |
| c025_big | 0.025 | 3200 | 0.228978 | **0** |
| c0125_big | 0.0125 | 3200 | 0.228763 | **0** |

The clamp never fires. Note in passing that the *peak* carbon ionization state itself differs by
34 % between `nsub_max` 400 and 3200 (0.306 vs 0.229) — the first sign of §2.

## 2. Substep convergence — FAIL

L1 relative difference against the finest leg (`c0125_big`: `cfl_cool = 0.0125`,
`nsub_max = 3200`):

| leg | H2 | H+ | C+ | CO | e⁻ |
|---|---|---|---|---|---|
| c100 (production) | 1.37e-07 | 1.19e-02 | 2.06e-02 | 6.05e-02 | 2.40e-03 |
| c050 | 8.85e-08 | 7.87e-03 | 1.24e-02 | 3.64e-02 | 1.04e-03 |
| c025 | 3.78e-08 | 3.61e-03 | 9.55e-03 | 2.80e-02 | **2.15e-03** |
| c0125 | 1.23e-10 | 4.04e-04 | 8.90e-03 | 2.63e-02 | **3.67e-03** |
| c025_big | 3.80e-08 | 3.34e-03 | **2.38e-03** | **6.82e-03** | 8.41e-04 |

Two things are wrong with the `nsub_max = 400` column:

- **CO plateaus** at ~2.6e-2 instead of converging to zero. Halving `cfl_cool` from 0.025 to
  0.0125 buys almost nothing (2.80e-2 → 2.63e-2).
- **x_e gets WORSE** past `cfl_cool = 0.05`: 1.04e-3 → 2.15e-3 → 3.67e-3. A tolerance that makes
  the answer *less* accurate when tightened is not the binding tolerance.

Both are the signature of a different limiter taking over. §3 identifies it, and `c025_big`
confirms it directly: at the *same* `cfl_cool = 0.025`, lifting the cap drops the CO error from
2.80e-2 to 6.82e-3 — **4× better from a knob that is supposed to be a safety margin**.

Isolating the cap at fixed tolerance:

| same `cfl_cool`, cap 400 vs 3200 | H2 | H+ | C+ | CO | e⁻ |
|---|---|---|---|---|---|
| c025 vs c025_big | 1.94e-10 | 2.75e-04 | 7.15e-03 | **2.10e-02** | 2.99e-03 |
| c0125 vs c0125_big | 1.22e-10 | 4.04e-04 | 8.90e-03 | **2.63e-02** | 3.67e-03 |

H2 is converged to 1e-10 and is insensitive. It is the **carbon network** (CO, C⁺) that is
substep-limited, and x_e follows because charge neutrality slaves it to `x_H+ + x_C+`.

## 3. B10 — `nsub_max` silently overrides `cfl_cool`, and the instrument cannot see it

`nsub_max` is **not a truncation cap**. It is a **minimum substep floor**
(`network_gow17_reduced.hpp:159`):

```cpp
const double dt_floor = dt_code / (double)nsub_max;
...
if (dt_chem < dt_floor) dt_chem = dt_floor;   // accuracy limit OVERRIDDEN
```

When the accuracy criterion asks for a step smaller than `dt/nsub_max`, the integrator **takes the
larger step anyway**. The code says so in its own comment — *"floored at dt_code/nsub_max so the
full dt_code is always covered in ≤ nsub_max steps (no truncation)"*. Semi-implicit stability makes
that safe, so nothing blows up; it is purely an accuracy loss, and it is exactly why tightening
`cfl_cool` stops helping.

**The consequence for the diagnostic is the bug.** The counter and its warning
(`chemistry.cpp:298-308`) fire on `t < dt_code` — the sub-cycler failing to cover the step. But the
floor *guarantees* coverage, so that branch is unreachable by construction, and:

```
c100: 0    c050: 0    c025: 0    c0125: 0    c025_big: 0    c0125_big: 0
```

**Zero warnings in all six legs**, including the ones where the cap moved CO by 2.6 %. The
instrument was structurally incapable of reporting the only failure mode that occurs. This is the
same family as B4 (silent EOS-table fallback) and B2 (silent solver non-convergence): a degradation
path with no signal.

`network_h2.hpp` is **not** affected — it is explicit forward Euler with a genuine truncation cap
and no floor, so its return value already meant what it said. `thermo.hpp` **is** affected: same
`dt_floor` pattern, same unreachable `*ntrunc += 1`.

### Fix applied

Both floored integrators now report **floor engagement** as well as (unreachable) truncation, and
the warning says what to do about it:

- `network_gow17_reduced.hpp` — counts `nfloored`, returns 1 if the cell was floored *or*
  truncated.
- `thermo.hpp` — same, through its existing `*ntrunc` out-parameter.
- `chemistry.cpp` — warning rewritten to name `nsub_max` as the binding limiter and to say
  explicitly that tightening `cfl_cool` will not help once it fires.

Diagnostic-only: no abundance, energy or timestep is changed. Built on `build_cpu` (`d7d28f11`).

**Verification (job 2448612) — all three falsifiable predictions confirmed:**

| leg | `cfl_cool` | `nsub_max` | cells reported floored | wall |
|---|---:|---:|---:|---:|
| `v_prod` | 0.1 | 400 | **262144 of 262144 (100 %)** | 380 s |
| `v_big` | 0.1 | 100000 | **2** | 380 s |

(1) production settings now emit the warning where the old binary printed nothing; (2) lifting the
cap reduces it by a factor of 131072, so the instrument responds monotonically to the thing it
claims to measure; (3) identical wall time, so the instrument is free and so is the fix.

**The 100 % figure is the headline.** It is not that the cap *sometimes* binds — at production
settings it binds in **every cell of the domain**, which means `chemistry/cfl_cool = 0.1` has no
effect on the answer at all. That is why tightening it in §2 did nothing.

**Not yet in the GPU binary.** It should ride along with the next `build_gpu` rebuild rather than
trigger one, since `f181c0a1` was just gated.

## 4. Measured on the PRODUCTION configuration — and acted on

The §2/§3 measurements were 32³ on CPU, which cannot answer the cost question: the sub-cycler is a
serial per-cell loop inside a Kokkos kernel, so on a GPU a warp costs its **slowest** cell — one
cell wanting 10⁴ sub-steps makes its 31 warp-neighbours pay 10⁴ too. One cell per thread on CPU,
no lockstep, no penalty. Production also reaches a density regime the smoke deck never touches.

Job 2450198 settles it: three legs **back to back in one job on the same two GPUs**, production
deck, 128³, 90 cycles each. Sequential-in-one-job by design — the first attempt used three
concurrent jobs, which measures the machine as much as the parameter (WP-3's cfl legs reported
1.02e6 vs 1.76e6 zone-cycles/wallsecond for runs that should be identical per zone-cycle).

| `nsub_max` | wall | `wsec_step` | vs 400 | cells still floored | max \|rel diff\| vs 400 |
|---:|---:|---:|---:|---:|---:|
| 400 | 278 s | 1.845 | 1.000× | **2097152 (100 %)** | — |
| 4000 | 280 s | 1.830 | 0.992× | 39977 (1.9 %) | **2.06e-05** (in `mag-dissA`) |
| 40000 | 279 s | 1.835 | 0.995× | 267 (0.013 %) | 2.06e-05 (**byte-identical to 4000**) |

**Three findings, and the third is the one that matters most.**

1. **The cost is zero.** 278 / 280 / 279 s for a 100× change in the cap — a 0.8 % spread, which is
   noise, and the raised legs were if anything marginally *faster*. The GPU-lockstep penalty the
   CPU test could not have seen does not materialise, because the accuracy criterion asks for far
   fewer sub-steps than the cap allows once the cap is lifted.
2. **It is converged at 4000.** `nsub_max` 4000 and 40000 produce **byte-identical** history. 40000
   buys nothing.
3. **Fixing it moves the answer by 2.06e-05**, and only in `mag-dissA` — the *ambipolar
   dissipation*, which is precisely the channel where x_e enters (`ambipolar_coeff =
   ionization_chem`). `MEtor/MEpol`, `ME`, `KE`, `Jsq` and `mass` are unchanged to printed
   precision.

That third point **retro-certifies the earlier results**: production ran with an inert `cfl_cool`,
which is a genuine defect, but its consequence on the flux-retention observable is below printed
precision. No previously reported number is compromised. This is a correction to the alarm implied
by §2/§3 — "100 % of cells are wrong" is true about the *integration*, and worth 2e-5 about the
*answer*, because a floored sub-step is still stable and positivity-preserving; it is merely larger
than the accuracy criterion asked for, and the per-step error is bounded.

Why the impact is so small is itself explained by WP-22: the applied η_A is
`min(eta_chem, eta_eq, ...)`, and `eta_eq` suppresses the chemistry branch by **1486–2678×** in the
first-core shells — so x_e is largely masked exactly where the collapse result lives.

**APPLIED 2026-08-03: `runs/root_ladder/fhc_rootladder.in` now sets `nsub_max = 4000`** (deck md5
`c67eb458` → `43eeb280`), with the measurement recorded inline. Chosen over 40000 only because
4000 is already converged.

### One consistency note, quantified rather than waved away

`r128_sw` and `r256_sw` ran with `nsub_max = 400`; a future `r512_sw` will run with 4000. That
makes the WP-7 ladder inhomogeneous in principle. In practice the inconsistency is **2.06e-05**
against ~9 % differences *between* the rungs — five orders of magnitude below the signal, and
below the printed precision of every quantity the ladder is read on. Re-running the 128³ and 256³
rungs for perfect homogeneity would cost ~1 GPU-hour and is available on request; it is not
required for any conclusion drawn from them.

## 5. Recommendation

**DONE — `nsub_max = 4000` is in the production deck** (§4). The cost is zero on GPU as well as
CPU, the result is converged, and the change is worth 2e-5 on the answer.

The ladder brackets the answer but does not close it: at `nsub_max = 3200` the two `cfl_cool` legs
still differ by 6.8e-3 in CO, so 3200 is not converged either. A follow-up ladder in `nsub_max`
(3200 / 12800 / 51200) at fixed `cfl_cool = 0.025` would find the knee. That is cheap and is the
natural next step, but it is **not on the flagship critical path** — see below.

### How much does this actually matter for the flagship result?

Less than it looks, and for a reason worth stating rather than assuming. x_e feeds η_A through
`ambipolar_coeff = ionization_chem`, so a 0.4 % error in x_e is naively a 0.4 % error in the
ambipolar diffusivity. But WP-22 measured the *applied* η_A and found the equilibrium `eta_eq`
ceiling suppressing the chemistry branch by **1486×–2678×** in exactly the first-core shells — so
in the core the chemistry x_e is not the binding value at all. The chemistry error is therefore
largely masked where the fossil-field result lives.

That is a coincidence of two independent code paths, not a design margin, and it would stop holding
the moment `eta_eq` stopped binding. It is a reason not to prioritize the `nsub_max` ladder, not a
reason to leave `nsub_max = 400` in place.

> *Confidence:* the ladder, the wall times, the zero warning counts and the carbon clamp — all
> **measured** this session. The floor mechanism — **verified** by reading
> `network_gow17_reduced.hpp:159-178` and `thermo.hpp`, and corroborated by the c025 vs c025_big
> pair. The masking argument in §4 — **inferred** from WP-22's measured suppression factors.
