# WP-16 part 3 — the Hall operator in 3D: a real instability, and why production is safe

**Status: the gap is closed, the bug is characterised, and production is measured to sit inside
the stable regime with margin.**

- **3D Hall is unstable at the production Ohmic floor** (`hall_ohmic_floor_code = 0.05`) when
  `η_H` is large. A circularly polarised eigenmode — which the Hall term must conserve exactly —
  grows by **five decades** and then crashes.
- **1D is completely unaffected**, at N = 128, 256, 512 and 1024. `Hall was only ever validated in
  1D`; the recorded historical "0.4 %" is reproduced here as 4.23e-03 by the 1D deck.
- **The fix is a floor/η_H ratio, not a resolution-dependent number**: raising the floor
  0.05 → 0.1 restores stability at both N = 64 and N = 128, with the two resolutions then agreeing
  to ~5 %.
- **The threshold is sharp: floor/η_H between 0.100 and 0.125.** One step across it moves the
  amplitude error from +2.65e5 to +1.6 %.
- **Production is safe by ~10×.** Measured over 78.7 M cells of `prod_v9`: `max|η_H| = 4.316e-02`,
  so `floor/|η_H| ≥ 1.16` everywhere against a threshold of ~0.11. **0.0000 % of cells exceed the
  0.05 floor.**

---

## 1. Scope — this is the half of WP-16 the plan actually asked for

`VALIDATION_PLAN.md` WP-16 reads *"Complete the Hall C-shock test […] Hall is the least-validated
of the three non-ideal terms."* `WP16_ambipolar_convergence.md` covers the **ambipolar** operator,
which had to be untangled first because the deck named there (`cshock_ad.in`) is ambipolar and
turned out to be invalid. This document is the Hall half.

Method is deliberately identical to the ambipolar half: measure against an **analytic** solution,
not by self-convergence. `src/pgen/diffusion.cpp` `iprob = 60` initialises a circularly polarised
Hall eigenmode (`B_x = B₀`, `B_y = amp cos kx`, `B_z = h·amp sin kx`) with a closed-form dispersion
relation, and reports `|ω_meas − ω_analytic|/ω_analytic` plus the in/out amplitude. The amplitude
is the stronger check: **the Hall term is non-dissipative**, so a circular mode must conserve
`⟨B_y² + B_z²⟩` regardless of whether the dispersion relation is right.

## 2. The finding — 1D works, 3D does not

Stock decks, **zero overrides**, all at their native N = 128 (job 2449288 Part A):

| deck | geometry | ω rel. error | amplitude in → out |
|---|---|---|---|
| `hall_whistler.in` | **1D** (nx2 = nx3 = 1) | 4.23e-03 | 1.00e-06 → 9.08e-07 (−9.2 %) ✓ |
| `hall_whistler_glm.in` | **3D** (nx2 = nx3 = 4) | 6.23e-01 | 1.00e-06 → **2.65e-01** ✗ |
| `hall_whistler_ct.in` | **3D**, CT | 6.23e-01 | 1.00e-06 → **2.65e-01** ✗ |

The two 3D variants fail **identically to every printed digit**, so this is not about the
divergence-control scheme — one is GLM, the other CT.

**Gap, not regression.** The 1D deck's 4.23e-03 reproduces the historically recorded *"Hall
whistler branch 0.4 %"*. That validation was 1D; the 3D case was evidently never checked against
the analytic ω.

### Two alternative explanations ruled out

**Not a timestep bug.** With the ladder corrected so all three directions refine together, dt
scales exactly as dx²: 5.859e-04 / 1.465e-04 / 3.662e-05 / 9.155e-06 for N = 32/64/128/256. The
Hall dt limiter is doing its job.

**Not resolution alone.** 1D is stable and accurate at every resolution tested:

| N (1D) | ω rel. error | amplitude change |
|---|---|---|
| 128 | 4.23e-03 | −9.24 % |
| 256 | 4.36e-03 | −9.24 % |
| 512 | 4.39e-03 | −9.24 % |
| 1024 | 4.39e-03 | −9.24 % |

whereas 3D degrades monotonically with resolution — amplitude change −9.3 % (N=32) → **+76 %**
(N=64) → **+2.65e5** (N=128) → crash (N=256), on both helicity branches.

*(Note the 1D ω error is flat, not converging. The 1D test therefore gives a stability-and-magnitude
check, not an order of accuracy — the residual is dominated by some fixed source. Establishing an
order for Hall remains open.)*

### A mis-specification of my own, corrected

The first attempt (job 2449287) refined `nx1` only. The deck's transverse directions are 4 cells
across a *fixed* extent, so `dx2 = dx3 = 0.0078125` always; since dt is set by `min(dx)` over all
directions, **dt came out identical (3.6621e-05) at N = 32, 64 and 128** — not a temporal
convergence ladder at all. `dx2 = dx1` exactly at N = 128, so the stock deck is self-consistent at
its native resolution; the override broke that, not the deck. Job 2449288 scales the transverse
extent as ±2/N so every direction refines together. The conclusion survived the correction.

## 3. Mechanism — the Ohmic floor is under-calibrated, and the requirement is a RATIO

Job 2449353 Part B holds the geometry fixed and scans `hall_ohmic_floor_code`:

| floor | N = 64: ω err / amplitude | N = 128: ω err / amplitude |
|---|---|---|
| **0.05** (production) | 3.80e-03 / **+76 %** ✗ | 6.23e-01 / **+2.65e5** ✗ |
| **0.1** | 9.45e-03 / −17.6 % ✓ | 9.95e-03 / −17.6 % ✓ |
| 0.2 | 2.46e-02 / −32.1 % ✓ | 2.51e-02 / −32.1 % ✓ |
| 0.4 | 6.59e-02 / −53.7 % ✓ | 6.63e-02 / −53.7 % ✓ |
| 0.8 | 1.63e-01 / −78.0 % ✓ | 1.64e-01 / −78.0 % ✓ |

Two things follow.

1. **The threshold sits between 0.05 and 0.1**, and once above it the answer is
   **resolution-converged** — N = 64 and N = 128 agree to ~5 % on ω and to three digits on the
   amplitude. So the requirement is a fixed **ratio** `η_floor/η_H`, not a resolution-dependent
   number. With the stress value `η_H = Q_H B/ρ = 0.5`, the threshold ratio is between **0.1 and
   0.2**.

   This *kills* the sub-hypothesis I started with — that the required floor scales as 1/dx. It
   should not have: Hall damping (`η_H k²`) and Ohmic damping (`η_O k²`) carry the same power of k,
   so a fixed ratio is exactly what theory predicts, and the resolution-dependent *onset* is simply
   where the marginal case at 0.05 tips over.

2. **The stabiliser is not free.** It is a real resistivity in the dispersion relation, so buying
   stability costs Hall fidelity: at floor = 0.1 the mode damps 17.6 % and ω error rises from
   4e-3 to 1e-2; at 0.8, damping is 78 % and ω error 0.16. Raising the floor "to be safe" is not a
   neutral act.

## 4. Is production affected? — measured, not argued

The test's `η_H = 0.5` is a deliberately large stress value (`hall_coeff = fixed`). Production runs
`hall_coeff = ionization`, so `η_H` is set by the ionization model and varies. Measured directly
from the applied `nonideal_eta` cache on `prod_v9` (t = 1.0946, 2402 blocks, 78.7 M cells):

```
max |eta_H|                                     = 4.316e-02
cells with |eta_H| > 0.05 (the floor)           = 0.0000 %
cells with |eta_H| > eta_O                      = 94.8 %
max eta_O                                       = 1.000e-01  (= eta_ohm_cap_code)
```

**The floor exceeds `|η_H|` in every cell of the domain** — `η_floor/η_H ≥ 0.05/0.0432 = 1.16`
everywhere, against a measured threshold ratio of 0.1–0.2. That is a margin of **6–12×** on the
ratio, and production's largest `η_H` is **11.6× below** the stress value at which the instability
was provoked.

### The threshold, measured directly

Job 2449887 holds the floor at production's 0.05 and scans `Q_H` (hence `η_H`), at N = 128:

| `Q_H` | floor/η_H | ω rel. error | amplitude change | |
|---|---|---|---|---|
| 0.50 | **0.100** | 6.23e-01 | **+2.65e+05** | ✗ unstable |
| 0.40 | **0.125** | 6.61e-03 | **+1.6 %** | ✓ |
| 0.30 | 0.167 | 1.08e-02 | −9.22 % | ✓ |
| 0.25 | 0.200 | 1.42e-02 | −9.21 % | ✓ |
| 0.20 | 0.250 | 1.90e-02 | −9.20 % | ✓ |
| 0.15 | 0.333 | 2.56e-02 | −9.18 % | ✓ |
| 0.10 | 0.500 | 3.46e-02 | −9.16 % | ✓ |
| 0.05 | 1.000 | 4.67e-02 | −9.13 % | ✓ |

**The onset is sharp and sits between floor/η_H = 0.100 and 0.125.** One step across it takes the
amplitude error from +2.65e5 to +1.6 %. That is a threshold, not a gradual degradation, and it
confirms the ratio picture from the other direction — the same boundary appears whether the floor
is raised at fixed `η_H` or `η_H` lowered at fixed floor.

Note also the accuracy cost is monotone in the ratio (ω error 6.6e-3 → 4.7e-2 from 0.125 to 1.0),
so the best place to sit is *just* above threshold, not far above it.

**Production margin: floor/|η_H| ≥ 0.05 / 4.316e-02 = 1.16 in every cell, against a threshold of
~0.11 — roughly 10× inside the stable regime.**

So: the instability is real and the code should be fixed, but **no production result is
compromised by it**.

Two caveats kept in view: (i) `η_H` is a *local* quantity and a future run reaching different
ρ/T/B could exceed the floor — the safety is empirical, not structural; (ii) production runs
`hall_floor_integrator = rkl2`, which super-time-steps the floor, whereas this test is pure
`unsplit`. The comparison is like-for-like on `η_H` magnitude but not on the floor's integration.

## 5. The fix, implemented and gated — `diffusion/hall_ohmic_floor_ratio`

The applied stabilizer is now per-cell:

```
eta_floor_cell = max(hall_ohmic_floor_code, hall_ohmic_floor_ratio * |eta_H_cell|)
```

`max()`, not replace, so raising the ratio can only ADD dissipation, never remove it. Default
`0.0` = disabled = the old constant. Implemented in `HallDiffusivity::EffectiveOhmicFloor`
(`src/hydro/diffusion/diffusion.hpp`).

### Every consumer was updated, not just the EMF

An EMF-only change would have produced a stable-but-wrongly-stepped run — the floor is a real
parabolic diffusivity and sets a `dx²/η_floor` constraint:

| file | site |
|---|---|
| `hall.cpp` | 3× cell-centred EMF (X1/X2/X3), GLM path |
| `hall.cpp` | `EstimateHallTimestep` — the parabolic constraint |
| `ct.cpp` | 3× edge EMF (E1/E2/E3), CT path |
| `diffusion.cpp` | fused **unsplit** dt estimator (`eta_O_tot = eta_O + floor`) |
| `diffusion.cpp` | fused **mixed RKL2** dt estimator (`eta_floor_par` / `eta_floor_strict`) |
| `hydro.cpp` | all three coefficient branches (`fixed`, `ionization`, `ionization_chem`) |

**One trap worth recording.** In the mixed RKL2 mode — production's integrator — the Hall kernel
is called twice: once `(eta_h_on=true, floor_on=false)` for the dispersive part, and again
`(eta_h_on=false, floor_on=true)` for the floor. On that second call `eta_H` was never evaluated,
so a naive implementation would have read `|eta_H| = 0`, silently fallen back to the absolute
floor, and been **inert in exactly the integrator that matters**. A `need_eta_h` gate now forces
the evaluation for the floor while leaving the dispersive coefficient gated on `eta_h_on`, so the
mode's split is preserved exactly.

### Gate results (job 2450082, `build_cpu` `549f7c28`)

**A — OFF-state, ratio unset, vs the pre-B11 binary `2ddea223`:**

```
A1 GLM unsplit Hall: PASS byte-identical
A2 CT Hall        : PASS byte-identical
```

**B — ON-state, does ratio = 0.2 fix the instability?**

| N | ratio | ω rel. error | amplitude in → out |
|---|---|---|---|
| 128 | 0.0 | 6.23e-01 | 1.00e-06 → 2.65e-01 (**+2.65e5**) ✗ |
| 128 | **0.2** | **8.30e-03** | 1.00e-06 → 8.24e-07 (−17.6 %) ✓ |
| 256 | 0.0 | — | **run aborts** (PARTHENON ERROR) ✗ |
| 256 | **0.2** | **8.42e-03** | 1.00e-06 → 8.24e-07 (−17.6 %) ✓ |

The 256³ case goes from *crashing* to well-behaved, and the two resolutions agree to 1.4 % on ω
and to three digits on amplitude — resolution-converged, as the ratio picture requires.

A detail that confirms the mechanism rather than just the outcome: the earlier *constant*-floor
scan at 0.1 gave ω error 9.95e-03, while the per-cell ratio at 0.2 gives **8.30e-03** — slightly
better. Expected, not noise: `η_H = Q_H·B/ρ` varies across the domain as the wave passes, so a
per-cell `0.2·|η_H|` applies less dissipation where `η_H` is small than a flat 0.1 does. The ratio
buys the same stability for less fidelity cost, which is the point of making it track the term it
stabilizes.

**C / D — production equivalence.** Production runs floor = 0.05 with `max|η_H| = 4.3e-2` measured
over 78.7 M cells, so `0.2·|η_H| ≤ 8.6e-3 < 0.05` and the `max()` must select the absolute floor in
every cell ⇒ enabling the ratio must change nothing. Both pass:

```
C production unsplit: PASS byte-identical
D production rkl2   : PASS byte-identical
```

D is the one that matters most — it is production's actual integrator and the path where a naive
implementation would have been inert (see the `need_eta_h` trap above).

### Operational note — enabling it on an EXISTING run needs the CLI, not the deck

A restarted AthenaPK run takes its parameters from the `.rhdf` + the command line, **not** from the
input deck (`runs/prod_v9/fhc.in` says so in its own header, and `runs/root_ladder/submit_root.sh`
mirrors every deck change as a CLI override for exactly this reason). So adding
`hall_ohmic_floor_ratio` to a deck affects **fresh starts only**; to enable it on a chain already
in progress it has to be passed as `diffusion/hall_ohmic_floor_ratio=0.2` on the `mpirun` line.
This is the same trap that made the `turb_ksample` key silently produce a superseded IC earlier in
this campaign.

The parameter itself is restart-safe: it is read in `Hydro::Initialize` via `GetOrAddReal`, which
runs on restart, so it is not subject to the B5 class of defect (Params registered only in
`ProblemGenerator`, which restart never calls).

Startup banner is now adaptive: it reports the per-cell formula when a ratio is set, **warns if the
ratio is below 0.15** (measured onset ~0.11), and otherwise keeps the "your floor is absolute"
notice.

## 6. Recommendation

**Set `diffusion/hall_ohmic_floor_ratio = 0.2` in the production deck.** It costs nothing today — gates C and D prove byte-identity at production's current `η_H` — but it
converts an *empirical* safety margin into a *structural* one. Right now the run is stable because
`max|η_H|` happens to sit below 0.05; with the ratio set it is stable because the code guarantees
it, at any density the run subsequently reaches.

Why 0.2 and not larger: the measured onset is ~0.11, so 0.2 carries ~2× margin, and the stabilizer
is not free — ω error rises monotonically with the ratio (6.6e-3 at 0.125 → 4.7e-2 at 1.0). Sit
just above threshold.

**Diagnostic added (job 2450078 gates it, and it is now superseded by the fix above).** A startup notice now states the criterion and how to
check it, because the three existing warnings fire only when the floor is *absent* — a floor that
is present and too small was silent. Gated the same way as every other change this session: rerun
the identical deck on the new binary and byte-compare the history against what the previous binary
produced.

```
new binary: 2ddea223 (was d7d28f11)
PASS: history BYTE-IDENTICAL -> the B11 notice is an OFF-state no-op.
notice present in the new run: 1
```

**A per-cell instrument would be better still.** The same reasoning as B2/B4/B10: a cell whose
`|η_H|` exceeds its applied floor is running an unstable configuration and nothing currently says
so. `cap_diag` already reports Ohmic/Hall/AD ceiling engagement; the Hall floor deserves the same
treatment.

> *Confidence:* every number above is **measured** this session — the stock-deck comparison, the
> 1D ladder to N = 1024, the floor scan at two resolutions, and the production `η_H` distribution
> read from the code's own applied-diffusivity cache. The mechanism (fixed ratio, same power of k)
> is **inferred** but is consistent with all three scans. The `rkl2`-vs-`unsplit` caveat in §4 is
> **unmeasured**.

---

# 2026-08-08 — the recommendation is now the DEFAULT

§6 recommended setting `diffusion/hall_ohmic_floor_ratio = 0.2` in the production deck. That
recommendation stood unexecuted: **0 of 295 decks set it**, so the ratio existed but nothing used
it, and the closing statement "production is measured to sit inside the stable regime with margin"
remained an *empirical* safety claim about one measured snapshot rather than a structural property
of the code.

**Changed:** the input default is now `0.2` at all four read sites in `src/hydro/hydro.cpp`
(`fixed`, `ionization`, `ionization_chem`, and the startup-banner probe). The C++ constructor
fallback in `diffusion.hpp` stays `0.0`; every call site passes the input value explicitly.

## Why this is free today, and why it is not cosmetic

Free: `max|η_H| = 4.316e-02` measured over 78.7 M cells of `prod_v9`, so `0.2·|η_H| ≤ 8.6e-03 <
0.05` and the `max()` selects the absolute floor in **every** cell ⇒ no change to production
numbers. This is not an argument, it is §5's gates C (unsplit) and D (rkl2 — production's actual
integrator, and the path where a naive implementation would have been inert), both **PASS
byte-identical**.

Not cosmetic: η_H is a *local* quantity. §4 already flagged caveat (i) — "a future run reaching
different ρ/T/B could exceed the floor; the safety is empirical, not structural". The default flip
is exactly what converts that caveat into a guarantee: once the ratio is on, the floor tracks the
term it stabilises, so a run that drives |η_H| up carries its own stabiliser up with it instead of
silently crossing the measured onset at η_floor/|η_H| ≈ 0.11.

## What this DOES change

The `hall_whistler*` validation decks run η_H = 0.5, so the floor becomes `max(0.05, 0.1) = 0.1`
and their recorded numbers move. In particular the historical **1D "4.23e-03"** in §2 is a
ratio = 0.0 result. Reproduce it with `diffusion/hall_ohmic_floor_ratio = 0.0`.

That is the intended trade: the stock decks become **stable in 3D by default**, and reproducing the
historical unstable configuration is now the opt-in. Given that §2's headline finding is that those
same stock 3D decks amplify by five decades and crash at N = 256 with zero overrides, defaulting
them to the stable branch is the correct polarity.

## Verification

`runs/audit_fix_regress/submit_default_flips.sh` — one binary, new defaults vs the old values
restored explicitly on the command line, on `A_multipole.in` (32³, `hall_coeff = ionization`,
`hall_ohmic_floor_code = 0.05`, radiation on). Includes a falsification leg
(`hall_ohmic_floor_ratio = 5.0`) so that a bit-identical result cannot be produced by the parameter
name simply being ignored.

## Verification of the default flip (jobs 2492381 / 2492597 / 2492626)

**Result: the flip is INERT at production-like eta_H, and the mechanism is demonstrably live.
Verified, but only after a falsification leg failed and had to be explained.**

`A_multipole.in`, 32^3, `hall_coeff = ionization`, `hall_ohmic_floor_code = 0.05`,
`diffusion/integrator = unsplit`, 12 cycles. Measured `max|eta_H|` on the deck itself:
**2.82e-02** (cycle 6) rising to **5.73e-02** (cycle 12).

| comparison | result |
|---|---|
| ratio 0.2 vs 0.0 (the flip) | **bit-identical** in `prim`, `grav.phi`, `rad.Er` |
| ratio 0.2 vs 5.0 | **bit-identical** |
| ratio 0.2 vs 1.0e4 | differs, rel 8.98e-01, all 3 670 016 cells |
| ratio 0.2 vs 1.0e6 | differs, rel 8.96e-01 |
| floor 0, ratio 0 vs 0.2 | differs (12 hst rows) |
| floor 0, ratio 0 vs 1e4 | differs (13 hst rows) |

### The falsification leg failed first, and that was informative

`ratio = 5.0` gives `EffectiveOhmicFloor = 5.0 x 5.73e-02 = 0.286`, comfortably above the 0.05
absolute floor, so it *should* have changed the answer — and it did not, in the fields, not merely
in the 6-significant-figure `.hst`. Three things were checked before accepting it:

1. **The parameter reaches the code.** The startup banner in that leg reads
   `per-cell max(0.05, 5*|eta_H|)`.
2. **The plumbing is complete.** The `ionization` branch passes the ratio into `HallDiffusivity`;
   the EMF kernels (`hall.cpp` 238/311/380, `ct.cpp` 1845/1894/1937) and **both** dt estimators
   (`diffusion.cpp` 177-179 and 290-293/305-308) apply the per-cell
   `fmax(eta_floor, ratio*|eta_H|)`. *(An earlier reading of this file that claimed the dt
   estimators used only the absolute `GetOhmicFloor()` was WRONG — it looked at the fetch lines
   124/228 and missed the per-cell application below them. §5's table is accurate.)*
3. **The zero-floor discriminator.** With `hall_ohmic_floor_code = 0` there is no absolute floor
   left to mask the ratio, and `ratio = 0` vs `0.2` then differs. So the ratio does reach the
   solution.

**Explanation.** The floor enters as `eta_floor * J`. On this deck the field is still nearly the
uniform `B0z` after 12 cycles, so `J ~ 0` precisely in the low-density envelope where `|eta_H|` is
largest — the EMF contribution is negligible at both 0.05 and 0.286 and cancels to the last bit.
What does change at `ratio = 1e4` is the **parabolic timestep**: `dt_par = cfl*dx^2/eta_floor`
with `eta_floor = 573` binds hard, dt collapses, and the whole trajectory moves — which is exactly
the *global* signature observed (all 3.67 M cells, not a localised set).

So `0.2` is inert here for the same reason it is inert in production — `0.2*|eta_H| << 0.05` — and
the guard only engages once `|eta_H| > floor/ratio = 0.25`, about **5x production's measured
maximum**. That is the runaway case the flip exists to catch.

### What this does NOT establish

The instability itself is not re-tested here; §3/§5 did that. This confirms only that the new
default perturbs nothing at present `eta_H` and that the knob acts when `eta_H` is large enough.
A deck that is both hot and strongly Hall-dominated would exercise the EMF path directly; none
exists in the suite.
