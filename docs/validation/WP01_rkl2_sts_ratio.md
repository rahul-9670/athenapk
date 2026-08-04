# WP-1 — RKL2 super-time-stepping ratio (`rkl2_max_dt_ratio = 1000`)

**Status: PASS — CLOSED 2026-08-04.** Four binding legs (caps 30/15/8/4, jobs 2453427 + 2453441)
plus a fifth at cap = 2 (job 2453562) as a check on the one loose end. Production's
`rkl2_max_dt_ratio = 1000` **never binds at any state this campaign has reached**, and forcing a
**9.2× tighter** super-timestep moves the flux-retention observable by **+0.001 %** and the
dominant dissipation channel by **≤0.08 %**. The only quantity that responds is the Ohmic
dissipation — **0.106 %** of the magnetic dissipation budget — which moves ~17 % non-monotonically.

**Verdict: production needs no change.** The scope limit stands as originally written and is
**not** closed: the cap-binding epoch (`prod_t4_full` near cycle 72410) was never checkpointed
among the surviving restarts, so it cannot be probed. Two rounds of probing established only that
the reachable states are not that state — the achieved ratios of 36.7 and 51.2 said so directly.
Crucially, at the one state where old and current CAN be compared, they **agree** (old 43–83,
current 51.2), so there is no regression and no anomaly to explain.

## Why this took three attempts to become testable at all

A cap can only be tested where it **binds**. The achieved STS ratio is what decides that, and it
was far too small in every earlier configuration:

| attempt | configuration | achieved max STS | verdict |
|---|---|---:|---|
| job 2450388 | 128³ uniform | **4.04** | **null test** — caps 1000/250/60 gave identical results; no cap in that range can bind on a ratio of 4 |
| job 2453201 | 128³ + AMR, deep | 3.75 → 6.45 → 39.4 → **95.7** → 82.0 | reached the regime, then hit the dt-wall (dt = 8.1e-7, 32 s/cycle) — stopped, ~18 GPU-h for no further answer |
| jobs 2453427 / 2453441 | restart from the deep state | **36.7** | **testable**: caps of 30/15/8/4 all bind |

The deep run's own history is why WP-1 is answerable now. Note the achieved 36.7 at the restart
is *lower* than the 95.7 the parent run reached later: restart `00001` sits at cycle 250, and the
deeper cycle-425 state was never checkpointed (`dn = 250` would next have fired at cycle 500).
**36.7 is what is reachable from the deepest restart that exists.**

## The test

Restart `deep_amr/run/parthenon.out2.00001.rhdf` — t = 1.0145874, cycle 250, **792 blocks** —
run to `tlim = 1.0150` under `build_gpu_v3` (`6b1fe753`), 4 GPUs, production physics.

**Compared at fixed TIME, not fixed cycles.** A smaller cap forces a smaller dt, so equal cycle
counts would compare different physical times. The verdict enforces this with a hard guard: any
leg whose last history row is not at `tlim` is reported INCOMPLETE and excluded. This was not
hypothetical — a mid-flight `cap = 8` leg at t = 1.01460 against the control's 1.01500 produced an
apparent **−59.75 %** on `mag-dissO` that was pure artefact of the time mismatch.

**The primary discriminant is the dissipation, not the flux ratio.** The window is
Δt = 4.13e-04, i.e. **0.041 % of the run**; over so short a span `MEtor/MEpol` barely evolves under
any cap, so its agreement alone would be weak evidence. RKL2 super-time-steps the **parabolic**
(ambipolar + Ohmic) terms specifically, so `mag-dissA` / `mag-dissO` are what the cap directly
controls. `mag-dissO = ∫ η_O J² dV` is an *instantaneous* volume integral (`mag_diag.cpp:132`),
not accumulated heating, so matched-time endpoints are a fair comparison.

## Results (all at matched t = 1.015000)

| cap | throttle | cycles | achieved | `mag-dissO` | vs control | `mag-dissA` | `MEtor/MEpol` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| **1000** | 1.0× | 14 | 36.7 | 6.626900e-04 | — (control) | — | 1.196092e-02 |
| 30 | 1.2× | 16 | 30.0 | 6.469360e-04 | −2.377 % | −0.027 % | +0.0002 % |
| 15 | 2.4× | 36 | 15.0 | 5.737140e-04 | −13.426 % | −0.017 % | +0.0013 % |
| 8 | 4.6× | 107 | 8.00 | 5.413090e-04 | **−18.316 %** | +0.036 % | +0.0013 % |
| 4 | 9.2× | 269 | 4.00 | 5.525170e-04 | −16.625 % | +0.080 % | +0.0013 % |

Every capped leg reports its achieved ratio pinned exactly at its cap, so engagement is positive,
not assumed.

### Null leg and determinism control

`cap = 60` (round 1) is **above** the achieved 36.7 and therefore cannot bind. It came back
**identical to the control except one number** — `grav-res`, 3.45999e-05 vs 3.45967e-05 (9e-5
relative, the cancellation-sensitive multigrid residual). Every state and dissipation column
matched exactly. That confirms both the null-leg prediction and that the harness is deterministic,
which the whole A/B rests on.

*(Round 1's engagement check was too weak to notice this on its own: it asked only whether the
achieved ratio EXCEEDED the cap, so 36.7 ≤ 60 printed "cap engaged" for a cap that restricted
nothing. The correct test — `cap < control's achieved max` — labels it NULL.)*

## Reading the Ohmic trend

A smaller cap means *more* substeps, i.e. a **better-resolved** parabolic term, so the control is
the **least** converged point and the sequence should be read as a convergence test. The full
four-point sequence:

| throttle | `mag-dissO` vs control | increment |
|---:|---:|---:|
| 1.2× | −2.377 % | −2.377 |
| 2.4× | −13.426 % | −11.049 |
| 4.6× | −18.316 % | −4.890 |
| **9.2×** | **−16.625 %** | **+1.691** ← turns back |

**The sequence is NOT monotone.** It deepens to −18.3 % at 4.6× and then *recovers* to −16.6 % at
9.2×. So this is not a clean first-order convergence sequence, and no reliable asymptote can be
extrapolated from it. The defensible statement is bounded, not extrapolated: **under throttling,
`mag-dissO` settles roughly 17 ± 2 % below the uncapped value**, with a non-monotone approach —
consistent with `J²` being dominated by small-scale current structure that is not smoothly
convergent in the super-timestep at this resolution.

> **Two corrections, both recorded because each was reported before it was checked.**
> 1. On two points alone (−2.38 %, −13.43 %) this was written up as *"the change grew five-fold ⇒
>    NOT converged"*. Wrong inference from insufficient data.
> 2. On three points it was then written up as *"increments halve per doubling ⇒ first-order
>    convergence, asymptote ≈ −24 %"*. **The fourth point falsifies that**: the increment changes
>    sign (+1.691). The −24 % extrapolation was never valid.
>
> A three-point trend does not establish an asymptote, and a monotone-looking start does not
> establish monotonicity. Both claims should have waited for the full ladder.

### The wobble is deterministic, not scatter

Worth stating because it changes how the non-monotonicity should be read: the `cap = 60` null leg
came back identical to the control on every column except the cancellation-sensitive `grav-res`,
which establishes that **this harness is deterministic**. Re-running any leg reproduces its number
exactly. So the −18.3 % → −16.6 % turnaround is **not statistical scatter** — it is a reproducible
property of how a given cap discretizes the parabolic operator. It cannot be averaged away, and a
fifth point (cap = 2, throttle 18.4×) is the only way to see whether it settles or oscillates
further. That leg is job 2453562.

## What this does and does not license

- **Does**: production `rkl2_max_dt_ratio = 1000` needs no change. It never binds at any state
  reached, and even a 4.6× tighter super-timestep leaves `MEtor/MEpol` flat to **1e-5 relative**
  and `mag-dissA` — 99.9 % of the dissipation budget — flat at the ±0.036 % noise level.
- **Does**: document that the STS ratio measurably and *non-monotonically* perturbs the **Ohmic**
  channel. `J²` is dominated by the smallest scales, so this says the super-timestep strides over
  small-scale current structure. At −18 % of a channel that is 0.106 % of the budget, that is
  **−0.019 %** of total magnetic dissipation.

> **A caveat on the script's own PASS wording.** The generated verdict prints *"the parabolic
> dissipation moves at most 0.080 %"* — that figure is `mag-dissA` only, because the summary
> statistic was built from `f('dissA')`. `mag-dissO` is equally parabolic and moved **18 %**, two
> orders larger. The PASS conclusion still stands (it rests on `MEtor/MEpol` and on dissA carrying
> 99.9 % of the budget), but the sentence as printed understates what was measured. Read the table,
> not that line.
- **The level-12/13 scope limit is now CLOSED — see the section below.** Earlier drafts of this
  document said the cap could not be tested where it pins, and that this was the honest residual.
  That residual has since been measured away.
- **Does NOT** cover more than Δt = 4.13e-04 (0.041 % of the run) from one restart point.

---

# The level-13 scope limit — CLOSED (probe job 2453567)

Earlier drafts said the cap could not be tested where it pins, and cited *"DEV_LOG 2026-07-10,
cap pinning every cycle at AMR level 12"*. Both the citation and the conclusion needed fixing.

## The citation was wrong

Reading `DEV_LOG.md` directly rather than from memory:

* **2026-07-11** (not 07-10): *"Production entered level 11 … the `rkl2_max_dt_ratio=`**400**` cap
  now binds every cycle (dt ~1.0-1.8e-7, 29 STS stages/half)"* — level **11**, cap **400**.
* **2026-07-12 evening**: *"at the level-12 state … **STS ratio pinned at 1000** (45 stages/half)"*,
  mean dt 3.09e-8 vs 1.35-1.5e-8 when cap-400-pinned.

Wrong date, and two different (level, cap) pairs conflated into one claim.

## The binding regime was real, and was measured

From the production logs, not assumed:

| run | depth | achieved max STS | cap=1000 binds? |
|---|---|---:|---|
| `prod_t4_full` (2026-07-21) | level **13** | **1.00e+03 = the cap**, on **1996 of 41880** cycles (4.8 %) | **YES** |
| `prod_v9` (2026-07-23) | level 12 | 644 | no |
| `deep_amr` (current, today) | t = 1.016 | 95.7 | no |

## What the probe found — and it inverts the concern

`prod_t4_full/parthenon.out2.00200.rhdf` (t = 1.08464, cycle 85500, 2486 blocks, **maxlevel 13**)
was restarted under the **current binary and current non-ideal settings**, 25 cycles, 5 GPUs:

```
cycles=26   achieved_max_STS = 5.12e+01   pinned_at_1000 = 0 of 50 records
```

**At the identical state where the July run pinned at 1000, current settings achieve 51.2** — a
factor of ~20 lower. The cap does not bind there any more.

This also **falsifies a prediction made in this document**: I argued B11's `hall_ohmic_floor_ratio`
raises η_O and would therefore push the required ratio *up* at equal depth. It went sharply down.

## RESOLVED: there was never an old-vs-new discrepancy — I compared two different states

The whole "current settings give 51.2 where the old run gave 1000" framing is **wrong**, and the
error was mine. Reading `prod_t4_full/run.log` at the **same cycle as the restart** (85500):

```
cycle=85500  t=1.0846417004626012  dt=1.3213e-08     STS ratio: 4.30e+01   (45 stages)
cycle=85501                        dt=2.5558e-08     STS ratio: 8.31e+01
```

The old run at that state achieved **43–83**. The current-binary probe at that identical state
achieved **51.2**. **They agree.** The current binary reproduces the old behaviour exactly where
it can be compared; there is no physics regression, no ionization mystery, nothing to explain.

**Where the 1000-pinning actually happened:** near **cycle 72410, t = 1.0839142** — a
*different, earlier epoch*, ~13 000 cycles before the restart I probed.

**Why it cannot be probed:** every surviving restart post-dates that epoch —

| restart | cycle | t |
|---|---:|---:|
| `out2.final` | 84576 | 1.0846231 |
| `out2.00199` | 85250 | 1.0846366 |
| `out2.00200` | 85500 | 1.0846417 |

The pinning-epoch checkpoints were destroyed when `prod_t4_full` was pruned 6.7 TB → 399 GB.
**So the original scope limit was correct**: the cap-binding state was never checkpointed and
remains unreachable. Two rounds of probing established only that the state I *could* reach is not
that state — which the achieved ratio of 36.7/51.2 said from the start.

## What DOES set the ratio (this part is established)

`hydro_driver.cpp:566-573` defines `ratio = 2·tau / mindt_diff`, i.e. **global dt ÷ diffusive dt**.
Decomposing the two ζ legs at fixed state:

| leg | dt | ratio | ⇒ dt_diff |
|---|---:|---:|---:|
| ζ = 1e-16 | 1.571894e-08 | 51.20 | **3.070106e-10** |
| ζ = 1e-17 | 1.571894e-09 | 5.12 | **3.070106e-10** |

`dt_diff` is **identical to 7 digits**, and the two dt values share the mantissa `1.5718943`
differing by exactly 10×. The cap diagnostics say why:

```
cap-VO/MO/DO  nonzero  -> the OHMIC cap engages
cap-VH/MH/DH  0.0      -> the Hall cap NEVER engages
cap-VA/MA/DA  0.0      -> the ambipolar cap NEVER engages
```

So η_O is pinned at `eta_ohm_cap_code = 0.1` and `dt_diff` is therefore **independent of ζ**,
while η_H is *uncapped* and scales as 1/x_e ∝ 1/ζ, making the Hall-limited global dt ∝ ζ.
Hence **ratio ∝ ζ at fixed state** — confirmed to 8 significant figures.

That is also why my prediction failed: I reasoned that lower ζ raises η and shortens the parabolic
timestep, but the parabolic timestep is *capped* and cannot respond. Lower ζ only shortens the
Hall-limited global dt, dragging the ratio down with it.

### Superseded: the ionization-inconsistency hypothesis

> **Read this before the subsection below, which is kept only as a record of a wrong turn.**
> The `ion_zeta` explanation was tested directly (job 2453604: the same probe with **both**
> `diffusion/ion_zeta` and `chemistry/zeta_cr_cgs` set to 1.0e-17, self-consistently) and it
> **FAILED**:
>
> | ζ | achieved max STS | dt |
> |---|---:|---:|
> | 1e-16 (current production) | 51.2 | 1.57e-08 |
> | **1e-17 (the old default)** | **5.12** | **1.57e-09** |
> | old `prod_t4_full` at the same state | **1000** | — |
>
> The prediction was that lowering ζ raises η, shortens the parabolic timestep and drives the
> required ratio back up toward 1000. Instead the ratio fell by exactly **10×**, *away* from the
> old behaviour, with dt falling by exactly 10× alongside it. That proportionality is a clue, but
> no mechanism here explains it, and none is asserted.
>
> **What stands empirically:** at the level-13 state, current settings give an achieved ratio of
> **51.2**, nowhere near the 1000 cap, and no configuration tried reproduces the old pinning.
> **What is NOT established:** why the July run pinned at 1000. The ionization inconsistency
> below is real and was really fixed, but it is *not* demonstrated to be the cause.
>
> This is the third premature mechanism claim in this document (see the two Ohmic-trend
> corrections). The pattern is the same each time: a plausible story written up before the
> decisive single-variable test was run.

### The (real, but not causally demonstrated) ionization inconsistency

`ionization_environment.hpp:60` — `ion_zeta` defaults to **1.0e-17**. The old `prod_t4_full`
submit sets only `eta_ohm_cap_code`, so it ran at that default. Current production sets
**1.0e-16** (`fhc_rootladder.in:150`, mirrored in `submit_root.sh:83` and `submit_deep.sh:87`) —
a **10× higher cosmic-ray ionization rate**, which raises x_e, lowers the non-ideal
diffusivities, lengthens the parabolic timestep and so reduces the required STS ratio. A 10×
change in ζ against a ~20× drop in achieved ratio is the right order.

*That falsification test was run (job 2453602 blocked by the consistency guard, then 2453604 with
both keys at 1e-17) and the explanation **did not survive** — see the box above. The guard itself
is worth recording: `ion_zeta` cannot now be set inconsistently with `chemistry/zeta_cr_cgs`, so
the old mixed state is unreachable by configuration. That is a real hardening, independent of
whether it explains the STS ratio.*

## Consequence for WP-1

The cap-binding regime documented in the DEV_LOG was a property of the **old ionization setting**,
not of depth per se. Under current production physics the achieved ratio is **51.2 at level 13**
and **95.7** at the highest value ever recorded in any current-configuration run — so
`rkl2_max_dt_ratio = 1000` carries **~10–20× headroom and has never bound**. The A/B legs in the
main table therefore cover the regime production actually operates in, and the earlier "safe in
this regime, not safe everywhere" caveat is resolved rather than merely acknowledged.

*Scope of the probe itself, stated plainly: the restart is gray (`n_group = 1`), 5 scalars, L=52 —
a different problem setup from current production (3 groups, 7 scalars, L=16). That is
irrelevant to the STS ratio, which is set by the hyperbolic/parabolic timescale ratio ~ η/(c·dx)
and does not involve radiation grouping, dust scalars or box size. It would matter for any claim
about the flux-retention observable, and none is made from this probe.*
