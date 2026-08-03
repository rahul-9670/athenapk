# WP-3 — temporal (CFL) convergence of the production configuration

**Status: PASS.** Production's `cfl = 0.3` is converged in time. Quartering it moves the
flux-retention observable by **+0.0737 %** — 0.46 % of WP-18's σ ≈ 16 % — and Richardson puts the
residual at `cfl = 0.075` at **+0.0196 %**.

- Observed order in time **p = 1.13**, which is the expected value, not a shortfall (§3).
- `mag-Jsq` moves **+1.06 %** — an order of magnitude more than everything else, consistent with
  every other study in this campaign finding it numerics-sensitive. Still not quotable.
- **Cost is the real reason to keep 0.3:** halving `cfl` doubles the cycle count exactly, and
  buys 0.05 % on the answer.

---

## 1. What was missing, and why the earlier attempt did not count

`cfl = 0.3` had never been tested against a smaller value on the production deck. One prior
attempt used a cheap CPU deck and reported **KE differing by a factor of 46** between cfl = 0.3
and 0.0375 — but that deck was out of equilibrium, and a factor of 46 on a configuration nothing
else shares transfers nowhere. It was never a production statement.

## 2. Design

Three legs at 256³ on the **new production physics** (binary `f181c0a1`, diode fluid BCs, relative
gravity residual), differing in `parthenon/time/cfl` and nothing else:

| leg | `cfl` | directory | cycles to t = 1.0 |
|---|---|---|---|
| production | 0.3 | `runs/root_ladder/r256_sw` (**reused**, not re-run) | 179 |
| half | 0.15 | `runs/wp3_cfl/cfl0.15` | 357 |
| quarter | 0.075 | `runs/wp3_cfl/cfl0.075` | 715 |

Reusing `r256_sw` as the reference leg is what makes the study cheap — it already exists as the
re-baseline rung.

**Verified apples-to-apples rather than assumed:** every leg's own job log reports binary
`f181c0a1`, deck `c67eb458`, `NX=256 NRANK=3`. The cycle counts 179 → 357 → 715 are exactly 2× and
4×, and the dt at t = 0.90 is 4.1049e-03 → 2.0509e-03 → 1.0251e-03 — factors of 2.0000 and 2.0007.
The knob reached the code.

Only `parthenon/time/cfl` is varied. `diffusion/cfl` (0.3) and `radiation/cfl` (0.4) stay at
production values: those govern the explicit non-ideal and M1 substeps, and moving them too would
confound "is the time integration converged" with "are the operator-split substep counts
converged", which are separate questions.

**Read at t = 0.90, never the t = 1.0 endpoint** (WP-7: the last 0.01 t₀ is a collapse singularity
on these uniform grids and each leg stalls at a different point).

## 3. Result

| leg | `MEtor/MEpol` @ t = 0.90 | vs production | as % of σ |
|---|---|---|---|
| cfl = 0.3 | 7.641359e-03 | — | — |
| cfl = 0.15 | 7.645221e-03 | +0.0505 % | 0.32 % |
| cfl = 0.075 | 7.646990e-03 | **+0.0737 %** | **0.46 %** |

```
observed order in time p = log2(|d(0.3->0.15)| / |d(0.15->0.075)|) = 1.13
Richardson residual at cfl = 0.075                                = +0.0196 %
```

**PASS.** The shift is monotone, converging, and more than two orders of magnitude below the
acceptance threshold. Production's `cfl = 0.3` is not a source of error worth spending compute on.

Other observables at t = 0.90, quarter vs production: `ME` −0.034 %, `KE` +0.198 %, `mass`
−0.037 %, **`Jsq` +1.058 %**.

### p = 1.13 is the expected value, not a shortfall

WP-14 established the **base MHD update** at p = 2.06–2.13 on a smooth problem in this exact
numerical configuration. The full production system adds self-gravity, M1 radiation, tabulated EOS,
chemistry and the non-ideal RKL2 stack — several of which are **operator-split and therefore 1st
order in time by construction**. A system built from a 2nd-order hyperbolic update and 1st-order
splits must land between 1 and 2, and 1.13 says the splitting error dominates.

That is exactly why WP-14 had to come first: without it, p = 1.13 would be an unattributable
number. With it, the shortfall is *located* — it belongs to the operator splitting, not to the
hydro scheme. Improving it means Strang-splitting the packages, not touching `cfl`.

### `Jsq` again

At +1.058 % it moves ~10× more than any other observable, and it moved 9–13 % in WP-2 and *diverged*
under refinement in WP-7 (99.6 → 264.8 → 752.4). Three independent studies, same conclusion:
`mag-Jsq` is a grid-scale quantity and must not appear in a physical claim.

## 4. Cost — the real argument for keeping 0.3

Cycle count scales exactly with 1/cfl: 179 → 357 → 715. Wall times were 2.94e3 / 3.40e3 / 4.39e3 s,
but **those are not a valid cost ratio** — the legs ran concurrently on different nodes and the
throughput figures show it (1.02e6 vs 1.76e6 zone-cycles/wallsecond between two legs that should be
identical per zone-cycle). The honest cost statement is the cycle count: **halving `cfl` doubles the
work and buys 0.05 % on the answer.**

## 5. Scope

- 256³ **uniform** grid. Production runs AMR, where the timestep is set by the finest level and the
  balance between hydro CFL and the RKL2 substep structure differs. The conclusion — that
  `cfl = 0.3` is far inside the converged regime — should hold, but it is measured on the uniform
  ladder.
- Under RKL2 the diffusion substep structure follows dt through `rkl2_max_dt_ratio`, so halving the
  hydro CFL does change the STS pattern. That is a property of the production configuration and is
  deliberately in scope here.
- The first core only, to t = 1.0. Nothing here licenses second-core densities — WP-2 and WP-22
  both draw that line independently.

> *Confidence:* every number **measured** this session from the three legs' own history files.
> Binary/deck/rank identity **verified** from the job logs. The attribution of p = 1.13 to operator
> splitting is **inferred** — consistent with WP-14 and with which packages are split, but not
> isolated package-by-package.
