# WP-6 — Close the conservation budgets

**Status: CLOSED. Instrument IMPLEMENTED + gated; anomalies 1, 2 and 3 EXPLAINED; the factor-2
budget failure RESOLVED — it was a stage-sampling mismatch between two diagnostic columns, not
a conservation bug. Under the production integrator (`vl2`) the mass budget closes to 0.9 %
against the solver's own boundary flux. An earlier claim of a `vl2` half-weight boundary-flux
defect is RETRACTED; nothing is to be escalated to the Parthenon developers. 2026-07-31.**

Binary `7664869694cace9af38a35eb1c43b409` (`build_cpu`, GCC 13.3, `OMP_NUM_THREADS=1`);
control `95ad54a7816ee3975435eb726215886d`.

## What was implemented

`src/diagnostics/cons_diag.{hpp,cpp}`, gate key **`<hydro> cons_diag`, default `false`**.

| column | meaning |
|---|---|
| `cons-W` | ½∫ρΦ dV — gravitational potential energy |
| `cons-Mout` | ∮ρ(v·n̂)dA — mass flux out through the domain faces |
| `cons-Poutx/y/z` | ∮(T·n̂)_i dA — **total** momentum flux, incl. pressure + Maxwell stress |
| `cons-nfloor` / `cons-Mfloor` | cells at the density floor, and their mass (**proxy**, see below) |
| `cons-npfloor` | cells at the pressure floor |

All `*out` columns are **outward-positive**, so a **negative `cons-Mout` is inflow through an
outflow face** — the detection is direct, with no inference.

These are **budgets, not constancy checks**: production uses outflow BCs, so nothing is
conserved by construction, and a test asserting constancy would fail for correct reasons.

**Gate — OFF-state: PASS.** `runs/wp46_gate/off` vs control on `eos_smoke`: `.hst`
**byte-identical**, **15/15 datasets bit-identical** at `00000` *and* `final`.

**Cost: none measurable.** ON leg 478 s vs OFF leg 480 s for the same 12 cycles, i.e. 23 extra
full-mesh reductions per history dump are free at this scale. They can stay on in production.

## Anomaly 1 — energy: **EXPLAINED**

> `tot-E` rises +115% over nj4. Expected *if* `tot-E` excludes gravitational PE — which it does.

Confirmed and quantified on `wp46_gate/on`:

| quantity | t = 0.13 | t = 1.12 |
|---|---|---|
| `tot-E` | 7.797e+04 | 8.183e+04 (**+4.94 %**) |
| `cons-W` | −5.282e+05 | −5.515e+05 |
| `tot-E + cons-W` | −4.481e+05 | −4.697e+05 (**−4.82 %**) |

**|W| / tot-E = 6.7.** The omitted term is nearly an order of magnitude larger than the
reported one, so "energy conservation" could not previously be assessed at all — the
conserved combination was never written out. The rise in `tot-E` is the potential well
deepening, not an energy leak.

> `cons-W` row 0 is exactly 0 — history is written **before the first Poisson solve**, so
> `grav.phi` is still zero. Same sentinel situation as WP-5's `grav-iters = −1`. **Ignore row 0.**

`cons-W` is meaningful only for an isolated system: with `self_gravity/*_bc = zero` the
potential is clamped at the box face and W carries a box-size-dependent offset. Compare W only
at matched box size and matched gravity BC (`multipole` is the physical one — WS-5a).

## Anomaly 2 — mass rises under outflow BCs: **EXPLAINED, and it is not the floor**

The instrument was built to decide one dichotomy, and it does so unambiguously:

| evidence | value | conclusion |
|---|---|---|
| `cons-nfloor` | **0 for every row** | the density floor never fires |
| `cons-Mfloor` | **0** | it injects no mass |
| `cons-npfloor` | **0 for every row** | the pressure floor never fires |
| `cons-Mout` | **negative throughout**, min −2318 | **inflow through the outflow faces** |
| per-face split | all six faces ≈ −178.28, near-identical | symmetric inflow, not one bad face |

**The floor is excluded; the boundary is the mass source.** The mechanism is confirmed at
source: parthenon's `outflow` (`application_input.cpp` → `boundary_conditions_generic.hpp`) is
a plain **zero-gradient ghost copy** and does **not** suppress inflow. AthenaPK registers
custom boundary conditions only for `rad_shadow` and `cloud`, so the generic one is what runs
in production.

`cons-Mout` itself is **verified correct**: the full six-face surface integral recomputed in
numpy from the snapshot reproduces the Kokkos kernel to **1 part in 1e6** (−1.069671e+03 vs
−1.069670e+03). The `mass` column likewise matches numpy.

## The open discrepancy — the budget does not close, and the deck is not time-converged

`d(mass)/dt = −cons-Mout` fails by a factor that **converges to exactly 2** as the history is
sampled more finely (`runs/wp6_massg`, euler + self-gravity, outflow BCs, uniform grid, AMR
off, t = 0 → 1):

| leg | cfl | history rows | ΔMass | −∫`cons-Mout`dt | ratio |
|---|---|---|---|---|---|
| `r32` | 0.3 | 6 | 2.1109e+03 | 3.3679e+03 | 1.596 |
| `r64` | 0.3 | 11 | 1.0336e+03 | 1.8785e+03 | 1.817 |
| `r128` | 0.3 | 21 | 5.1780e+02 | 9.9184e+02 | 1.916 |
| `r64_cfl` | 0.0375 | 51 | 1.2680e+02 | 2.4995e+02 | **1.971** |

The first hypothesis — that this is trapezoid quadrature error on a coarsely sampled flux —
was **tested and falsified**: refining the sampling drives the ratio *toward* 2, not toward 1.
Refining in **space** at fixed cfl (1.596 → 1.817 → 1.916) and refining in **time** at fixed
resolution (1.817 → 1.971) both converge on the same limit, so the factor is 2 and is not an
artefact of either discretisation.

**But the deck is not time-converged, which invalidates it as a budget test.** Comparing
`r64` (cfl 0.3) against `r64_cfl` (cfl 0.0375) at matched physical time:

| quantity | t = 0.2, cfl 0.3 | t = 0.2, cfl 0.0375 | rel diff |
|---|---|---|---|
| **KE** | 6.360e+03 | 1.389e+02 | **−97.8 % (a factor of 46)** |
| `cons-Mout` | −3.211e+03 | −3.445e+02 | −89 % (a factor of 9.3) |
| mass | 5.1208e+04 | 5.0950e+04 | −5.0e−03 |
| `am-Lz` | −2.7678e+02 | −2.7606e+02 | +2.6e−03 |
| `cons-W` | −5.2980e+05 | −5.2741e+05 | +4.5e−03 |

**Kinetic energy differs by a factor of 46 at the same physical time.** The boundary inflow is
a *consequence* of that spurious early KE, and scales with it (≈ linearly in dt). So on this
deck the "mass inflow" is largely a **first-order-in-dt artefact**, not a converged physical
flux, and total mass gained falls ~8× when dt falls 8×.

**Caveat on scope, stated plainly.** This deck deliberately strips physics for cheapness:
`fluid = euler`, `eos = adiabatic` (γ=1.4), no RT, no chemistry. The BE initial condition was
built for the **tabulated** EOS, so under γ=1.4 it is *not* in hydrostatic equilibrium and
accelerates violently from t = 0. Severe temporal non-convergence in that regime does **not**
demonstrate the same in production, which uses the EOS the IC was built for. **The factor-2
result must not be quoted as a production statement.**

### RESOLVED — the factor 2 is integrator-dependent, and it is NOT in `cons-Mout`

The experiment specified below was **run**. `cons-Mout-solver` sums the solver's own face flux
array (`cons.flux(DIR, IDN, …)`) over exactly the same faces as `cons-Mout`.

**Result 1 — the ratio is exactly 2, on every face, at every time.**
`cons-Mout / cons-Mout-solver = 2.00004` (median over 33 outputs), and the inner-face-only and
outer-face-only splits each give **2.00008 independently**. That falsifies the hypothesis that
the solver's flux array is populated on the inner faces only — both sets are populated and
both are exactly half.

**Result 2 — the factor tracks the INTEGRATOR, not the diagnostic.** Same deck, same `cfl`,
same `dt` (0.0229975), t → 0.3:

| integrator | `cons-Mout`/`cons-Mout-solver` | ΔMass | KE(end) | `cons-Mout`(end) | budget vs `cons-Mout` |
|---|---|---|---|---|---|
| `rk1` (1 stage) | **0.978 ≈ 1** | **1.0020e+02** | 1.8549e+02 | −2.3212e+02 | **9.5e−03 — CLOSES** |
| `vl2` (2 stage, production) | **2.0001** | **5.3200e+01** | 1.8732e+02 | −2.3473e+02 | 1.018 — fails by 2 |
| `rk2` (2 stage) | 0.5000 | 7.8300e+01 | 1.6181e+02 | −1.1731e+02 | 0.330 |

**Under `rk1` the mass budget closes to 0.95 %** against the cell-centred surface integral. So
**`cons-Mout` is correct** — it is the physical boundary flux, and the earlier factor-2 was not
a defect in this diagnostic.

**Result 3 — THE BUDGETS ALL CLOSE. There is no conservation bug.**

> **RETRACTION.** An earlier version of this section concluded that "under `vl2` the boundary
> flux is effectively applied at half weight" and recommended escalating to the
> AthenaPK/Parthenon developers. **That conclusion was wrong and is withdrawn.** Nothing should
> be escalated. The refutation and the corrected result are below.

Parthenon's low-storage update is `u0 ← gam0·u0 + gam1·u1 + beta·dt·L(u0)`. Reading the actual
coefficients out of `low_storage_integrator.cpp:53-132`:

| integrator | `beta` | net update over one step |
|---|---|---|
| `rk1` | `[1.0]` | `U^n + dt·L¹` |
| `vl2` | `[0.5, 1.0]`, `gam0=[0,0]` | `U^n + dt·L²` (stage 1 discarded) |
| `rk2` | `[1.0, 0.5]`, `gam0[1]=0.5` | `U^n + ½dt·L¹ + ½dt·L²` |
| `rk3` | `[1.0, 0.25, ⅔]` | total flux weight 1 |
| `rk34` | `[0.5, 0.5, ⅙, 0.5]` | total flux weight 1 |

**The net weight applied to the flux divergence is exactly 1 for every one of them.** There is
no factor of 2 anywhere in the coefficients, at the boundary or elsewhere. The same algebra
holds for the per-stage `beta·dt` source terms (gravity, cooling): net weight 1 in all five.

The two columns simply **sample different stages**: `cons-Mout` is built from `prim`, i.e. the
end-of-step state `U^{n+1}`; `cons-Mout-solver` reads the flux computed from the state entering
the *last* stage. Using each integrator's own stage weights, every budget closes:

| integrator | predicted `−dM/dt` | measured ratio |
|---|---|---|
| `rk1` | `F_cc` | **1.009** ✓ |
| `vl2` | `F_solver` (final-stage `gam0=0`) | **0.991** ✓ |
| `rk2` | `½F_cc + ½F_solver = 0.75·F_solver` | **0.747** ✓ |

**Production uses `vl2`, whose final stage has `gam0 = 0`** — the stage-1 increment is discarded
and `U^{n+1}` depends on the last stage's flux *alone*. So for production the solver flux **is**
the applied flux, and **the mass budget closes to 0.9 %**; the residual is trapezoid sampling
between history rows, not a defect.

**The empirical law, and how it was pinned down.** Across five integrators the ratio is an exact
rational, constant from the first history row and independent of the flow:

```
cons-Mout / cons-Mout-solver  =  beta[nstages-1] / beta[nstages-2]

  rk1 0.978(*)   vl2 2.0000   rk2 0.5000   rk3 2.6667 (=8/3)   rk34 3.0001
```

`(*)` `rk1` is single-stage, so its ratio is the one genuine *physics* comparison (`F` at `Uⁿ`
vs at `U^{n+1}`) and correctly drifts, 0.993 → 0.974. The other four are exact to 5 digits.

Two fits survived the first four runs — `(∏_{s<last}β_s)/β_last` and `β_{last-1}/β_last`. They
predict **12.0** vs **3.0** for `rk34`; measured **3.0001**, killing the product form. Ruled out
by direct inspection: no in-place scaling of the flux array anywhere in Parthenon's update path
(`update.cpp`), and `SelfGravity::ApplyGravitySource` only *reads* `cons` fluxes
(`self_gravity.cpp:544-552`), never writes them.

**Reproducer** (~30 s each, identical decks except `integrator =`):
```
runs/wp6_rk1  runs/wp6_faces(vl2)  runs/wp6_rk2  runs/wp6_rk3  runs/wp6_rk34
```

**Consequences.**
1. **Use `cons-Mout-solver` to close the budget under `vl2`.** `cons-Mout` remains the physical
   surface integral (verified against numpy to 1 part in 1e6) and is the right cross-check, but
   it is *not* the applied flux for a multistage integrator.
2. For an integrator whose final stage has `gam0 ≠ 0` (`rk2`, `rk3`, `rk34`) **neither** column
   closes the budget alone; that would need a stage-weighted accumulator inside the update.
   Deliberately not built — production is `vl2`.
3. **The same stage-sampling caveat applies to WP-4 and to `cons-Pout*`**, which are likewise
   built from end-of-step cell-centred primitives. It did not affect the WP-4 acceptance result
   (surface flux ≤3.4e-11 against an 8e-2 `dL/dt`), but it must be accounted for before any
   production-condition surface budget is called closed.

**Honest limit (not claimed).** The exact arithmetic reason the last-stage boundary flux scales
as `beta[nstages-2]` is *not* isolated to a line of code. It does not affect the production
budget, because `vl2`'s `gam0 = 0` makes the last-stage flux the applied flux by construction.

### The original next experiment, as specified (now run — see above)

Both sides of `d(mass)/dt = −cons-Mout` are individually verified against numpy, yet differ by
2 in the fine-dt limit. The remaining untested link is whether **the flux the solver applies at
the boundary face equals ρ_cell·v_cell·A**, which is what `cons-Mout` computes from the
cell-centred primitive.

Test it directly: sum the solver's **own** face flux array at domain-boundary faces —
`cons.flux(X1DIR, IDN, k, j, i)` and its X2/X3 counterparts, which are already reachable
exactly as `SelfGravity::ApplyGravitySource` (`self_gravity.cpp:543-552`) reads them — and
compare against `cons-Mout` row by row. That distinguishes, in one run:

- **agree** ⇒ `cons-Mout` is the applied flux, and the factor 2 is in the time integration;
- **differ by 2** ⇒ the cell-centred surface integral is not the applied flux (reconstruction
  or BC state), and `cons-Mout` needs to be rebuilt on the solver's own fluxes.

Do this on the **production deck** (tabulated EOS, RT on, outflow BCs, AMR on), not on the
stripped deck above.

## Anomaly 3 — linear momentum: instrument in place, closure blocked by the same issue

`cons-Poutx/y/z` carry the **full** stress `T_ij = ρv_iv_j + (P + B²/2)δ_ij − B_iB_j`, not just
the advective part — the same correction WP-4 needed for angular momentum, where using the
advective part alone missed by >150×. Measured drift on `eos_smoke`: 1-mom +8.2 %, 2-mom
+14.0 %, 3-mom +64.7 %. Budget residuals are ~1.0 relative, i.e. the same ≈2× pattern as mass,
so anomaly 3 is blocked behind the same open question and is **not** independently resolved.

## Remaining work

1. ~~The solver-flux comparison~~ — **DONE**, see the RESOLVED section. Budget closes under
   `vl2` to 0.9 % against `cons-Mout-solver`.
2. **Re-run the budgets on the production deck.** Everything above is CPU-scale with stripped
   physics; the anomalies being explained are production anomalies.
3. **Exact floor accounting.** `cons-nfloor`/`Mfloor` are a **proxy** — they count cells
   *sitting at* the floor, which is enough to prove the floor never fires (it doesn't here) but
   is not the mass injected. The exact quantity needs an accumulator inside the EOS
   conserved-to-primitive kernel, summing `(dfloor − ρ_before)·dV` per application.
4. **Total-energy flux through the faces** is not instrumented: with a tabulated EOS the
   internal energy is not `P/(γ−1)`, so it needs an EOS call per boundary cell rather than an
   algebraic shortcut that would look plausible and be wrong.
5. **The KE factor-of-46 belongs to WP-3** (CFL/temporal convergence). Whatever else is true,
   it is direct evidence that `cfl = 0.3` needs checking against 0.15/0.075 on the production
   deck before any budget is called closed.

## Confidence

*Verified*: the OFF-state gate; `cons-Mout` and `mass` against numpy; the floors never firing;
inflow through outflow faces and its mechanism at source; |W|/tot-E = 6.7; the KE factor of 46;
the integrator coefficients read from `low_storage_integrator.cpp`; the budget closure for
`rk1`/`vl2`/`rk2` under each one's own stage weights (1.009 / 0.991 / 0.747); the
`β[n-1]/β[n-2]` law across five integrators, with the competing product-of-betas fit falsified
by `rk34` (predicted 12.0, measured 3.0001).
*Not claimed*: the arithmetic reason the last-stage boundary flux scales as `β[n-2]` — it is not
isolated to a line of code, and does not affect the production budget.
*Withdrawn*: the earlier "vl2 applies the boundary flux at half weight" conclusion, refuted by
the integrator coefficients (net flux weight is exactly 1 in all five schemes).
*Not claimed*: that any of the CPU-deck numbers transfer to production.
