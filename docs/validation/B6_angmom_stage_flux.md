# B6 — stage-consistent angular-momentum surface flux (`am-FTsolver*`)

**Status: CLOSED. The original gate's prediction is falsified twice over. The angular-momentum
surface diagnostic is VALIDATED — the budget closes to ~1 % — but `am-FTsolver*` and `am-FT*` are
equivalent to 0.02 %, so B6 changes nothing.** Re-gate: job 2453283.

## Bottom line (read this first)

1. **The budget closes.** On a deck that carries **100 % of L_z** out through the faces with
   T_grav = T_mag ≡ 0, `d(L_z)/dt = −(surface flux)` holds to **0.5–1.0 %** with either column.
   That validates the angular-momentum surface diagnostic itself — which the FHC deck **cannot**
   do, because B7's discretization drift swamps the flux there by ~10⁴.
2. **The two columns are equivalent**: `am-FTz/am-FTsolverz = 1.0002` under *both* integrators.
3. **B6 therefore buys nothing measurable — and costs nothing.** `am-FTsolver*` stays as the more
   principled term (it reads the flux the update actually applied), but no result moves.

## The defect B6 addresses

`am-FT*` is built from **end-of-step primitives**, so under a multistage integrator it samples
u^{n+1} while the flux that was actually applied came from the state entering the last stage.
`am-FTsolver*` instead reads the solver's own momentum fluxes, `cons.flux(dir, IM1..IM3)`, which
already carry the pressure and Maxwell stress (`angmom_diag.cpp:153-186`). By construction it
cannot drift from what the update applied.

## What the first gate predicted, and what it got

Gate `gate_b6.sh`, job 2452597, deck `runs/b2b4_gate/fhc.in`, `build_cpu` `0c81219a`.

Prediction, by analogy with WP-6's mass result: `am-FT*/am-FTsolver*` ≈ **2.0** under `vl2`
(= `beta[nstages-1]/beta[nstages-2]`) and ≈ **1.0** under `rk1`.

| integrator | `am-FTz/am-FTsolverz` | predicted | mass ratio, same run |
|---|---|---|---|
| vl2 | **1.098** (range 1.03 … 2.0e8) | 2.0 | **2.0039** (2.0000–2.0078, std/med 0.001) |
| rk1 | **1.217** (range −24.7 … 9.7e8) | 1.0 | 0.928 |

The mass column reproduces WP-6 **exactly** — 2.0000–2.0078 with std/med = 0.001 — which confirms
the deck, the integrator and the WP-6 mechanism are all as expected. The angular-momentum ratio
does not.

**The `rk1` leg is decisive.** With a single stage there is no earlier stage to mis-sample, so a
pure stage-sampling difference *must* give exactly 1.0. It gives 1.217 with a range spanning four
sign changes and eight orders of magnitude. Whatever separates `am-FT` from `am-FTsolver` for
angular momentum, it is **not** stage sampling.

The vl2 series also decays monotonically — 3.15, 1.52, 1.28, 1.17, 1.11 — i.e. it is a transient,
not the time-constant signature WP-6 found for mass ("constant to five digits, independent of the
flow"). The reported median 1.098 was a snapshot of that decay over only 12 history rows.

**Mechanism (inferred):** `am-FT*` evaluates the stress at the **cell centre** and multiplies by
the face area; `am-FTsolver*` uses the actual **face** value from the Riemann solver. For angular
momentum the leading term is the pressure torque r × P n̂, which very nearly cancels over a cube,
so the surviving integral is dominated by that centring difference. Mass has no such
cancellation, which is why only mass shows the clean stage factor. This is the same
cancellation-sensitivity that produced the false-positive FAIL in the `build_gpu_v4` gate.

## Budget closure could not be gated on that deck either

Re-analysing the same history files (`b6_budget.py`) for closure of
d(L)/dt = −(surface flux) + T_grav + T_mag gave residual/scale ≈ 1.0 for **both** columns —
i.e. no closure at all. The magnitudes say why:

| term (vl2, z-axis) | median abs |
|---|---|
| **dL_z/dt** | **6.68** |
| am-FTz | 4.1e-02 |
| am-FTsolverz | 3.7e-02 |
| am-Tgravz | 1.3e-01 |
| am-Tmagz | 4.2e-02 |

Angular momentum is changing **50–180× faster than any accounted term**. That is **B7** —
Cartesian finite volume does not conserve angular momentum — and here it swamps the surface flux
entirely (L_x drifts 173 → 183, **+5.8 %**, consistent with B7's 7.9–11.7 % on its smoke deck).

In that regime the choice between `am-FT` and `am-FTsolver` is both unmeasurable *and*
irrelevant: both are ~1e-4 of dL/dt. **On the FHC deck, B6 cannot matter.**

## The re-gate (job 2453198)

Gating requires a deck where the surface angular-momentum flux is the **dominant** term:

- **No gravity, no B** ⇒ T_grav = T_mag ≡ 0, so the budget is exactly d(L_z)/dt = −(flux) and
  any residual is discretization with nothing to hide behind.
- **A dense blob offset in y** (radius 2 at c_y = 4), drifting at v_x = 1 and leaving through the
  +x face. It carries L_z = −V·M_blob·c_y out of the box, so the flux is large and does **not**
  cancel.
- A uniform wind will **not** work — the algebra was checked and every face pair cancels exactly
  (±x contribute ±ρUV·Vol, ±y cancel them, ±z vanish), because uniform translation conserves L.
  A wind deck would have gated 0 = 0. This is the same class of null test as WP-17 attempt 1 and
  the first B1 blast falsification.

Acceptance: whichever column makes |dL_z/dt + flux| / |dL_z/dt| smaller wins.
If `am-FTsolverz` wins under vl2 by a wide margin and the margin **collapses** under rk1, B6 is a
staging fix as designed. If it wins by the **same** margin under both, B6 is a **centring** fix
and the wording everywhere must change. If neither wins, `am-FTsolver*` buys nothing measurable.

## The re-gate result (job 2453283)

525 cycles per leg to t = 12.0, both integrators. **L_z: −120.1330 → −0.0003** — the blob and all
its angular momentum left the box, exactly as designed.

Judged over the rows where the flux is actually large (|dL_z/dt| > 5 % of peak):

| | `−dL_z/dt ÷ am-FTz` | `−dL_z/dt ÷ am-FTsolverz` | `am-FTz/am-FTsolverz` |
|---|---:|---:|---:|
| **vl2** (2 stages) | 0.99997 | 1.00005 | **1.00025** |
| **rk1** (1 stage) | 1.00004 | 1.00028 | **1.00020** |

The original prediction was 2.0 under vl2 and 1.0 under rk1. Measured: **1.0002 under both.**
Falsified a second time, now on a deck purpose-built to measure it.

### A metric error of my own, corrected

My first pass at this verdict reported "**am-FTz wins by 2.17×**". That was an **artifact of the
metric, not a result.** It took the median of |dL/dt + flux| / |dL/dt| over all 121 history rows —
but the blob has fully left by t ≈ 10, after which dL/dt ≈ 0 and the normalised residual is 0/0
noise that dominates the median. Restricted to rows where the flux is significant, the two columns
are indistinguishable. Normalising by a quantity that passes through zero manufactures
differences — the same cancellation trap that produced the false-positive FAIL in the
`build_gpu_v4` gate.

## Standing conclusion regardless of the re-gate

`am-FTsolver*` reads the flux the solver actually applied, so it is the correct term to put in an
angular-momentum budget on principle. But the claimed *mechanism* (stage sampling, factor 2.0) is
falsified, and on the production deck the term is four orders of magnitude below the B7
discretization drift, so **B6 changes no production result.**
