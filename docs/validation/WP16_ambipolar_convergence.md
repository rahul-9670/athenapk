# WP-16 — Convergence of the non-ideal operators (ambipolar; Hall in part 3)

> **SCOPE NOTE.** `VALIDATION_PLAN.md`'s WP-16 reads *"Complete the Hall C-shock test […] Hall is
> the least-validated of the three non-ideal terms"*. This document covers the **ambipolar**
> operator, which is what `runs/validation/cshock_ad.in` actually tests and what had to be
> untangled first (§1 — that deck is invalid). The **Hall** half is part 3, job 2449287, measured
> by the same method against the analytic whistler / ion-cyclotron dispersion relation
> (`src/pgen/diffusion.cpp` iprob=60, both helicity branches). WP-16 is not closed until that
> lands.

**Status (ambipolar): PASS.** The ambipolar operator converges at **second order** against a closed-form
solution, `p = 1.993 → 1.999` over a 16× resolution range (32 → 512).

The route to that answer is not the one WP-16 was originally set up to take. The C-shock test the
work package named turned out to be invalid, and diagnosing *why* is itself a result, so both
halves are recorded here.

---

## 1. The C-shock test as written cannot answer WP-16

`runs/validation/cshock_ad.in` initialises two uniform states either side of `x = 0` and asserts
in its own header comment that "AD smooths the discontinuity into a continuous C-type shock […]
Run to steady state; the relaxed profile (v_x, B_y) is compared to the semi-analytic C-shock ODE
(`runs/validation/cshock_ode.py`) and self-converges at 2nd order under resolution refinement."

Three things are wrong with that, all verified this session.

### 1a. The initial condition is not a shock

Evaluating the steady 1-D MHD flux vector (Heaviside–Lorentz, `P_mag = B²/2`, `Bx = 1` on both
sides) on the deck's own left and right states:

| conserved flux | left    | right    | mismatch |
|----------------|---------|----------|----------|
| mass           |  2.0000 |  2.0010  |   0.05 % |
| x-momentum     |  5.0000 | 21.3347  | **327 %** |
| y-momentum     | −1.0000 | −3.0000  | **200 %** |
| induction      |  2.0000 |  2.0010  |   0.05 % |
| energy         | 11.0000 | 43.1331  | **292 %** |

Only the mass and induction fluxes balance. A steady shock requires all five. The states were
evidently picked so that `ρv_x` and `v_x B_y` match and the rest was not checked. What the deck
actually poses is an arbitrary discontinuity, which immediately opens into a full Riemann fan.

### 1b. Nothing sustains an upstream state, so the domain drains

Both x-faces use `outflow`. There is no inflow boundary and no piston, so the fan simply advects
out. Running to `tlim = 400` (10× the original `tlim = 40`, job 2448325) makes this unmistakable —
at `t = 400` **all three resolutions are spatially uniform**:

| N   | v_x (min = max) | B_y (min = max) | ρ (min = max) |
|-----|-----------------|-----------------|---------------|
| 128 | −1.0211         | 1.4651          | 1.6810        |
| 256 | −1.2272         | 1.2594          | 1.5512        |
| 512 | −1.4011         | 1.0794          | 1.4458        |

`max − min < 1e-6` on every grid. The "self-convergence order" measured there is a comparison of
three *constants*, which is why it came out `p = 0.51` (v_x), `−0.03` (B_y), `0.18` (ρ). The three
constants differ by 37 % in `v_x` because the fan leaves the box at a resolution-dependent rate.

The structure amplitude decays monotonically throughout:

```
 t      n=128 Δv_x   n=256 Δv_x   n=512 Δv_x
   0      1.3330       1.3330       1.3330
  40      0.1828       0.2608       0.3552
  80      0.0400       0.0575       0.0766
 160      0.0090       0.0042       0.0036
 320      0.0001       0.0000       0.0000
```

This also retires the original diagnosis. The `tlim = 40` measurement (`p ≈ 0.35` v_x, `0.65` B_y)
was blamed on the three resolutions relaxing at different rates, with the falsifiable prediction
that running far longer would drive `d/dt → 0` and lift the order toward 2. `d/dt` **did** fall to
`1.2e-5` — but because the solution went uniform, not because it relaxed, and the order got
*worse*, not better. The prediction is falsified; the hypothesis was wrong.

Note the second row: at `t = 40` the amplitude is already down to 14–27 % of initial **and is
larger on finer grids** (0.183 / 0.261 / 0.355). Numerical diffusion, not ambipolar diffusion, was
setting the decay rate at the resolutions the original test used.

### 1c. The reference solution does not exist

`runs/validation/cshock_ode.py` — the semi-analytic comparison the header commits to — is not
present in the tree and there is no record of it ever having been. Only `cshock_analyze.py` exists.

**Conclusion.** `cshock_ad.in` exercises the AD term (its earlier qualitative result, that the
smoothed layer thickens with `ambipolar_coeff_code`, is real) but it is not a convergence test and
no order should ever have been quoted from it. It is left in place, with this document as the
pointer, rather than deleted.

---

## 2. What WP-16 actually asks, answered

The question behind the work package is: **does the ambipolar operator achieve the scheme's formal
second order?** That is answered far more cleanly by the existing quantitative eigenmode test than
by any shock, because a shock is a discontinuity and *first* order there is the correct behaviour
of any shock-capturing scheme — a low order at a C-shock would not have implicated the AD term
even if it had been measured properly.

`src/pgen/diffusion.cpp` `iprob = 50` (`inputs/diffusion_ambipolar.in`) initialises a uniform guide
field `B_x = B₀` with `B_y = amp·sin(kx)`, which evolves as a damped Alfvén mode with a
**closed-form amplitude**; the pgen computes `A_pred` and reports
`rel err = |A_meas − A_pred| / |A_pred|`. This is a true error against an analytic solution, not a
self-convergence difference, and the test is periodic so no boundary influence enters.

Parameters: `B₀ = 1`, `ρ₀ = 1`, `k = 2π`, `Q_A = 5` ⇒ `η_A = Q_A B₀² = 5`, `η_A k² = 197`,
`v_A = 1` (overdamped, so the mode decays monotonically), `t_fin = 0.01`, `Q_A` fixed
(`ambipolar_coeff = fixed`), `diffusion/integrator = unsplit` so the measured decay is
unambiguously the AD operator and not RKL2 substepping.

Ladder (job 2448497, `build_cpu`, 1-D, seconds per rung):

| N   | dx       | rel. error | observed p |
|-----|----------|------------|------------|
|  32 | 0.03125  | 6.250e-03  |     —      |
|  64 | 0.01562  | 1.570e-03  |   1.993    |
| 128 | 0.00781  | 3.940e-04  |   1.994    |
| 256 | 0.00391  | 9.870e-05  |   1.997    |
| 512 | 0.00195  | 2.470e-05  |   1.999    |

`p = log₂(e(h)/e(h/2))`, monotone and converging on 2 from below across four successive halvings.

**WP-16 PASS: the ambipolar diffusion operator is second-order accurate.**

### Scope, stated explicitly

- This is the **AD operator in the `unsplit` integrator**, in the linear (`amp = 1e-6`) regime,
  with a *constant* coefficient (`ambipolar_coeff = fixed`). Production runs
  `ambipolar_coeff = ionization_chem`, where `η_A` is a function of ρ, T and B. The spatial
  discretization under test is the same; the coefficient evaluation is not covered here.
- The RKL2 super-time-stepping path used in production is **not** measured by this ladder.
- Complementary, already on record: `iprob = 51` covers a spatially varying `η_A`, and the
  base MHD update's order is WP-14 (`WP14_order_of_accuracy.md`, `p ≈ 2.06–2.13`).

### Remaining gap

There is no validated **shock-structure** test for AD in this tree — a real steady C-shock would
need (i) a downstream state that satisfies all five jump conditions, (ii) a fixed-inflow upstream
boundary or a shock-frame formulation so the structure persists, and (iii) the missing ODE
reference solution. That is a genuine piece of work and it is *not* on the path to the flagship
result: the C-shock regime is not what the FHC collapse runs sit in, and the operator's order is
now established by a cleaner test. Recorded as a known gap, deliberately not closed.
