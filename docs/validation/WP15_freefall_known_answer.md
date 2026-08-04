# WP-15 — collapse against a known answer: pressure-free free-fall

**Status: PASS.** AthenaPK reproduces the exact pressure-free collapse solution, with a residual
that **converges at first order in dx (p = 1.02)** and is fully attributable to representing a
spherical density discontinuity on a Cartesian grid.

The headline number arrived with the wrong sign, and getting from there to here required
falsifying three candidate causes. That process is the substance of this document.

---

## 1. Why this WP exists

**Every other convergence study in this campaign is self-referential.** WP-7 refines the grid,
WP-3 the timestep, WP-8 the Jeans number — each asks whether the code agrees with *itself*. None
can detect an error the code makes consistently at every resolution. WP-14 and WP-16 do compare
against analytic solutions, but for the base MHD update and the ambipolar operator **in
isolation** — not for gravitational collapse, which is what the flagship result is about.

This is the only test in the campaign that compares *collapse* against an external answer.

## 2. The known answer

A pressure-free uniform sphere collapses homologously with an exact Newtonian solution:

```
r/r0 = cos²β      ρ/ρ0 = sec⁶β      t/t_ff = (2/π)(β + sinβ cosβ)
t_ff = sqrt(3π / (32 G ρ0))
```

With `four_pi_G = 1` (so G = 1/4π) and ρ0 = 1: **t_ff = π√(3/8) = 1.923825**. Arrival times:
t(2×) = 1.073029, t(10×) = 1.616643, t(100×) = 1.836188. These are properties of the equations,
not of AthenaPK.

The measured quantity depends on the Poisson solve, the gravity source, the hydro advection **and
their stage coupling** all being right simultaneously — precisely the combination WP-13's
stale-RHS bug corrupted, and which no self-convergence study could have caught, because a
one-stage-lagged RHS converges beautifully to the *wrong* answer.

## 3. The finding, and the measurement that made it legible

The first reading used *crossing times* — when does ρ_max reach a given contrast — and gave
**negative** errors: −0.909 % (2×), −0.430 % (10×), −0.532 % (100 %). The collapse was running
*ahead* of free-fall. This document's own acceptance criterion, written before the run, flagged
negative as the serious direction: *"nothing physical makes a pressurised sphere collapse FASTER
than free-fall."*

Crossing times were the wrong instrument. They require interpolating between dumps, and they
compress the whole history into three numbers. Comparing **density at each dump time** — no
interpolation, no fitting — shows the actual structure:

| t | ρ error (64³) |
|---|---|
| 0.05 | +0.03 % |
| 0.50 | +0.64 % |
| 1.00 | +1.18 % |
| 1.55 | +3.12 % |
| 1.80 | +8.19 % |

Monotone, positive, growing from t = 0. That single change turned three ambiguous numbers into a
curve with a shape, and every subsequent step depended on it.

## 4. Three candidate causes, two falsified by measurement

### (a) The ambient is not massless — EXCLUDED

`poisson_test` gives interior and exterior the same *pressure*, so a low-density ambient is
**hot**: at ρ_out = 1e-6 its sound speed is 1.183 versus the interior's 0.001183, and it destroyed
an early attempt outright (dt fell 4.2e-2 → 6.4e-7, KE jumped 3.5 → 1305, every cell to dfloor).
Raising ρ_out to 1e-2 fixed that but gave the ambient 14.3 % of the sphere's mass — and infalling
ambient augments the enclosed mass, which shortens t_ff as 1/√M. Right size, right sign.

Swept across **3.3× in ambient mass**:

| t | ρ_out = 1e-2 | 5e-3 | 3e-3 |
|---|---|---|---|
| 0.101 | +0.123 % | +0.123 % | +0.123 % |
| 1.001 | +1.178 % | +1.175 % | +1.173 % |
| 1.801 | +8.190 % | +8.183 % | +8.179 % |

Identical to three decimals. **Excluded.**

### (b) Pressure — EXCLUDED

Halving the pressure (1e-6 → 5e-7) moved the error by < 0.01 % at every time. Pressure is also
excluded on physical grounds: it can only **delay** collapse, so it cannot produce a positive
density offset in any amount.

### (c) The discretized sphere carries the wrong mass — EXCLUDED, and this one is instructive

A stair-stepped sphere does not have the ideal sphere's mass. Measured directly from the ICs:

| leg | dx | mass excess vs ideal |
|---|---|---|
| 64³ | 0.2500 | **+0.4984 %** |
| 128³ | 0.1250 | **+0.0743 %** |

+0.5 % excess mass gives ≈ −0.25 % on t_ff, the right size. But it converges at **p ≈ 2.75**,
which predicts the density error should fall **6.7×** from 64³ to 128³. It falls **2.0×**. The
rates do not match, so the IC mass excess — while real — is not what produces the error.

## 5. The answer: surface discretization, p = 1.02

| t | 64³ | 128³ | ratio |
|---|---|---|---|
| 0.501 | +0.637 % | +0.325 % | 1.96× |
| 0.551 | +0.625 % | +0.302 % | 2.07× |
| 0.601 | +0.715 % | +0.341 % | 2.10× |
| 0.651 | +0.672 % | +0.401 % | 1.68× |
| 0.701 | +0.890 % | +0.412 % | 2.16× |

Mean ratio over the later half: **2.03× for a 2× refinement ⇒ p = 1.02.**

**This is the result.** A resolution-*independent* offset (p ≈ 0) would have meant the
gravitational source term is genuinely too strong — a first-order defect in the flagship physics.
p = 1.02 means the opposite: the code **converges to the exact solution**, and the residual is a
property of the discretization, not of the gravity.

First order is the signature of a **geometric surface error**. A sphere on a Cartesian grid has a
stair-stepped boundary whose misrepresented volume fraction scales as dx; the gravitational field
of that body differs from the ideal sphere's at O(dx), and the collapse rate follows.

**No contradiction with WP-14's p ≈ 2.06–2.13.** That measured the base MHD update on a *smooth*
solution. A density discontinuity is not smooth, and first order at a discontinuity is correct
behaviour for a shock-capturing scheme — the same argument WP-16 makes for C-shocks.

## 6. Scope and honest limits

- **Hydro + self-gravity only.** No MHD, no radiation, no chemistry, no AMR. This validates the
  gravity–hydro coupling, which is what the free-fall solution can speak to; it says nothing about
  the non-ideal terms (WP-16) or the base MHD order (WP-14).
- **The 128³ leg was stopped at t ≈ 0.70** by its wall limit (500 cycles at dt_ceil = 0.002 needs
  ~3.2 h at 128³). The ratio is stable over t = 0.4–0.7 and the order is read there; the deeper
  points exist only at 64³.
- **p = 1.02 is from a two-point ladder.** A third resolution would tighten it. The conclusion
  (converges, ≈ first order, not resolution-independent) is robust to that; the exponent's second
  digit is not.
- The error reaches +8 % in density by t = 1.80, i.e. 94 % of t_ff. Close to the singularity a
  uniform grid cannot represent the core at all, and nothing there should be read as physics.

## 7. Six setup faults, none of them code defects

WP-15 reused `runs/validation_ws5a/sphere_multipole.in`, a deck written for `nlim = 1`. It carried
no sane value for anything that matters over many cycles:

1. no `<parthenon/output*>` hst block — the first run produced no history at all
2. `output0/dt = 1e-10` — 3214 dumps, 570 MB, run became pure I/O
3. `cfl = 1.0e-3` — 300× below production
4. ρ_out = 1e-6 ⇒ ambient sound speed 1000× the interior's ⇒ blow-up
5. **no gravity-based dt limiter exists**, so a pressure-free problem has no dt constraint at all:
   dt came out equal to tlim and the whole run took one step
6. `dt_max` is an **error trigger** in Parthenon, not a clamp — the clamp is `dt_ceil`
   (`driver.hpp:98`); `dt_max`/`dt_min` throw

Point 5 is worth keeping: the very property that makes the analytic solution exact (negligible
pressure) removes the only thing constraining the timestep. Production never sees it because it
has a real EOS and real pressure.

**Generalised rule, earned four times over this session:** when adopting a deck built for a
different purpose, audit *every* numerical setting it carries — output cadence, CFL, tlim, floors,
ambient state — not just the physics block.

> *Confidence:* every number **measured** this session. The attribution to surface discretization
> is **inferred** from p = 1.02 plus the exclusion of the alternatives — consistent with all the
> data, but the mechanism was not isolated by, say, comparing against a grid-aligned cube.
