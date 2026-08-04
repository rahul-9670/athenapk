# B7 — Cartesian-FV angular-momentum non-conservation, CLOSED (on the one axis it can be)

**Status: CLOSED with explicit scope.** The budget closure is **resolved on L_y**, where the
numerical residual is **2.5× larger than the physical change itself**. No closure is claimed on
L_x or L_z — their `dL` is unresolvable at the history file's printed precision, which is a
*separate* failure mode from the one the test was designed around.

Job 2454557, `build_gpu_v4` (`869c1d34`), 4 GPUs, restart `deep_amr/run/parthenon.out2.00001.rhdf`
(t = 1.0145874, cycle 250, 792 blocks), production configuration, **history written every cycle**.

## What was already known, and why it was not enough

The deep AMR production run measured the *drift* under production conditions — L_x **+2.18 %**,
L_y **−9.79 %**, L_z **+2.44 %** over a full t_ff — and that stands. What it could not do is close

```
d(L_i)/dt = -(surface flux)_i + T_grav,i + T_mag,i + R_i
```

because its history spacing was **1.78e-03**, far coarser than the surface flux varies. A
trapezoid over those rows is not a valid integral, so the >100 % "unexplained" fraction computed
from it was a **quadrature artefact** and was deliberately withheld rather than reported as a B7
confirmation.

## The fix was cadence

Same restart, same configuration, history every cycle: **120 rows at median spacing 1.00e-05**,
i.e. **178× finer** than before. Blocks grew 792 → 939 as AMR refined during the window.

Getting per-cycle output required working around a Parthenon constraint: it **throws** if both
`dt` and `dn` are enabled on one output block, so `dt` must first be disabled with a negative
value — `output0/dt=-1.0 output0/dn=1`. That same fault silently gave WP-1 round 1 zero cycles on
three legs while SLURM still reported `State=COMPLETED, ExitCode 0:0`. The script verifies rows
track cycles rather than assuming (`120 rows / 120 cycles` — confirmed).

## Result

| axis | actual dL | accounted (flux+T_grav+T_mag) | residual R | R/\|dL\| | dL in quanta | usable |
|---|---:|---:|---:|---:|---:|---|
| x | −3.00000e-03 | −5.43448e-03 | +2.43448e-03 | +81.1 % | **3** | **no** |
| **y** | **+3.10000e-03** | **+1.07850e-02** | **−7.68500e-03** | **−247.9 %** | **31** | **yes** |
| z | −1.50000e-02 | −8.73438e-04 | −1.41266e-02 | −94.2 % | 15 | **no** |

**On L_y:** the physical budget predicts a rise of **1.079e-02**; the actual rise is
**3.10e-03**. The missing **7.69e-03** is Cartesian finite-volume non-conservation — **2.5× the
size of the change it is competing with** over this window. That is B7 in quantitative form.

## Two independent validity checks, both required

**1. Quadrature convergence** — recompute the residual at half the sample rate; if the cadence is
adequate the answer must not move:

```
R_x: full +2.43448e-03   half-rate +2.35390e-03   change 3.31 %   OK
R_y: full -7.68500e-03   half-rate -7.55452e-03   change 1.70 %   OK
R_z: full -1.41266e-02   half-rate -1.30943e-02   change 7.31 %   OK
```

Converged on all three. This is precisely the test the previous attempt would have failed.

**2. `dL` resolvability** — and this one rules out two axes that check 1 passed. The history file
prints **6 significant figures**, so a small `dL` riding on a large `|L|` is quantized:

| axis | \|L\| | printing quantum | dL | uncertainty from quantization alone |
|---|---:|---:|---:|---:|
| x | 157.5 | 1e-03 | 3 quanta | **±33 %** |
| y | 20.1 | 1e-04 | 31 quanta | ±3.2 % |
| z | 107.5 | 1e-03 | 15 quanta | ±6.7 % |

A converged flux integral does not rescue an unresolvable `dL`. The first version of this analysis
checked only the quadrature and would have reported all three axes — including x, whose result is
meaningless at ±33 %. **These are independent failure modes and both must pass.**

L_x can never be resolved this way regardless of run length: |L_x| is too large relative to its
change. L_z would need roughly twice the span. Recovering either would mean computing L directly
from `.phdf` dumps at full double precision rather than from the history file — a viable route,
not taken here because L_y already answers the question.

## No extrapolation

The per-t_ff rate is deliberately **not** quoted. This window spans 1.4e-03 in t; naively scaling
its rate would predict a far larger drift than the **−9.79 %** actually measured for L_y over the
full run, so the non-conservation rate is plainly **not constant in time**. The residual is a
statement about this window.

## Run termination

The job ended at 120 of a requested 500 cycles with `exit=134` — a GPU OOM
(`failed to allocate 1.068 MiB, label="bnd_flux::rad.Fr2_g1"`), not a wall-clock truncation. The
cause is **block-distribution imbalance**, not capacity:

```
device 0: 49.3 GiB    device 1: 56.6 GiB    device 2: 79.2 GiB (saturated)    device 5: 0.0 GiB (idle)
```

One rank carries ~60 % more than another while a fourth GPU goes unused. The same signature killed
the WP-1 `cap=150` leg twice. It did not affect this result — the 120 rows collected before the
crash already satisfy both validity checks — but it is a real scaling defect on deep AMR
hierarchies and is worth its own investigation.
