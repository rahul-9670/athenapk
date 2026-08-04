# WP-17 — Sink-particle accretion rate against an analytic answer

**Status: PASS — `Mdot / Mdot_B = 1.01 ± 0.02` at N = 128 (job 2453375).** The first measured sink
accretion rate in this project; see the RESULT section at the end for the numbers and for the
two-point convergence ladder. The narrative below is kept in order, so the six defects that had to
be cleared first are on record — every one of them produced a run that exited 0.

> Note: the "Attempt 2" section below describes a configuration that was itself superseded twice
> (rho 1 → 1e-3, and racc_cells N/64 → 4). The RESULT section is authoritative.

Also established along the way: **the original WS-1 Bondi gate never produced a number.**
Job 2347716 was `CANCELLED AT ITS TIME LIMIT`.

**This was a self-inflicted repeat, and the project record had already warned about it.** The
WS-1 note states plainly, of increment 2: *"the plan's literal Bondi-rate gate is UNACHIEVABLE
with accretion off (no gas sink -> infall piles up and chokes; observed Mdot built to 0.41
Mdot_B then collapsed) -> that gate DEFERRED to after increment 5 (accretion); bondi.in +
compare_bondi.py ready to re-run then."* Attempt 1 re-ran `bondi.in` **unmodified** — i.e. the
exact configuration that was documented as unachievable — and reproduced the documented choking
almost exactly. Increment 5 has since landed, so the deferred gate is now runnable; the deck
just had to be switched to `accretion = true`, which is what "re-run then" meant.

Lesson, in the spirit of rule 1: reusing a validation deck means reading what its own project
record says about *why it was parked*, not only what it configures.

## Scope, stated honestly

`VALIDATION_PLAN` asks for "strict Shu-rate validation". The Shu (1977) inside-out collapse rate
Ṁ = 0.975 c_s³/G applies to a **singular** isothermal sphere, and no pgen in this tree builds one
(`collapse_be` makes a bounded-Emden sphere, `poisson_test` a uniform one). That IC is real work
and is **not** done here.

What is available, and is the same class of test, is **Bondi accretion** — steady spherical
accretion onto a point mass, with the exact rate

```
Mdot_B = 4 pi lambda_c (G M)^2 rho_inf / c_inf^3,   lambda_c(gamma=1.4) = 0.625
       = 49.087   for rho=1, c_inf=1, G=1 (four_pi_G=4pi), M=2.5  => r_B=2.5, r_sonic=0.5
```

The Shu geometry differs (self-gravitating envelope, expansion wave), but the accretion
**algorithm** under test is identical.

## Why attempt 1 measured nothing

Three legs (N = 64/128/256) all exited 0, running 342 / 1084 / 2616 cycles over 2.9 node-hours.
Exit code 0 is, once again, the worthless signal.

**Defect 1 — `sinks/accretion` defaults to `false` and `bondi.in` never set it**
(`sinks.cpp:81`). The sink exerted gravity but never removed gas. The result is unambiguous in
the N=128 dumps:

| shell | ρ(t=0) | ρ(t=6) | mean v_r(t=6) |
|---|---|---|---|
| r < 0.5 | 1.000 | **15.257** | **+0.060 (outward)** |
| 0.5–1.0 | 1.000 | 5.998 | **+0.106 (outward)** |
| 1.0–2.0 | 1.000 | 3.039 | +0.083 (outward) |
| 2.0–4.0 | 1.000 | 1.808 | −0.019 |
| 4.0–7.0 | 1.000 | 1.224 | −0.190 (inward) |

Gas piled up at the centre with no sink to remove it, central pressure rose, and the inner flow
**reversed** while the envelope still fell in. The measured Ṁ at r = 1 was therefore the
difference of an inflow and a pressure-driven rebound:

```
N=64   Mdot =  2.383  = 0.049 x Mdot_B
N=128  Mdot =  1.695  = 0.035 x Mdot_B
N=256  Mdot = -0.897  = -0.018 x Mdot_B     (sign change)
```

≈ zero, non-monotone, observed order p = −0.05. Nothing about the accretion algorithm was tested.
The sink **gravity** was, incidentally, confirmed working — ρ rose 1.0 → 15.3.

**Defect 2 — the accretion rate was unreadable.** Sink swarm values carry `Metadata::Restart` but
there is **no `Metadata::Output` flag for swarm values** in this Parthenon (verified: it does not
compile, and `metadata.hpp` has no such enumerator). The `.phdf` dumps therefore held only the
built-in `{id, x, y, z}` — sink `mass`, and hence dM/dt, the most direct accretion rate, was
absent from the science output. Parthenon selects swarm variables in the **output block**
(`outputs.cpp:316-336`), not via metadata:

```
<parthenon/output0>
swarm_variables = mass, vx, vy, vz, Lx, Ly, Lz, t_created
```

A pointer to this is now recorded at `src/sinks/sinks.cpp:43`; no functional code changed.

**Defect 3 — 29 % core utilisation.** One task, 8 OMP threads. Attempt 2 uses 16 MPI ranks.

## Attempt 2

Deck `bondi_wp17.in`, submit `submit_wp17b.sh`, `build_cpu` `63e954fc`, `tlim = 10` (up from 6:
t = 6 is only 2.4 sound-crossings of r_B, and with the sink now draining the centre a steady
state can actually form — the analysis prints Ṁ(t) per dump so the plateau is **checked**, not
assumed).

`racc_cells` is in **cells**, so the physical accretion radius shrinks as the grid refines — and
that is the whole ladder:

| N | dx | r_acc = 4dx | vs r_sonic = 0.5 |
|---|---|---|---|
| 64 | 0.2500 | 1.000 | **larger** — the sink swallows the sonic point; this leg *cannot* reproduce Bondi |
| 128 | 0.1250 | 0.500 | exactly at r_s, marginal |
| 256 | 0.0625 | 0.250 | inside r_s — the first leg that can host a genuine sonic transition |

Ṁ is read at **r = 1.5**: outside r_acc at every N (max 1.0) and inside r_B = 2.5. A rate that
does *not* converge with resolution would mean the accretion is set by the numerics rather than
the flow — the failure mode WP-17 exists to find.

**A further scheduling fault, for the record:** the first attempt-2 submission left the deck's
fixed 32³ meshblock, giving 8 blocks at N=64 under 16 ranks, and the script's comment asserted
this was "harmless — Parthenon just leaves ranks idle". It is not: the leg aborted immediately
with `### FATAL ERROR in CalculateLoadBalance` (exit 134, 0 cycles). Parthenon requires
nblocks ≥ nranks. `meshblock = N/4` now gives 64 blocks at every N, which also holds the
decomposition fixed across the ladder. The N=64 leg was relaunched separately (job 2453204).

---

# RESULT — N=128: Mdot = 1.01 +/- 0.02 x Mdot_B (job 2453375)

**The first measured sink accretion rate in this project.** The deferred WS-1 gate is answered:
the sink accretes at the analytic Bondi rate to within 1 %.

Deck `bondi_wp17.in`, `build_cpu` 63e954fc, 32 ranks, t = 0 -> 6, 13 dumps.

| t | Mdot | M_sink | Mdot_B(M(t)) | ratio |
|---|---:|---:|---:|---:|
| 1.51 | 0.03164 | 2.5179 | 0.049794 | 0.636 |
| 2.51 | 0.04339 | 2.5455 | 0.050892 | 0.852 |
| 3.50 | 0.05156 | 2.5824 | 0.052378 | 0.984 |
| 4.00 | 0.05425 | 2.6035 | 0.053236 | 1.019 |
| 4.50 | 0.05591 | 2.6257 | 0.054148 | 1.033 |
| 5.00 | 0.05660 | 2.6491 | 0.055117 | 1.027 |
| 5.50 | 0.05646 | 2.6735 | 0.056137 | 1.006 |
| 6.00 | 0.05579 | 2.6985 | 0.057192 | 0.975 |

**Read the RATE, not the ratio.** The measured Mdot is **flat to 1.44 %** over t = 4.5-6.0
(0.05591 / 0.05660 / 0.05646 / 0.05579) -- a genuine steady state. The *ratio* drifts down late
only because `Mdot_B ∝ M²` and the sink gained 4.5 % mass over that window, so the denominator
grows while the numerator holds. Quoting the endpoint (0.975) or the peak (1.033) alone would
both mislead; the honest number is the mean over the steady window:

> **Mdot / Mdot_B = 1.01 +/- 0.02** (t = 4.0-6.0)

The rho = 1e-3 design worked: sink mass grew **7.9 %** over the full run against the 12 % budgeted,
so the accretor stayed quasi-static and Bondi's fixed-point-mass assumption held.

## Convergence: a TWO-point ladder, not three

`r_acc = 4dx` shrinks with resolution, so the ladder was meant to be r_acc = 0.5 / 0.25 / 0.125 at
N = 128/256/512 against r_sonic = 0.5.

* **N = 256** (job 2453312) — **CONVERGED**. Measured against N=128 at matched times:

  | t | N=128 | N=256 | diff |
  |---|---:|---:|---:|
  | 0.50 | 0.286 | 0.289 | +1.05 % |
  | 1.01 | 0.487 | 0.489 | +0.41 % |
  | 1.51 | 0.636 | 0.638 | +0.31 % |
  | 2.01 | 0.755 | 0.759 | +0.53 % |
  | 2.51 | 0.852 | 0.857 | +0.59 % |
  | 3.00 | 0.928 | 0.934 | +0.65 % |
  | 3.50 | 0.984 | 0.991 | +0.71 % |

  **Mean offset +0.61 %, spread 0.73 %.**

  What matters is *what changed* between these two runs: `r_acc = 4dx` **halved**, 0.500 → 0.250,
  straddling the sonic radius (0.5). At N=128 the accretion sphere reaches exactly to the sonic
  point; at N=256 it sits well inside it. The rate moves **under 1 %** anyway.

  That is the answer to the question WP-17 exists to ask -- *is the rate set by the flow or by the
  numerics?* **By the flow.** A rate controlled by the accretion radius would have shifted
  substantially when r_acc halved. Together with N=128's absolute result
  (Mdot = 1.01 +/- 0.02 x Mdot_B) the test now has both correctness AND resolution-independence.

  The leg is CAPPED at t ≈ 4.4 by its 12 h wall limit, and that is final: `scontrol update
  TimeLimit=` returns **Access/permission denied** (users may lower a limit, not raise it), and
  `bondi_wp17.in` configures no restart output, so it cannot be resumed either. Do not re-attempt
  either route.
  This does **not** weaken the result. The convergence statement above comes from the MATCHED-TIME
  comparison over t = 0.5-3.5, where both legs have data; it never required N=256 to reach its own
  plateau. N=128 supplies the absolute rate (1.01 +/- 0.02), N=256 supplies resolution-independence
  (+0.61 %). A longer N=256 would add confirmation, not a new claim.
* **N = 512** (job 2453376) was **CANCELLED**. Measured rate: t = 0.134 after 3.66 h, i.e. **164 h
  needed against a 12 h limit**; it would have reached t = 0.44, far short of the t >= r_B/c = 2.5
  the flow needs merely to settle. 1.34e8 cells on CPU is not affordable for this test. A GPU port
  of the run would be the way to get a third rung; the sinks package is GPU-smoke-tested but its
  full-physics accretion path on GPU is not validated (WS-1 inc6), so that is not a free swap.

**N = 64 is excluded on physics, not cost:** at racc_cells = 4 it would have r_acc = 1.0, twice the
sonic radius, so the sink would swallow the sonic point and no Bondi solution exists outside it.

## Why this took six attempts

Every one of these produced a run that exited 0:

1. `sinks/accretion` defaults **false** and `bondi.in` never set it -- gravity without accretion;
   gas choked and the inner flow reversed (+0.06 outward while the envelope fell in at -0.19).
2. Sink `mass` absent from dumps -- there is no `Metadata::Output` for swarm values; it must be
   requested via `swarm_variables` in the output block (`outputs.cpp:316-336`).
3. rho = 1 gave a **49x timescale inversion** (sink e-fold 0.051 vs settling 2.5) -- runaway, sink
   mass 2.5 -> 910 by t = 3 with Mdot climbing monotonically to 31x Mdot_B.
4. Lowering rho to 1e-3 tripped the **Truelove/Jeans** `rho_sink` default, whose threshold scales
   as 1/dx² -- accretion went **inert at every resolution** (262x / 1047x / 4189x the ambient).
   Fixed with an explicit `rho_sink_code = 3e-5`.
5. My own `r_acc` fix starved it: the removal fraction is a quadratic ramp `((racc-r)/racc)²` that
   vanishes at the sphere edge, so racc_cells = 1 puts every cell at the edge (ramp 0.018, sink
   ~55x too weak). racc_cells must stay >= 4; the code default is well chosen.
6. My plateau watcher fired "RESULT IN HAND" on a still-**rising** curve because the last three
   points happened to fall inside a 5 % spread. A spread test is not a plateau test.

Plus the pre-existing fact that the original WS-1 Bondi gate (job 2347716) was cancelled at its
time limit and never produced a number at all.
