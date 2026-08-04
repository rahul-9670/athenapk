# WP-17 — Sink-particle accretion rate against an analytic answer

**Status: attempt 1 was a NULL TEST (job 2452595). Attempt 2 in flight (jobs 2453190 + 2453204).**

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
