# D1 — per-rank GPU memory imbalance on deep AMR hierarchies

**Status: CHARACTERISED, THREE hypotheses FALSIFIED, root cause NOT established, and there is NO
known mitigation.** A real defect that killed four runs on 2026-08-04. This document exists so the
next person does not re-derive the dead ends — including the "just use more ranks" one, which I
asserted here earlier and then disproved (§ Hypothesis 3).

## The observation

On the level-13 restart (`prod_t4_full/parthenon.out2.00200.rhdf`, 2486 blocks) and on the
deep_amr restart (792 blocks), per-GPU memory is grossly uneven and one rank OOMs while others
have tens of GiB free:

```
physical device 0 : 49.3 GiB     (rank 0, blocks 0:197)
physical device 1 : 56.6 GiB     (rank 1, blocks 198:395)
physical device 2 : 79.2 GiB     (rank 2, blocks 396:593)   <- OOM, "rank 2 exited on signal 6"
physical device 5 :  0.0 GiB     (rank 3, blocks 594:791)   <- sampled after the abort
```

**A 1.6× spread on identical block counts.** Failure mode is always the same:
`Kokkos ERROR: Cuda memory space failed to allocate <a few MiB>` — the card is saturated, and the
allocation that happens to fail is tiny and arbitrary (`bnd_flux::grav.phi`, `prim`,
`bnd_flux::rad.Fr2_g1`).

**Casualties:** WP-1 `cap=150` twice (jobs 2454256 at 5 ranks, and the leg inside 2454224), and
BOTH B7 closure runs — job 2454557 at 4 ranks and job 2454734 at 5 ranks — each truncated at
exactly 120 of 500 cycles.

## Hypothesis 1 — block-distribution imbalance. FALSIFIED.

Parthenon prints the assignment, and it is **exactly equal**:

```
Blocks assigned to rank 0: 0:197        Blocks assigned to rank 1: 198:395
Blocks assigned to rank 2: 396:593      Blocks assigned to rank 3: 594:791
```

**198 blocks each.** Equal counts cannot produce a 1.6× memory spread. (I initially reported this
as the cause — it is not.)

Relevant context from the source: `mesh.cpp:1043-1049` shows the **timing-based balancer is
disabled outright** (`PARTHENON_FAIL("Timing based load balancing is currently unavailable.")`),
so only the count/cost-equal `default` is available; `parthenon/loadbalancing/tolerance` defaults
to 0.5 and `interval` to 10.

## Hypothesis 2 — the pinning wrapper aliases devices. FALSIFIED.

Two ranks sharing one GPU would explain both the saturated device and the idle one. Tested
directly (`runs/deep_amr/probe_gpumap.sh`, job 2454732): each rank prints SLURM's allocated list
and what `wrap_mod.sh` selects from it.

```
rank 0  slurm_list=[0,1,2,3]  selected=0
rank 1  slurm_list=[0,1,2,3]  selected=1
rank 2  slurm_list=[0,1,2,3]  selected=2
rank 3  slurm_list=[0,1,2,3]  selected=3
```

**All distinct — the wrapper is correct.**

> **A trap inside the diagnostic itself.** The first version of this probe queried
> `nvidia-smi --query-gpu=uuid | head -1` after setting `CUDA_VISIBLE_DEVICES`, and every rank
> returned the *same* UUID at index 0 — which looks exactly like device aliasing. It is not:
> **`nvidia-smi` ignores `CUDA_VISIBLE_DEVICES`** and enumerates every device it can see
> (confirmed: `nvidia-smi_sees=5 devices` while the job was given 4, with the env pinned to one).
> `head -1` therefore always returned physical GPU 0. **Do not use `nvidia-smi` to ask what a rank's
> CUDA context sees** — it answers a different question. A correct test needs the CUDA runtime.

## Also note: RRZ device indices are PHYSICAL, SLURM's are the job's view

The job was given `CUDA_VISIBLE_DEVICES=0,1,2,3` and the RRZ report lists devices `0,1,2,5`. Job
device 3 is physical 5. Reconciling the two is necessary before attributing memory to a rank; the
`0.0 GiB` on physical 5 is rank 3 sampled *after* rank 2's abort killed the job, not an idle rank.

## What remains — the leading hypothesis, untested

Equal block count does **not** imply equal memory. Parthenon allocates communication buffers per
*neighbour relationship*, and on an AMR hierarchy a rank spanning refinement boundaries also
carries prolongation/restriction buffers. Block ranges are assigned as contiguous **Z-order**
intervals, and those intervals can differ substantially in how many refinement boundaries and
neighbours they touch. So the suspect is **buffer memory scaling with subdomain surface and
refinement-boundary count, not with block count**.

Testing it needs per-rank buffer accounting, which Parthenon does not currently report. The
cheapest next step is a `Kokkos` allocation dump or a per-rank `View` inventory at startup,
attributing bytes by label — the failing labels are already known to be `bnd_flux::*`, which is
itself consistent with this hypothesis.

## Hypothesis 3 — "more ranks buys headroom". FALSIFIED 2026-08-04.

I originally documented this as the working mitigation, reasoning that per-rank memory falls as
1/nranks. **It does not help at all.** Two runs of the identical B7 closure:

| run | ranks | outcome |
|---|---|---|
| `b7_closure` (job 2454557) | 4 | OOM at **cycle 120**, `bnd_flux::rad.Fr2_g1` (1.068 MiB) |
| `b7_closure2` (job 2454734) | 5 | OOM at **cycle 120**, `bnd_flux::cons.coarse` (2.93 MiB) |

**The same cycle, to the cycle, with 25 % more GPUs.** A per-rank-load problem would have moved
the failure point. This one did not budge.

## What that implies — and it now points at a mechanism

Both runs follow an identical mesh history: **792 → 883 blocks (row 16) → 939 blocks (row 61)**,
then die at row 120. The failing allocations are always `bnd_flux::*`, and on the second run
specifically **`cons.coarse`** — the COARSE buffer used for AMR prolongation/restriction across
refinement boundaries.

That is consistent with the surviving hypothesis and sharpens it: the exhausted memory is
**AMR coarse/prolongation buffer space, sized by the mesh's refinement structure rather than by
each rank's block count**. Rank count changes how blocks are distributed; it does not change how
many coarse-fine boundaries the hierarchy contains. Hence adding ranks moves nothing.

It also explains the earlier per-device spread (49.3 / 56.6 / **79.2** GiB on equal block counts):
ranks holding Z-order intervals that straddle more refinement boundaries carry more coarse-buffer
memory, independent of how many blocks they own.

**There is currently NO known mitigation.** More ranks does not work. The practical consequence is
that deep-AMR runs on this hierarchy terminate at ~cycle 120 from the deep restart, whatever the
decomposition. Both B7 closure runs got their science anyway — the quadrature had converged well
before the crash — but a run needing to go further will simply stop there.

**Next step for whoever picks this up:** instrument the `bnd_flux::*coarse` allocations directly
(a Kokkos allocation dump by label at each regrid) and check whether their total scales with the
coarse-fine boundary count rather than with blocks/rank. If so, the fix is in Parthenon's buffer
sizing, not in job geometry.

## The instrument (2026-08-05) — implemented, not yet run

`<hydro> d1_meminfo` (default **false** => no output, no cost, bit-identical) makes every rank
print, whenever the block count changes (i.e. after each regrid):

```
[D1] cycle=N rank=R nblocks=... coarse_fine_nbrs=... same_level_nbrs=... dev_free_GiB=... dev_total_GiB=...
```

`coarse_fine_nbrs` counts this rank's neighbours sitting at a **different** refinement level;
`dev_free/total` come from `cudaMemGetInfo` under `KOKKOS_ENABLE_CUDA`. Implementation:
`src/hydro/hydro.cpp`, `PreStepMeshUserWorkInLoop`.

**The prediction to test.** If the surviving hypothesis is right, `dev_free_GiB` tracks
`coarse_fine_nbrs` and **not** `nblocks`. That is precisely what would explain why 4 and 5 ranks
die at the same cycle: rank count changes how blocks are distributed, not how many coarse-fine
boundaries the hierarchy contains. If instead `dev_free_GiB` tracks `nblocks`, the hypothesis is
dead and the search restarts.

To collect: add `-hydro/d1_meminfo=true` to the CLI of a deep-AMR restart and grep `[D1]` out of
the job log. It needs a GPU run; none has been launched.

---

## FIRST INSTRUMENT OUTPUT (2026-08-06, job 2468612, 4-rank leg) — the prediction fails, and a better hypothesis appears

The instrument ran for the first time. 8 `[D1]` samples at cycles 250 and 269. Two results, one
methodological correction, and a new leading hypothesis.

### 1. Per-rank coarse-fine sizing is DEAD

Parthenon load-balances to **equal block counts**, which turns each cycle into a clean
fixed-`nblocks` contrast — the opposite of the collinearity the original plan assumed. At cycle
250 every rank holds exactly 198 blocks while `coarse_fine_nbrs` spans 1381–1682:

| cycle | nblocks | cf_nbrs spread | consumed spread | amplification | r |
|---|---|---|---|---|---|
| 250 | 198 on all 4 ranks | +21.8 % | **+0.9 %** | **0.040** | −0.284 |
| 269 | 220–221 (0.45 %) | +19.8 % | **+0.8 %** | **0.040** | −0.643 |

A ~20 % swing in a rank's own coarse-fine boundary count moves its memory by <1 %, *non-monotonically*
(the highest-`cf` rank consumes the least), with the same amplification 0.040 at both cycles.

**Scope — do not overstate this.** It kills the **per-rank** form only. A *global* form (buffers
sized by the whole mesh's coarse-fine count, therefore identical on every rank) is invisible to a
within-cycle contrast and remains open. Only the 4-vs-5-rank leg comparison separates global-cf
from blocks.

### 2. The unmodelled term: memory grows 7.4× faster than blocks

| cycle 250 → 269 | change |
|---|---|
| blocks/rank | +11.5 % |
| coarse_fine_nbrs | +12.7 % |
| **consumed memory** | **+84.8 %** |
| amplification vs blocks | **7.4×** |

Consumed memory nearly **doubles** (18.4 → 34.0 GiB) over 19 cycles while neither per-rank
predictor moves more than 13 %. **Neither hypothesis was framed to explain this**, and it is the
more likely story: something scaling with mesh **depth over time** — level count,
prolongation/restriction buffer growth, or a leak — rather than with how existing work is divided.
That also explains the observation that killed hypothesis 3: if memory tracks depth rather than
division, adding ranks cannot help, and 4 and 5 ranks OOM **at the same cycle** exactly as seen.

Next test: `d1_meminfo` currently prints only when a rank's block count *changes*. Make it print
unconditionally every N cycles to resolve the growth curve between regrids, and add a Kokkos
allocation dump by label so the growing allocation is named rather than inferred.

### 3. Methodological correction to `d1_analyse.py`

The original cross-leg test compared **leg means** of consumed memory. That is invalid here:
consumed memory is dominated by *when* a sample was taken (18.4 GiB at cycle 250, 34.0 at 269),
and samples are emitted only when a rank's block count changes — so the sampling cycles are
themselves a function of rank count, which is the very thing the test varies. A leg mean compares
sampling history, not the hypothesis. The script now compares at **matched cycles** and refuses a
verdict if the legs share none. It also gained the fixed-`nblocks` within-cycle discriminator
above, which needs only one leg.

---

## COMPLETED LEG (2026-08-06, job 2468612) — **THE RUN DID NOT OOM**, and this supersedes the partial-data section above

The 4-rank leg finished: `exit=0`, `last_cycle=380` (stopped by `nlim`, **not** by an allocation
failure), `oom_lines=0`, peak **61.0 of 79.7 GiB**. It ran straight through cycle 369, the point at
which the b7_closure runs died. Three corrections follow, two of them to what *this document* said
a few hours earlier on partial data.

### 1. The founding observation no longer reproduces

| | original OOM (2026-08-04) | this leg (2026-08-06) |
|---|---|---|
| blocks/rank | 198 (792 total) | **378–379 (1513 total)** |
| per-rank consumed | 49.3 / 56.6 / **79.2** GiB → rank 2 died | **43–50 GiB**, no death |
| spread across ranks | 61 % | **16 %** |

Nearly **double the blocks, less memory per rank, and four times better balanced.** The mesh is
not the variable that changed.

### 2. Every OOM casualty had DIAGNOSTIC PACKAGES on; this leg had none

| script | diagnostics |
|---|---|
| `submit_wp1_bind150.sh` | `cap_diag mag_diag` |
| `submit_wp1_bind150_r7.sh` | `cap_diag mag_diag` |
| `submit_wp1_binding.sh` | `cap_diag mag_diag` |
| `submit_b7_closure.sh` | `cap_diag mag_diag cons_diag angmom_diag solver_diag` |
| `submit_b7_closure2_r5.sh` | the same five |
| **`submit_d1_leg.sh` (survived)** | **none** |

Four for four. If those packages allocate per-block fields, they are a large and entirely
avoidable addition to the footprint — and D1's premise (AMR coarse/prolongation buffer sizing)
is the wrong tree. **This is correlational**, with two confounds: the surviving leg also used a
newer binary (`84a6d248`, carrying the 2026-08-05/06 audit batch) and 4 ranks rather than 5.
`submit_d1_diagab.sh` (job 2468980) removes both — same binary, same restart, same rank count,
same `nlim`, diagnostics the *only* difference.

### 3. Two claims made earlier today from partial data are WITHDRAWN

- **The "7.4× unmodelled term" is falsified.** With all 16 samples the sequence is
  **7.4× → 1.3× → 0.5×**. The later intervals are at or *below* linear in blocks. 7.4× was one
  early interval, not a trend. (It was labelled "leading hypothesis to test next, not a result" —
  the hedge was right, and the test has now killed it.)
- **"Per-rank coarse-fine sizing is dead" was drawn from shallow cycles only.** At 198–235
  blocks/rank the amplification is 0.040–0.086 with inconsistent sign (r = −0.28, −0.64, +0.55),
  but at **378 blocks/rank it is 0.531 with r = +0.990**. The dependence **emerges with depth**.
  The hypothesis is not dead — it is depth-dependent, and the shallow cycles cannot refute it.
  `d1_analyse.py` now reports this verdict per depth and refuses to pool. Caveat stated plainly:
  the deep regime is **one cycle**, so it is a strong hint, not an established scaling.

**Lesson for the next person:** every number in the superseded section came from a run that was
still going. Amplifications and correlations computed mid-run here are provisional — the sign of
the correlation flipped and the amplification fell by 15× once the leg finished.
