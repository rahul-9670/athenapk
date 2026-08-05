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
