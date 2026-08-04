# D1 — per-rank GPU memory imbalance on deep AMR hierarchies

**Status: CHARACTERISED, two hypotheses FALSIFIED, root cause NOT established.** It is a real
defect that killed three runs on 2026-08-04. The practical mitigation (more ranks) works. This
document exists so the next person does not re-derive the two dead ends.

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
the B7 closure truncated at 120 of 500 cycles (job 2454557).

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

## Mitigation that works today

**Use more ranks.** Memory per rank falls roughly as 1/nranks while the *imbalance ratio* stays
put, so headroom grows: at 5 ranks the peak was 79.2 GiB of 80; at 7 ranks the same state projects
to ~57 GiB. WP-1's `cap=150` and the B7 closure rerun were both resubmitted at 7 ranks for exactly
this reason.

This is a workaround, not a fix: it wastes GPUs to buy headroom for one overloaded rank, and it
will fail again on a deeper hierarchy.
