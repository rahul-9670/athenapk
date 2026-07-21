# AthenaPK GPU Optimization Plan — first-core state (prod_v8)

**Status:** MEASURE + PLAN pass, read-only on source. Written 2026-07-21.
**Target:** `runs/prod_v8/` — 5×H100 on g003, `athenaPK_eos_v8_auditfixes` (md5 c8d5c2a8),
first-core formation, `wsec_step ≈ 17 s`, `zone-cycles/wsec ≈ 3.1e6`.

> **Reading guide.** Every optimization is tagged **(A)** = *result-preserving engineering
> win* (apply behind a bitwise / tight-tolerance gate) or **(B)** = *numerics/physics-affecting
> lever* (must be A/B-validated against a frozen reference and reported with BOTH throughput
> AND simulated-time — never quoted as free). Numbers are labeled **measured** (read this
> session), **derived** (computed from measured inputs), or **inferred** (from code structure /
> historical probes; confidence stated). The kernel/API/occupancy micro-breakdown is **PENDING**
> the Nsight job (see §7) — the cluster had **zero free GPUs** during this pass so the profile
> could not run yet. This document is structured so the PENDING numbers drop into §3/§7 without
> reworking the lever ranking, which rests on the step *anatomy* (measured counts), not on
> per-kernel fractions.

---

## 1. What was measured this session (evidence base)

All from the LIVE run (`runs/prod_v8/run.log`, cycles ~1160–1180) and from a **copy** of restart
`parthenon.out2.00004.rhdf` (cycle 1000) placed in `runs/gpu_profile_v8/restart_c1000.rhdf`.
The live run, its binary, inputs, restarts and logs were **not** touched.

| Quantity | Value | Source |
|---|---|---|
| `wsec_step` (wall/step) | **16.1–18.8 s, median ~17.3 s** | run.log, measured |
| `zone-cycles/wsec` | **2.5e6–3.4e6, median 3.15e6** | run.log, measured |
| `wsec_AMR` (regrid/step) | **1.42 s = 8.2% of step** | run.log, measured |
| `dt` | **2.49e-7 – 3.84e-7, median 3.02e-7** code | run.log, measured |
| STS "Taking N steps" | **17 / 19 / 21** (mode 19) | run.log, measured |
| STS invocations **per hydro step** | **2** (80 "Taking" lines / 40 cycles) | run.log, measured |
| ⇒ diffusion sub-stages per step | **≈ 2 × 19 = 38** | derived |
| MeshBlocks | **1807** (32³ each) | restart Info, measured |
| MaxLevel | **10** (of 19 allowed; numlevel=20) | restart Info, measured |
| blocks/level | L0:485 … L8:132, L9:94, L10:16 | restart Levels, measured |
| total zones | **5.92e7** (1807 × 32768) | derived |
| blocks/rank | **≈ 361** (1807 / 5) | derived |
| finest `dx_min` | **6.10e-5** code | derived (16/256/2¹⁰) |
| radiation `chat` | **c_code/creduc = 1.58e6/1000 ≈ 1579** code | `src/radiation/radiation.cpp:100`, measured |
| radiation `nsub` per step | **≈ 20** at this state | derived: `ceil(dt·chat/(cfl·dx_min))` |
| gravity Poisson solves / step | **2** (stage-consistent, audit #2) | `hydro_driver.cpp` §gravity, measured |

**GPU-utilization ground truth (user, live `nvidia-smi dmon`):** 5 GPUs oscillate ~20% ↔ 84–99%
sm-util, in phase — peak compute is healthy; the time-average is diluted by **overhead valleys**
(sync / comm / AMR / gravity-solve). This plan is about the valleys and the *count* of whole-mesh
sweeps, not peak kernel efficiency.

---

## 2. Per-step anatomy (why 17 s)

One `vl2` hydro step at the cycle-1000 state issues, over the whole mesh (5.92e7 zones):

1. **Strang STS half #1** (`AddSTSTasks`, `hydro_driver.cpp` stage==1, 0.5·dt): **~19 RKL2
   sub-stages**. Each sub-stage = `CalcDiffFluxes` (parabolic Ohm+AD, Hall floor) → RKL2 update →
   **full-container boundary exchange** → **full-mesh `FillDerived`**.
2. **2 VL2 hydro stages**: HLLD MHD flux + flux-correction reflux + unsplit sources (incl. unsplit
   Hall EMF) + full boundary exchange.
3. **Self-gravity**: **2 BiCGSTAB+MG Poisson solves** (audit #2 stage-consistent), one per VL2
   stage. Each solve = N_iter × (MG V-cycle of many *tiny* per-level kernels + **2 MPI all-reduces**
   for the Krylov inner products) + a phi halo exchange.
4. **Radiation M1**: **~20 sub-cycles** (final stage), each = `CalculateRadFluxes` → reflux
   (audit A1) → `ApplyRadUpdate` → **rad boundary exchange**; then one `MatterCoupling` implicit
   solve over the full dt.
5. **Chemistry** (gow17-reduced semi-implicit) + **EOS** tabulated-hydrogen `FillDerived`.
6. **Strang STS half #2** (0.5·dt): another **~19 RKL2 sub-stages** as in (1).
7. **AMR regrid**: measured **1.42 s (8.2%)**.

The step is dominated by the **count of whole-mesh sweeps**: **≈38 diffusion sub-stages + ~20
radiation sub-cycles + 2 hydro stages + 2 Poisson solves ≈ 60+ whole-mesh passes**, each dragging
its own boundary exchange. The 20%-util valleys the user sees are the exchanges + the two
all-reduce-bound Poisson solves interleaved between compute bursts. **AMR is NOT the bottleneck**
(measured 8.2%) — consistent with the previously-FALSIFIED `amr_check_interval` speedup.

---

## 3. Kernel / API / load-balance micro-breakdown — PENDING (job 2377882)

The `nsys` profile (§7) will fill this table. Structure reserved so the ranking in §5 does not
change when it lands:

| Bucket | Expected top contributors (inferred) | Measured GPU-time % | Confirm/refute |
|---|---|---|---|
| Diffusion STS | `PrecomputeNonidealEta`, `Ambipolar/Hall/Ohmic X{1,2,3} fluxes`, RKL2 update, `FillDerived` | PENDING | is ~38× the driver? |
| Radiation | `CalculateRadFluxes`, `ApplyRadUpdate`, `MatterCoupling` | PENDING | nsub = instances/step |
| Gravity | MG smoother/restrict/prolong (many tiny), BiCGSTAB axpy/dot, `MPI_Allreduce` | PENDING | all-reduce wall? |
| Hydro | HLLD flux, reconstruction, flux-correction | PENDING | |
| Boundary exchange | pack/unpack + `MPI` + device syncs | PENDING | valley cost |
| CUDA API | `cudaDeviceSynchronize`, `cudaMalloc/Free` (post-setup churn?), stream/event sync | PENDING | malloc after setup? |
| Load balance | per-rank total GPU time (r0..r4 spread) | PENDING | SFC imbalance |

**Historical measured anchors** (prior probes, from run notes — treat as inferred for v8):
`PrecomputeNonidealEta` was **65.7% of GPU kernel time** before `rkl2_freeze_eta` (which then gave
**2.18× wsec_step**); `eta_ohm_cap_code=0.1` gave **6.1× sim-time/wall**. These two big wins are
**already applied** in prod_v8 — the low-hanging eta fruit is spent, which is why the remaining
levers target the *sweep count* and the *per-sweep overhead*.

---

## 4. HONEST HEADLINE — what actually moves time-to-second-core

**Time-to-second-core = (remaining sim-time / dt) × wsec_step.** Two independent factors:

- **`wsec_step`** (wall per step) — lowered by **(A) engineering** levers. These are pure win:
  they cut wall time proportionally and **preserve results**. Ceiling for this pass ≈ **1.3–1.8×**
  (see §5), concentrated in *not repeating whole-container work 38× per step*.
- **`dt`** (⇒ number of steps) — raised only by **(B) physics** levers. `dt` here is **NOT**
  hydro-CFL-limited; it is pinned by the **unsplit dispersive Hall/whistler term** (dt ∝ dx²) plus
  the parabolic STS budget. The big historical dt wins (`eta_ohm_cap`, freeze_eta) are already in.
  Remaining dt levers (Hall cap, larger `creduc`, Poisson tolerance) are physics and were mostly
  **priced-and-rejected** already (Hall cap: binding |η_H|~7e-3 ≪ any safe cap → 0 gain).

**The single biggest *safe* lever on wall-clock-to-second-core is cutting the per-STS-stage
overhead** (Lever 1): ~38 full-container boundary exchanges + ~38 full-mesh `FillDerived` per step,
when the parabolic stencil touches only `cons`(B) and the frozen-eta stages don't need a fresh EOS
inversion. This is (A), result-preserving, and the code's own `TODO(pgrete)` flags exactly it.

**Raising GPU% is NOT the goal.** Filling the 20% valleys by, e.g., overlapping the gravity
all-reduce would raise average sm-util but only helps time-to-second-core to the extent it shortens
`wsec_step`. Report every lever in `wsec_step` / sim-time-per-wall, not in sm-%.

**Skeptic's note:** the `amr_check_interval` lever was measured and FALSIFIED (bit-identical, zero
wall saving — AMR cost is coupled to the step, not gateable slack). Treat every payoff below as a
hypothesis until the §7 profile (and, for (B), an A/B against a frozen reference) confirms it.

---

## 5. Ranked levers

### Lever 1 — Slim the per-STS-stage boundary exchange + FillDerived  **(A)** ★ top priority
- **Evidence:** ~38 diffusion sub-stages/step (measured), each does
  `AddBoundaryExchangeTasks(..., base, ...)` on the **full base container** (cons, prim, grav.phi,
  rad.Er/Fr, nonideal_eta, 5 scalars) and a full-mesh `Update::FillDerived`. The parabolic RKL2
  stencil only needs the **magnetic field** ghosts; with `rkl2_freeze_eta=true` the eta cache is
  frozen so the EOS-inversion part of `FillDerived` is redundant in the frozen stages.
- **Code:** `src/hydro/hydro_driver.cpp` — `AddSTSTasks`: L285–289 (first sub-stage exchange +
  FillDerived) and L351–357 (loop sub-stage exchange + FillDerived). Both carry
  `// TODO(pgrete) optimize (in parthenon) to only send subset of updated vars`. Register
  allocation L232–234 (`MY0`, `Yjm2`) also over-allocates prim+cons (`TODO(pgrete) ... allocate
  only required vars`).
- **Payoff (inferred, medium confidence):** exchanging only the evolved MHD `cons` (or just B) and
  skipping the redundant prim-inversion in frozen stages removes a large fraction of the ~38
  per-stage exchange+FillDerived passes. If exchange+FillDerived is ~30–45% of `wsec_step`
  (PENDING §3 to confirm), a 2–3× reduction of that slice ⇒ **wsec_step 17 s → ~11–13 s
  (1.3–1.5×), same sim-time/step ⇒ 1.3–1.5× faster to second core.** Result-preserving.
- **Validation gate:** bitwise-identical gas state vs frozen reference for ≥100 steps
  (`cons` + `grav.phi` + `rad.Er` checksums). The narrowed exchange must include EXACTLY the ghost
  vars the parabolic flux + next reconstruction read (B, and the x_e scalar read one cell into the
  ghost by the first RKL2 non-ideal flux — see the audit-A3 comment at `hydro_driver.cpp` §fill).

### Lever 2 — Narrow the STS registers (`MY0`, `Yjm2`) to the evolved vars  **(A)**
- **Evidence / code:** `hydro_driver.cpp:232–234` adds `MY0`/`Yjm2` as full copies of `base`
  (prim+cons+rad+grav+scalars) though RKL2 only recurses on `cons`. The two `DeepCopy`s at
  L215–221 and the init copies likewise move full prim+cons.
- **Payoff (inferred, low-medium):** cuts device memory footprint (helps the L≥13 headroom problem
  the run notes flag) and removes prim/rad/grav bytes from every RKL2 `DeepCopy`/update. Wall win
  is secondary to Lever 1 but it is free memory and bandwidth. **Quantify from §7 `gpumemtimesum`
  + the DeepCopy kernel time.**
- **Gate:** bitwise gas state; check peak `nvidia-smi` mem/card drops.

### Lever 3 — Gravity Poisson: cut solves and/or all-reduce latency  **(A) for engineering parts / (B) for tolerance**
- **Evidence:** 2 BiCGSTAB+MG solves/step (audit #2), each with ~2 `MPI_Allreduce`/iteration
  (`src/self_gravity/self_gravity.cpp:216–240`, note at L218–227). BiCGSTAB inner products are the
  classic multi-GPU latency valley; the many small MG-level kernels are launch-bound.
- **(A) engineering:** (i) **warm-start** the 2nd (corrector) solve from the predictor's φ (the run
  already keeps φ across cycles post-v6 named-pack fix) to cut N_iter; (ii) **CUDA-graph / fused
  launch** of the fixed MG V-cycle to kill per-tiny-kernel launch latency; (iii) evaluate the
  `solver=MG` pure-multigrid path (no global inner products — self_gravity.cpp:218–227 documents it)
  as an all-reduce-free preconditioner-solver — **but** it needs an adequate AMR smoother, so this
  is a numerics change ⇒ treat as **(B)** and A/B it.
- **(B) tolerance:** loosening the Poisson residual tolerance changes the gravity field ⇒ collapse
  timing/energetics. Do NOT touch without an A/B vs frozen reference.
- **Payoff:** PENDING §7 — need the measured all-reduce wall and N_iter/solve. If the 2 solves are
  ~15–25% of the step (plausible given the valleys), warm-start + graph capture could recover a
  meaningful slice **result-preservingly** (i, ii) — the tolerance/MG-solver parts are (B).

### Lever 4 — Radiation sub-cycling (~20/step)  **(B)** — physics, not free
- **Evidence:** `nsub = ceil(dt·chat/(cfl·dx_min))` (`radiation_moments.cpp:460`), `chat=c/creduc`,
  `creduc=1000`. **Measured/derived nsub ≈ 20** at cycle 1000; each sub-cycle is a full-mesh flux +
  reflux + exchange. This is ~20 of the ~60 whole-mesh passes/step.
- **Levers (all (B), must A/B + report sim-time):**
  (a) **raise `creduc`** (e.g. 1000→2000–3000): `chat` still ≫ gas/collapse speeds (v_gas ~ few
  code units ≪ 1579), so RSLA likely stays valid; nsub ∝ 1/chat ⇒ **~2–3× fewer rad sub-cycles**.
  Changes the radiation equilibration timescale ⇒ validate front-speed + matter-coupling against a
  frozen reference (`chat=c_code/creduc` trap: creduc is a *reduction factor*, not a speed).
  (b) **level/local sub-cycling** — only sub-cycle rad on blocks whose `dx` demands it (coarse
  blocks take fewer sub-steps), instead of the global `dx_min`. Larger change; also (B).
  (c) asymptotic-preserving / implicit transport — research-scale (B).
- **Payoff:** (a) is the cheapest big rad win; if radiation is ~25–35% of the step (PENDING §3),
  2× fewer sub-cycles ⇒ ~12–17% off `wsec_step`, **but** it perturbs radiation dynamics — quote
  throughput AND the sim-time/energetics A/B, never as free.
- **Scaling caveat (inferred):** if `dt` is whistler-limited (dt ∝ dx²) while `dt_rad ∝ dx`, then
  nsub ∝ dx and **shrinks** as the core deepens — i.e. radiation may matter *less* at L≫10. Confirm
  the dt-limiter in §7 before investing here; measure, don't assume.

### Lever 5 — Post-setup device malloc/free + redundant syncs  **(A)**
- **Evidence:** t1's OOM was allocator growth under AMR churn (run notes) — a symptom of
  per-cycle device alloc/free. `gpumem.log` is already collected in prod. PENDING §7
  `cudaapisum` to quantify `cudaMalloc/cudaFree`/`cudaDeviceSynchronize` time after setup.
- **Payoff:** if there is measurable post-setup malloc/free or redundant device/stream syncs
  between the ~60 sweeps, pooling/removing them is pure (A) valley-fill. Quantify first.

### Lever 6 — Kokkos team/vector tuning + larger blocks  **(A)**, low priority
- **Evidence:** blocks are 32³. PENDING `ncu` occupancy/registers on the top kernels (§7) will say
  whether the flux/eta kernels are occupancy- or register-limited and whether 64³ blocks (fewer,
  larger launches → less launch overhead, better AMR-boundary/compute ratio) help. `njeans=8`
  refinement interacts with block size — changing nx1/2/3 alters the AMR block *count* and the
  decomposition, so a 64³ test is A but must be checked bitwise + for memory.

---

## 6. Lever summary

| # | Lever | Class | Touches physics? | Expected `wsec_step` win | Confidence |
|---|---|---|---|---|---|
| 1 | Slim per-STS-stage exchange + FillDerived | **A** | No | **1.3–1.5×** | medium (needs §3) |
| 2 | Narrow STS registers MY0/Yjm2 | **A** | No | mem + minor wall | low-med |
| 3 | Gravity: warm-start + graph capture | **A** | No | slice of ~15–25% | PENDING §7 |
| 3′| Gravity: MG solver / looser tol | **B** | **Yes** | (as above) | A/B required |
| 4 | Radiation: raise `creduc` / level-subcycle | **B** | **Yes** | ~12–17% (if nsub big) | A/B required |
| 5 | Kill post-setup malloc/free + syncs | **A** | No | valley-fill | PENDING §7 |
| 6 | Kokkos team/vector + 64³ blocks | **A** | No | occupancy-dependent | PENDING ncu |

**Net safe (A) ceiling this pass ≈ 1.3–1.8× `wsec_step`** (Levers 1+2+3-eng+5+6), i.e. a
proportional 1.3–1.8× cut in wall-clock-to-second-core with results preserved. (B) levers 3′/4 are
additional but must each be A/B-validated and quoted with sim-time, never stacked as "free."

---

## 7. PENDING measurement — Nsight job 2377882 (queued)

**Why pending:** during this pass **all GPU nodes were 8/8 allocated** (g001/g003/g004 ALLOCATED,
g002 gres/gpu=8 alloc) — **zero free GPUs cluster-wide**. The profile needs 5×H100 (the 1807-block
first-core mesh will not fit on fewer cards → OOM), so it can only backfill onto g003 after the
prod slot ends. `squeue --start` estimate: **~03:53** (~8 h out). The job is queued (does not touch
prod files; profiles a *copy* of the restart).

**Artifacts staged in `runs/gpu_profile_v8/`:**
- `restart_c1000.rhdf` — copy of prod restart 00004 (cycle 1000, representative first-core, 38 STS
  sub-stages/step). Originals untouched.
- `libkp_nvtx.so` (+ `kp_nvtx.cpp`) — Kokkos-Tools→NVTX connector (Kokkos 4.7.2, `KOKKOS_TOOLS_LIBS`)
  so nsys shows physics **labels** (`PrecomputeNonidealEta`, `Ambipolar/Hall X{1,2,3} fluxes`,
  `CalculateRadFluxes`, MG kernels) instead of mangled functor symbols.
- `submit_nsys.sh` + `wrap_nsys.sh` — 5-rank nsys job, restart + **3 steps** (nlim=1003; 1001
  warmup + 2 clean), outputs disabled, **physics CLI mirrors prod submit.sh** (rkl2 / freeze_eta /
  eta_ohm_cap / ion_zeta / coalesced comms). Traces cuda,nvtx,mpi; per-rank `.nsys-rep`.
- `submit_ncu.sh` (to author) + `wrap_ncu.sh` — ncu on **rank 0 only**, `--set basic --nvtx`,
  `--launch-skip 300 --launch-count 40` (skip warmup). After nsys names the top 5–10 kernels,
  retarget with `--nvtx-include "<label>/"` and `-c 2` each.
- `parse_nsys.sh` — runs `nsys stats` (cuda_gpu_kern_sum, nvtx_sum, cuda_api_sum,
  cuda_gpu_mem_time_sum) to emit the §3 table. **radiation nsub = `CalculateRadFluxes` instances /
  steps profiled; load balance = per-rank total kernel GPU time (r0..r4).**

**When the profile lands, fill §3 and reconfirm/deny:** (i) is diffusion STS the top bucket at ~38×
the driver? (ii) exact nsub; (iii) gravity all-reduce wall + N_iter/solve; (iv) any post-setup
`cudaMalloc/Free`; (v) per-rank load spread; (vi) ncu occupancy/registers on the flux/eta kernels.
Then re-rank §5 if the measured fractions contradict the inferred ones — and remember the
`amr_check_interval` lesson: **measure the win, don't assume it.**

**Scheduling note / trade-off:** the queued profile (`nsys_v8`, 2377882) and the prod chain's
successor (`prod_v8`, 2376457, dependency-held) both want g003 when 2376456 ends. If the profile
grabs the node first it delays the next prod slot by at most its ~10–15 min actual runtime (1 h
wall cap, releases early). To fully protect the prod slot, `scancel 2377882` and resubmit when a
GPU node is genuinely idle.
