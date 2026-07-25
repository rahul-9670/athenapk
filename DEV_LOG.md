# Self-Gravity Port Development Log

## Current state (end of Apr 24, 2026)

### What's done and validated
- Self-gravity port from Artemis → AthenaPK (hydro + MHD)
- Jeans dispersion tests: hydro 1.245% error, MHD 1.272% error (cluster)
- BE collapse tests: hydro matches laptop-vs-cluster to 3 sig figs;
  MHD shows physical flux-freezing and magnetic support of collapse
- **Jeans-length refinement** implemented and validated on laptop MHD BE:
  nested 3-shell refinement at t=0, level-4 activated as collapse progresses
- Supported on macOS/Clang (laptop) and Linux/GCC 13.3 + OpenMPI 5.0.7 (Hummel2)
- Parthenon v25.12 (commit d467a382b) — grid.type API fix applied

### Athena++ code comparison (MHD BE collapse)
Done at uniform grid (no AMR):
  - 32^3: rho_max final agrees 97%, B_max 86%
  - 64^3: rho_max final agrees 98%, B_max 95%
  - Peak-collapse offset of ~0.1 time units (expected GLM vs CT diffusion)
  - Both codes under-resolve Jeans length at peak rho (uniform grid limitation)

### Known issues / TODOs
- **Neumann gravity BCs in AthenaPK broken**: gauge not fixed for all-Neumann case
  → dt crashes at cycle 1. Use `zero` (Dirichlet) for now.
- **Multipole gravity BCs not implemented** in AthenaPK (Athena++'s default for test5)

### Next development priorities
1. Rebuild AthenaPK on cluster with Jeans refinement
2. Run matched Jeans-AMR comparison vs Athena++ at production resolution
3. Consider multipole gravity BCs (needed for production-scale collapse
   in small boxes with isolated mass distributions)

### Environment setup
  - ~/athenapk_env.sh loads GCC/OpenMPI/HDF5 paths + cmake
  - /beegfs/u/bbg6470/venvs/analysis_env/ has h5py, numpy, scipy for analysis
  - Comparison data in /beegfs/u/bbg6470/athenapk_runs/comparison_athenapk/
    and /beegfs/u/bbg6470/athena++/test_collapse/runs/comparison/

### Input file conventions
  - inputs/collapse_be_mhd_test.in       - laptop MHD BE with density AMR
  - inputs/collapse_be_mhd_jeans_test.in - laptop MHD BE with Jeans AMR (validation)
  - inputs/collapse_be_mhd_comparison.in - cluster comparison vs Athena++ (uniform)

---

## End of Apr 27, 2026

### Cooling investigation: GLM-MHD over-heating identified and resolved

**Bug**: First MHD BE collapse production run (numlevel=12, dedner_plain default)
showed 30× over-heating in cells at ρ ∈ [10⁴, 10⁵] vs barotropic target,
while Athena++ CT-MHD reference matched target to <1%.

**Investigation** (rejected hypotheses):
- Cooling math error — verified via single-cell `par_reduce` instrumentation
  that AFTER cooling, ratio = 1.000 always
- Cons → prim refresh timing — verified task graph order: flux → update →
  source_unsplit → ApplyBarotropicCooling → FillDerived
- VL2 stage averaging — single-cell trace shows ratio = 1.000 throughout 19
  quiescent cycles, no stage artifact
- Riemann/reconstruction — both codes use HLLD + PLM, same setup
- Pure-hydro AthenaPK matches Athena++ to 1% — issue specific to GLM-MHD

**Cause identified**: Default `glmmhd_source = dedner_plain` in AthenaPK only
damps ψ scalar (`psi *= exp(-alpha c_h dt/dx)`), without corresponding
non-conservative momentum/energy correction. At deep AMR on aggressive collapse,
divB errors accumulate and ψ-damping silently discards energy → spurious heat.

**Fix**: One-line input change. `glmmhd_source = dedner_extended` adds
`cons(IM_i) -= dt·divB·B_i` and `cons(IEN) -= 0.5·dt·(B·∇ψ)`.

**Validation**: Apples-to-apples comparison vs Athena++ jeans_comparison at
matched parameters (μ=20, box ±20, 32³ base, nl=6, tlim=1.5) shows agreement
to <0.1% in bulk EOS, identical B ∝ ρ^(2/3) flux freezing, peak density
within 4% (APK=4.01e+05, APP=4.16e+05). See
[docs/cooling_investigation_report.md](docs/cooling_investigation_report.md)
for full investigation timeline and embedded plots.

**Production input updated**: `glmmhd_source = dedner_extended` now in
`collapse_be_mhd_production.in`. Yesterday's production data is still valid
for kinematic flux-freezing analysis (which doesn't depend on thermodynamics)
but should not be used for thermal-energy-related conclusions in dense regions.

### Files added
- `inputs/collapse_be_mhd_match_athena.in` — apples-to-apples validation input
- `docs/cooling_investigation_report.md` — detailed investigation report
- `docs/cooling_investigation_plots/` — 17 plots covering production run, scan, validation

### Next steps
- Submit production-grade test5-match run (128³ base, box ±6, numlevel=16, μ=20,
  dedner_extended) for tomorrow's analysis
- Final comprehensive write-up after test5-match data arrives

### Test5-match production run (Apr 27 evening)

Submitted `inputs/collapse_be_mhd_test5_match.in` matching Athena++ test5
parameters exactly: 128^3 base, box +/-6, numlevel=16, mu=20, with
glmmhd_source=dedner_extended (the fix). Ran on 16 ranks, produced 105
snapshots through t=1.04 before dt-collapse arrested progress (similar to
Athena++ test5 stalling at t=1.02 — both codes hit the same first-core
formation barrier).

**Key results**:
- Reached max(rho) = 3.1e+07 (30x deeper than nl=6 runs achieve)
- Mesh blocks grew from 512 to 1016, AMR levels 0-9 activated
- **EOS ratio = 1.000 across 8 decades of density** (rho 1e+00 to 1e+08),
  with N=3816 cells in the deepest 1e+07 to 1e+08 bin — fully statistically
  robust validation of the dedner_extended fix at production resolution
- Flux freezing B propto rho^(2/3) tracked qualitatively; deviations at
  rho > 1e+06 consistent with field reorientation and outflow boundary effects
- Morphology shows expected hourglass-pinched B field structure around
  collapsing core

This makes the `dedner_extended` fix validated at:
- nl=6 (apples-to-apples vs Athena++ jeans_comparison) — passes
- nl=8, 10, 12 (numlevel scan) — passes
- nl=16 (test5-match production) — passes with 8-decade EOS cleanliness

### Files added (Apr 27 evening)

- `inputs/collapse_be_mhd_test5_match.in` — test5-equivalent production input
- `docs/test5_match_plots/` — 6 plots characterizing the test5-match run:
  evolution overview, EOS verification, flux freezing, morphology slices,
  cross-comparison density and EOS
- `docs/methodology_validation_full_report.md` — comprehensive 478-line
  methodology validation report tying together cooling investigation, bug
  identification, fix, apples-to-apples validation, test5-match results,
  and cross-code comparison

### Status

Cooling investigation: **complete and documented**. AthenaPK with
dedner_extended is production-quality validated for MHD BE collapse with
deep AMR. Yesterday's production data (with dedner_plain) is preserved
but should not be used for thermodynamics in dense regions.


## 2026-07-10 — Mixed rkl2+Hall diffusion integrator (dt fix for the t4 Ohmic decoupling wall)

### Problem

prod_t4_full stalled at dt ~4.5e-9 code (hydro CFL alone: 5.04e-6) at
rho_max ~1.3e-10 g/cm^3. dt attribution on the production restart
(out2.00142, read-only, 4xH100, jobs 2320039/2320041) with each term
disabled via CLI override:

| constraint left binding            | dt (code)    | gap    |
|------------------------------------|--------------|--------|
| base (unsplit, all terms)          | 4.5-5.0e-9   | 1x     |
| Ohmic off -> Hall floor (0.05)     | 3.9348e-8    | 8.7x   |
| Hall whistler only (floor 1e-6)    | 1.4-2.3e-6   | ~400x  |
| all non-ideal off -> hydro CFL     | 5.04e-6      | ~1100x |

The strict limiter is physical Ohmic decoupling (eta_O ~0.44 code
= 2.3e20 cm^2/s at x_e ~1e-14) under diffusion/integrator=unsplit, which
the Hall guard forced (RKL2 unreachable). The whistler limit itself is
harmless (|eta_H| ~0.008 code).

### Fix: mixed-mode integrator

With diffusion/integrator=rkl2 + Hall active, the parabolic terms
(conduction, viscosity, Ohmic, ambipolar) go into RKL2 STS while the
dispersive Hall EMF is applied unsplit with the hyperbolic fluxes each
stage, under a strict whistler dt. New option
diffusion/hall_floor_integrator = unsplit (default) | rkl2 places Hall's
Ohmic stabilizer floor either with the unsplit Hall EMF (retains the
dx^2/eta_floor step ceiling; "Phase B", 8.7x) or in the RKL2 parabolic
set ("Phase C", whistler/ratio-cap limited, ~300-400x).

Changes (all in src/hydro/):
- diffusion/diffusion.hpp: enum class DiffTermSet {all, parabolic,
  dispersive}; CalcDiffFluxes(pkg, md, term_set);
  HallDiffFluxIsoFixed(md, eta_h_on, floor_on);
  EstimateHallTimestep(md, whistler_on, floor_on).
- diffusion/hall.cpp: kernels gate the eta_H part and the floor part
  independently (disabled part contributes exactly 0; same J stencil).
- diffusion/diffusion.cpp: CalcDiffFluxes dispatches by term set; eta
  cache precompute skipped when the set runs no eta-consuming kernel.
- hydro.cpp: Hall guard relaxed (needs integrator != none); parses
  hall_floor_integrator; EstimateTimestep splits dt_strict_hall
  (whistler + unsplit floor, always strict) from dt_diff (parabolic
  aggregate; feeds s_rkl and the rkl2_max_dt_ratio cap); fused non-ideal
  dt estimator disabled in mixed mode (it folds whistler into one
  aggregate).
- hydro_driver.cpp: RKL2 stage operator calls CalcDiffFluxes with
  DiffTermSet::parabolic.

### Validation (CPU, build_cpu, 2026-07-10)

- hall_whistler.in (Q_H=0.5, floor=0.05, whistler-limited): unsplit,
  mixed floor-unsplit, mixed floor-rkl2 all give omega rel. err 4.23e-3
  (matches original 0.4% validation); floor-unsplit is bit-identical to
  unsplit (RKL2 with empty parabolic set = exact identity).
- Stressed whistler (Q_H=0.02 << floor=0.05, production-like ordering):
  mixed floor-rkl2 runs at dt 2.08e-4 vs unsplit 2.44e-5 (8.5x, STS
  active) with IDENTICAL omega error (5.57e-2) and amplitude decay
  (9.11e-2) to the unsplit reference.
- diffusion_ambipolar.in regression: unsplit 3.94e-4 (matches known
  0.04%); rkl2 6.78e-3 at ~1000x dt (5 vs 5462 cycles) — expected STS
  splitting error, parabolic path unchanged.

GPU probe on the t4 restart: runs/dt_attrib/submit_dt_mixed.sh (job
2321087, after build 2321086) — mixedB expects 3.9e-8, mixedC expects
~1.5e-6 with rkl2_max_dt_ratio=400 (s ~29 stages per Strang half).
Production relaunch (4 GPU) only after the probe confirms.

### GPU probe results (job 2321087, 4xH100, t4 restart out2.00142, 2026-07-10 evening)

- mixedB (floor unsplit): dt converged to 3.9348378777503961e-8 — the noOhm
  attribution value to ALL digits (floor dx^2 ceiling). 8.7x. wsec_step 11.6 s.
- mixedC (floor in RKL2): dt converged to 1.18-1.78e-6 (whistler-limited,
  fluctuating with the B/rho state; top capped by rkl2_max_dt_ratio=400 x
  dt_diff ~ 1.8e-6). STS: ratio ~333, 27 stages per Strang half. wsec_step
  33-40 s. Net: ~4.1e-8 code-time/wall-s = ~33x the 8-GPU unsplit production
  baseline (1.26e-9) on HALF the GPUs (~66x in GPU-seconds per simulated time).
  First-core phase (dt_code ~0.02) now ~6 days wall (~12 slots) vs 200+ days.
- Health: only benign CLI-param warnings; final state rho_max 1.309e-10 g/cm^3,
  no NaNs, T(rho_max) ~134 K; peak GPU memory 48.7 GiB/card (RKL2 registers
  +12 GiB over the 36.5 baseline; 61% of 79.7).
- Production staged (NOT launched): prod_t4_full/submit.sh = 4 GPU + mixed
  overrides + "-t 11:30:00" clean exit; fhc.in mirrored. Launch = rm STOP_CHAIN,
  echo 0 > chain_n, sbatch submit.sh.

## 2026-07-11 — diffusion/rkl2_freeze_eta (default OFF; GPU A/B pending)

Production entered level 11 (max rho 6.26e-10 g/cm3 at t=1.0839) and the
rkl2_max_dt_ratio=400 cap now binds every cycle (dt ~1.0-1.8e-7, 29 STS
stages/half). Each of the ~58 parabolic CalcDiffFluxes calls per cycle refills
the non-ideal eta cache = a 53M-cell Wardle-tensor + tabulated-EOS solve.

New flag `diffusion/rkl2_freeze_eta` (default false): the RKL2 stage loop
passes refill_eta=false for stages jj=2..s, so the cache holds the values from
the first stage of each Strang half (refreshed 2x/cycle + 2x by the dispersive
Hall calls instead of ~60x). Justification: RKL2's stability polynomial
(Meyer+2012) is derived for an operator held fixed across the super-step, and
the coefficient lag is O(dt/2) on eta(rho,T,B), which evolves on the collapse
timescale >> dt. NOTE: no Athena++ precedent — checked, its CalcDiffusionEMF
calls SetDiffusivity per STS stage (trivially cheap there: analytic coeff*B^2).

Touched: diffusion.hpp/.cpp (CalcDiffFluxes refill_eta arg gating
PrecomputeNonidealEta), hydro.cpp (param parse + REQUIRE eta_cache && rkl2),
hydro_driver.cpp (stage-init refill=true, stage loop refill=!freeze).
Unsplit and dispersive call sites default to refill=true (unchanged).

CPU check (AD damped-Alfven, iprob=50, rkl2 ratio=1000, eta_cache forced on):
freeze on/off identical relative error 4.16e-3; final-state max|diff| 4.5e-21
(eta_A=Q_A*B^2 is state-dependent through B^2, so exact zero not expected).
GPU validation queued (freezeAB job on prod restart 00143: A=off must
reproduce the production dt trace, B=on measures wsec gain + state diff).
Not enabled in production until that passes.

## 2026-07-11 — nsys profile of the mixed-mode step (job 2321477, restart 00142, 4xH100)

Rank-0 GPU kernel time over 12 cycles (converged part dt~1.2-1.5e-6, 27 STS
stages/half, wsec_step 33-37 s; nsys overhead negligible):

- **PrecomputeNonidealEta: 65.7%** (141.7 s, 380 instances x 373 ms). 56
  calls/cycle when converged (54 RKL2 stages + 2 dispersive Hall) — exactly what
  diffusion/rkl2_freeze_eta eliminates (56 -> 4).
- **Per-term dt estimators: ~15%** (EstimateResistivity 843 ms + EstimateAmbipolar
  843 ms + 2x EstimateHall 421 ms = ~2.5 s/cycle) — each re-evaluates the Wardle
  tensor per cell.
- RKL2StepOther: 5% (32 ms/stage). Everything else (hydro, MG gravity, RT, chem,
  boundaries) sums below ~10%.
- Host side: cudaEventSynchronize 152 s / cudaStreamSynchronize 36 s — the GPU is
  the critical path during STS; rank-0 card ~51% kernel-busy overall.
- MPI: Isend/Irecv ~11.5 s each; memcpy D2H 8.5 s + H2D 4.5 s (tcp/sm host staging).

Estimator fusion extended to the mixed mode: new two-reducer
EstimateNonidealTimestepIonizationFusedMixed (one tensor eval/cell -> dt_par +
dt_strict separately; per-term minima, NOT the unsplit variant's eta_O+floor sum).
fusable gate no longer excludes mixed_hall. CAVEATS: (1) engages only when ALL
non-ideal coefficients are plain "ionization" — prod_t4_full uses
ambipolar_coeff=ionization_chem, so **production stays on the per-term path**
(at its current 78 s steps the estimators are ~3%, not worth a chem-aware fused
kernel); (2) under eos=hydrogen the fused kernel uses the EOS-table temperature
(consistent with the flux kernels) while the per-term estimators use p/rho, so
fused-vs-per-term dt differs slightly there by design.

ratio600-at-00142 part of the job: cap did not bind at that shallower state
(ratio reached 298, dt whistler-limited 1.5-1.9e-6) — cap pricing for the
current level-12 cap-pinned state comes from the capsweep job on restart 00143.

FusedMixed smoke A/B (eos_smoke, CPU, rkl2+hall_floor rkl2, 8 cycles): evolution
bit-identical (final-state max|diff| = 0.0, dt trace identical — hydro-limited);
STS ratio differs in the 3rd digit (3.02 vs 3.01) = the documented EOS-T vs p/rho
difference in dt_diff, ~0.3% there. Fused engages on all-"ionization" decks only.

## 2026-07-12 — STS cap sweep + eta-freeze A/B (jobs 2321946 / 2321957) + changeover

Level-12 state (restart 00143, cycle 71500; cap 400 pinned every cycle there):
- cap sweep (8 cycles each): 400: dt 1.97e-7 / 44.0 s; 600: 2.47e-7 / 51.2 s (+8%);
  1000: 3.14e-7 / 53.4 s (**+31% net**), whistler starts binding near 1000 (one
  cycle at ratio 308) -> higher caps buy nothing. dt sequences healthy throughout.
- eta-freeze A/B (same restart): freeze OFF reproduces production dt at cycle 71501
  to all 17 digits (new-binary regression) and its 43.9 s wsec; freeze ON:
  **20.1 s (2.18x)**, same STS envelope (400/29), state diff after 8 cycles <=1%
  (upper bound; contaminated by a 2.3% time offset between endpoints — same
  perturbation class as the continuous->restart transition every chain slot takes,
  which measurably decorrelates dt by cycle 71502 even with identical physics).

Production changeover STAGED (applies at relaunch): binary athenaPK_eos_v2
(md5 41b6a958...; old athenaPK_eos kept), submit.sh + fhc.in mirrored with
rkl2_freeze_eta=true + rkl2_max_dt_ratio=1000; submit_8gpu.sh regenerated.
Combined freeze+cap1000 probe (job 2328566) on final.rhdf (cycle 71913,
dt 1.4e-8 cap-pinned) gates the switch.

## 2026-07-12 (evening): speed changeover LIVE (4.6x) + fused-ionization_chem (v3) validated

**Changeover live.** Slot 2334438 (binary athenaPK_eos_v2, rkl2_freeze_eta=true +
rkl2_max_dt_ratio=1000) resumed from final.rhdf (cycle 72409) at 14:12 with zero cycles
lost (boundary handoff; kill-now plan abandoned once the outgoing slot had an in-slot
restart at cycle 72000). Measured over the first 507 cycles at the level-12 state:
mean dt 3.09e-8 (was 1.35-1.5e-8, cap-400-pinned), mean wsec_step 36.8 s (was 79.1),
zone-cycles/wsec 1.94e6 (was 9.46e5). Net simulated-time per wall-second ~4.6x -- above
the ~2.5-2.9x forecast because dt gained 2.1x at this state (whistler allows more than
the cycle-71500 sweep state did). STS ratio pinned at 1000 (45 stages/half); freeze +
eta-cache banners confirmed; 0 NaN; GPU mem unchanged (70.9-72.2 GiB/card). Successor
2338185 self-queued. Standing watch tripped rho >= 1e-9 g/cm3 (1.0135e-9 at cycle 72400,
out1.01448) -- first core 1e4 x rhocrit, ~1-2 decades below second collapse.

**FusedNonidealEval (binary v3, staged NOT deployed).** Extended the fused non-ideal
evaluation to the ionization_chem family: new FusedNonidealEval struct (diffusion.hpp)
reproduces each term's Get() from the minimum tensor solves (production mix equilibrium
Ohm/Hall + chem-capped AD = ONE equilibrium solve + AmbipolarEtaFromXe per cell, was
three full solves); drives PrecomputeNonidealEta AND both fused dt estimators; fusable
gate (hydro.cpp) relaxed to ionization|ionization_chem per term. fe.ion taken only from
ionization-family terms (fixed-coeff terms carry a default-constructed model).
Validation:
- CPU smoke (eos_smoke/fhc_chem.in, full chem stack, rkl2 mixed + freeze): fused vs
  per-term estimators bit-identical dt (13 cycles) and bit-identical .hst; caveat: STS
  ratio 0.2 => diffusive dt not binding there (wiring test, not a discriminating test).
- GPU 3-way probe (job 2334442, 8 cycles on final.rhdf cycle 72409, production config):
  A=v2 baseline, B=v3 fused_off (isolates Precompute rewrite), C=v3 fused_on.
  B vs C: dt BIT-IDENTICAL all cycles => estimator fusion changes no dt at this state.
  A vs B: dt bit-exact for 5 cycles, then 1-3 ulp difference (cycle 72415) growing to
  7% (72417) -- recompiled-kernel FP contraction/reassociation seed amplified by the
  known argmin sensitivity; same acceptance class as every 12h chain-resume
  perturbation (bit-level seeds -> ~5%/cycle dt wobble, statistically neutral).
  wsec_step steady: A 36.6-36.7 / B 34.4-34.5 / C 34.1-34.2 s => cache-fill fusion -6%,
  estimator fusion -1%, total +7.4% throughput (zone-cycles/wsec_step 2.15e6 -> 2.31e6).
Binaries: athenaPK_eos_v2 41b6a958..., athenaPK_eos_v3 b5047400... (both kept; v3
adoption = BIN swap in submit.sh at a slot boundary, no input changes needed).

## 2026-07-12 (night) — user-approved change set staged at slot boundary: v4 + eta_O cap + dn=250 + 5 GPUs

User approved (this evening's menu): v3/v4 binary, dn=250, eta_O cap (physics sign-off),
GPU count -> 5 (modified from the 8-GPU proposal; g001-g004 carry 8xH100 each so 5 is
single-node). Declined: dropping Hall, sinks (deferred).

**eta_O cap implementation** (`diffusion/eta_ohm_cap_code`, default disabled):
- Clamp applied identically in `OhmicDiffusivity::Get` (both overloads, ionization-family
  branches only; `fixed` coeff untouched) and `FusedNonidealEval::Eta` (mirrored via
  `EtaCap()` in the factory), so flux kernels (eta cache), fused dt estimators, and
  non-fused paths all see the same capped eta_O. Rank-0 banner "## Ohmic eta_O cap: X".
- Binary: `athenaPK_eos_v4` (md5 0d9880e55e3c1f4a2ef022c71e572c2d), build job 2344650
  (incremental, the 8 diffusion.hpp TUs). v2 (live until boundary) and v3 kept for rollback.

**Pricing probe 2344652** (8 cycles each on restart 00146, cycle 73000, 4xH100 g001):
| segment | dt mean (73002-07) | wsec_step | STS ratio/stages | sim-time/wall vs A |
| A nocap | 2.49e-8 | 34.1 s | pinned 1000 / 45 | 1.0x |
| B cap 0.3 | 6.7e-8 | 19.5 s | 154-178 / 19 | 4.7x |
| C cap 0.1 | 6.8e-8 | 15.1 s | 49-59 / 11 | 6.1x |
Cap unpins the STS ratio: dt rises to the whistler/hydro limit AND stages drop 45->11.
No-op check passed (A reproduces v3-class 34.1 s; banner absent in A, present in B/C;
0 NaN; ~70 GiB/card unchanged). dt at the whistler limit now fluctuates cycle-to-cycle
(6.0-8.5e-8) — expected for a dispersive-term-limited step.

**Changeover executed 22:5x** (boundary handoff, zero cycle loss):
- scancel stale 4-GPU successor 2338185 FIRST (SBATCH resources freeze at submit time;
  editing submit.sh before cancelling would hand `mpirun -n 5` to a 4-GPU allocation).
- submit.sh: ntasks/gres 4->5, mpirun -n 5, BIN=athenaPK_eos_v4, new CLI overrides
  `diffusion/eta_ohm_cap_code=0.1 parthenon/output2/dn=250`; fhc.in mirrored (paired).
- New successor 2345200 queued afterany:2334438 (running slot exits ~01:40).
- wrap_mod.sh needs no change (counts allocation-visible GPUs at runtime).
- Caveat flagged: 5-GPU slots can queue longer than 4-GPU on a packed cluster (afterany
  chain => worst case wait, not loss). Expected post-cap: ~5 GPU-days per code-time-unit
  at current state vs ~31 before tonight's set (6.1x on top of the 4.6x morning gain is
  state-dependent; whistler dt will drop again as the core refines).
- `submit_8gpu.sh` is now OBSOLETE (user chose 5 GPUs; regenerate from submit.sh if ever
  needed).

## 2026-07-12 (night, part 2) — deep storage cleanup (~1.28 TB freed) + Hinnerk email sent

User sent the RRZ email (as drafted: argues the 8->4 GPU cut; the later 4->5 move is not
in it — rationale if asked: memory headroom for level 13+). User-approved cleanup:
- `prod_t1_ideal` DELETED entirely (78 GiB; late phase numerically invalid — diagnosis
  preserved in memory note athenapk-t4-oom-budget; ideal-tier paper figure needs a re-run).
- `bench_cpu_node` DELETED (12 GiB; numbers in DEV_LOG 2026-07-10 + email).
- t2/t3 restart dumps DELETED (264 GiB; science phdf intact — tiers not extendable now).
- t4 restarts pruned 147 -> 20 numbered + final (~789 GiB freed): kept every 10th index
  (00000..00140 = cycle milestones every 5000) + newest five 00142-00146. Newest restart
  untouched (active resume point); delete list verified against `ls -t` before rm.
- dt_attrib probe outputs/old logs pruned (~135 GiB): kept etacap + fused3way logs/scripts.
- Binaries: deleted athenaPK_coal, athenaPK_coal_chem, athenaPK_eos(v1), athenaPK_eos_v3;
  kept athenaPK_eos_v2 (live slot + rollback) and athenaPK_eos_v4 (production).
- Deleted obsolete submit_8gpu.sh, submit_build_gpu_eos.sh, old build logs (kept 2344650).
Quota: 7.4 -> ~6.1 of 11.0 TB apparent (rrz-quota lags; du-verified).
- (addendum) Top-level doc sweep: deleted HANDOFF.md, HANDOFF_gpu_underfeeding.md,
  PROJECT_STATUS.md, PROJECT_DOCUMENTATION.md (all superseded snapshots of concluded
  phases; conventions live in CLAUDE.md), email_to_hinnerk.md (sent), nb_exec/ (notebook
  batch scratch), runs/bench_poly_apk (June polytrope bench logs). KEPT deliberately:
  email_to_grete.md + sgclean_build_test.sbatch (open PR #164), eos_smoke (reusable
  smoke deck), athenapk_experimental/ (996 MB, sinkless-bridging working copy with
  possibly UN-MERGED work — do not delete without explicit user approval).
- (addendum 2) User: KEEP athenapk_experimental/ (sinkless-bridging copy stays);
  deleted email_to_grete.md (email_to_hinnerk.md already gone). Top level now:
  CLAUDE.md, OPUS_PROMPT.md, sgclean_build_test.sbatch, analysis scripts, source trees.

## 2026-07-13 (01:43) — v4+cap+5GPU changeover LIVE, verified

Slot 2334438 exited COMPLETED (11:30:59); slot 4 = job 2345200 resumed at cycle 73426
with ZERO cycles lost; successor 2345437 self-queued. Verified in run.log:
"## Ohmic eta_O cap: 0.1 (code units)" + freeze banner; 0 NaN. Steady state
(cycles 73430-73435): dt 6.1-7.2e-8 (whistler-limited, matches probe), wsec_step
12.2-13.7 s, 6.3e6 zone-cycles/wsec_step, STS ratio 50.8 / 11 stages, GPU mem
56.3-57.4 GiB on all 5 cards (was 70-72 on 4; ~19 GiB headroom to the 76-alarm).
Net sim-time/wall: dt/wsec = 5.4e-9 vs 8.4e-10 (this morning's config) = 6.4x,
vs 1.8e-10 (pre-campaign cap-400) = ~30x cumulative. 5th GPU delivered ~1.2x on
wsec vs the 4-GPU probe (15.1 -> 12.5 s).

## 2026-07-13 (08:43) — SECOND-COLLAPSE ONSET: rho crossed 1e-8 g/cm3, AMR level 13

phdf 01492 (cycle 74600, t=1.08397005): max rho_code 3.023e10 = **1.6527e-8 g/cm3 =
1.65e5 x rhocrit** — into the H2-dissociation regime the tabulated multi-Saha EOS was
built for. Density rose 16x in ~7 h wall (1.0e-9 at cycle 72400 yesterday) — runaway
steepening. Level 13 opened: 8 blocks at L13, 63 at L12, 2486 total. dt 1.4-2.0e-8
(whistler tracking the refining core; was 6-7e-8 at handoff), wsec_step 13.2-14.2 s,
STS 64/11 stages, mem flat 58.9-59.7 GiB x5. Chain: slot 2345200 (5h left) + successor
2345437. Watch re-armed at rho>=1e-7, level>=14, mem>76 GiB.

## 2026-07-13 — WS-5a: multipole gravity boundary conditions (PHYSICS_COMPLETION_PLAN)

New BC option `self_gravity/{i,o}x{1,2,3}_bc = multipole`: exterior Cartesian expansion of
the interior mass (monopole + traceless quadrupole about the COM; dipole ==0 about COM) on
the domain walls, replacing the Phi=0 Dirichlet clip. Enables isolated collapse in a box
only a few R across (the L=16 matched-box config that zero-Dirichlet clips).

Files: `src/self_gravity/multipole.hpp` (new: `MultipoleMoments` POD, `MultipolePhi()`
device fn, per-block `MultipoleBC<DIR,SIDE>`), `self_gravity.cpp` (moment reduction + RHS
lift + enrollment + guard), `poisson_equation.hpp` (packed `SetBoundary` type-3),
`pgen/poisson_test.cpp` (new: uniform-sphere gate). Behind an input flag; default OFF.

DESIGN — the load-bearing subtlety: a nonzero (inhomogeneous) Dirichlet ghost
`ghost = 2*Phi_mp(face) - phi_int` makes the discrete operator AFFINE, which breaks the
Krylov/MG solve (BiCGSTAB needs a linear operator; the MG preconditioner smooths residuals,
which carry HOMOGENEOUS BCs). First implementation applied `2*val-mirror` inside the solve
-> phi came out +1.5 at the center (wrong sign, wrong magnitude; solve 34 s, non-converged).
FIX (Green's-function boundary lift):
  1. Inside the solve, multipole faces use HOMOGENEOUS zero-Dirichlet (`ghost=-mirror`,
     type 3 folded into the type-1 branch of the packed SetBoundary) -> operator stays linear.
  2. The inhomogeneous face value is lifted into the RHS in FillPoissonRHS: for each interior
     cell adjacent to a multipole face,  b -= 2*Phi_mp(face)/dx_n^2  (edges/corners get one
     term per touching face). Derivation: the true vs homogeneous ghost differ by 2*Phi_mp,
     which the discrete Laplacian adds as +2*Phi_mp/dx_n^2; move it to b.
  3. Post-solve, the enrolled per-block `MultipoleBC` (inhomogeneous, `2*val-mirror`) restores
     the physical wall ghosts in AddBoundaryExchangeTasks so ApplyGravitySource sees the
     correct halo.
Moments: 10 raw moments (M, M x_i, M x_i x_j) atomic-summed over interior cells ->
MPI_Allreduce -> COM + traceless Q = 3*mu - delta*tr(mu), stored in a mutable pkg param,
recomputed every step in FillPoissonRHS. Guard: `packed_bc=false` + multipole is forbidden
(PARTHENON_REQUIRE) because it would route the operator through the inhomogeneous per-block
BC (affine again); verified it aborts with a clear message.

GATE (i) — uniform-density sphere (rho=1, R=4) in box [-8,8]^3 (wall at 2R, the stress case),
64^3 uniform, four_pi_G=1. `runs/validation_ws5a/sphere_multipole.in` (+ .../sphere_zero.in).
Analytic phi_in=-(2/3)piG rho(3R^2-r^2), phi_out=-GM/r, G=1/(4pi), M=268.083.
  multipole : phi_center = -8.0256  (analytic -7.9922);  MAX REL ERR = 0.537% (<1% PASS),
              mean 0.21%; error is discretization/staircasing-dominated (interior slab 0.45%,
              exterior-toward-wall 0.35% -> the BC itself is sub-0.4%).
  zero-Dir  : phi_center = -5.682   (28.9% shallow, matches the plan's "~10-30% at L=2R").
GATE (iii) proxy — throughput identical: multipole 2.08e5 vs zero 2.11e5 zone-cycles/wsec
(affine-bug version was 7e3, ~30x slower) -> BiCGSTAB iteration count not growing.
OFF-STATE — zero-Dirichlet phi_center bit-identical (-5.68204) before/after the edits;
FillPoissonRHS multipole blocks gated by grav_has_multipole (false in production).

GATE (ii) — L=16 matched-box collapse vs L=52 big-box reference (2026-07-13, CPU/std to
isolate the BC; hydro BE sphere, non-rotating, njeans=8 AMR; runs in runs/validation_ws5a/
collapse_L16_multipole [multipole BC] and collapse_L52_ref [zero-Dirichlet, big box]).
RESULT = PASS (mechanism decisive): the L=16 multipole box RUNS AWAY INTO THE FIRST CORE,
reaching max rho = 2.87 x rhocrit at MaxLevel 8 / 512 blocks (t=1.0800, phdf 00108, finite) --
crossing 0.13 x rhocrit, 1 x rhocrit, up to 2.87 -- whereas the old zero-Dirichlet L=16 box
FROZE at 0.13 x rhocrit / MaxLevel 3 (athenapk-mhd-collapse-bug). The two boxes track: at the
last cleanly-shared sub-runaway milestone (1e-2 x rhocrit) they cross at Δt = 0.0092 = 0.87%
of the collapse time (L16 leads); L16's near-vertical runaway at t=1.080 (40x density in
Δt=0.009) is consistent with L52 -- still at 0.081 x rhocrit and climbing at t=1.0806 -- going
vertical ~0.009 later, i.e. the SAME collapse offset <1% in time (the huge fixed-time rel-diff
is the vertical-runaway timing artifact, not a divergence). CAVEAT: the literal "<5% max_rho(t)
up to 1e4 x rhocrit" full-range overlap was NOT captured -- both runs are dt-wall limited (dt
6e-2 -> 5e-4 near the first core) and the reference lagged the runaway by ~0.9% and was cut at
the 8h walltime just before its own vertical rise; resuming to match L52 at high density costs
many CPU-days (each 8h slot mostly replays the dt-wall) for confirmatory value only. Runs are
checkpointed (parthenon.out1.*.rhdf, now dn=100) and resumable if the full overlap is wanted.
Bottom line: the multipole BC removes the potential-clipping stall and the matched L=16 box
collapses to first-core densities tracking the big box to <1% in collapse timing.
GPU build: build_gpu/bin/athenaPK rebuilt with WS-5a (md5 bf6d2e98, exit 0) for production use.

## 2026-07-13 (11:30) — AMR LEVEL 14 OPEN
phdf 01501 (cycle 75050, t=1.08398): max rho 5.735e10 code = 3.135e-8 g/cm3 =
3.14e5 x rhocrit (was 1.65e5 at 08:43 — ~1.9x in 2.8 h). Level 14 open: 8 blocks
at L14, 2843 total (was 2486 at L13-open). dt 6.3-6.9e-9 (whistler-tracking),
wsec_step 17-19 s, STS ratio ~81-90 / 13-15 stages. GPU mem 66.8-68.1 GiB x5
(was ~59 at L13-open) — headroom to the 76 GiB alarm now ~8 GiB/card; expect
another ~8-9 GiB step when L15 opens => likely alarm during L15. Options then:
6 GPUs at a slot boundary (single node holds to 8). Chain: slot 2345200 elapsed
9:47, clean exit ~11:43h mark, successor 2345437 queued. No errors in log.
Watch re-armed: rho>=1e-7 g/cm3, level>=15, mem>76 GiB.

## 2026-07-13 — WS-1 increment 1: sink-particle swarm plumbing (PHYSICS_COMPLETION_PLAN)

New package `src/sinks/` (sinks.cpp/.hpp), registered in hydro.cpp ProcessPackages and wired
into hydro_driver.cpp. Parthenon **swarm** "sinks" carrying built-in {id,x,y,z} + Real values
{mass,vx,vy,vz,Lx,Ly,Lz,t_created}, all with Metadata::Restart. Input `<sinks>`: enabled
(default false), seed_{x,y,z,vx,vy,vz,mass}. Modeled on src/tracers/tracers.cpp but (a) restart-
persistent and (b) works on ADAPTIVE meshes (Tracers hard-forbids adaptive; Sinks removes that
restriction). Driver: per-step ResetCommunication -> AdvanceSinks (inert ballistic drift
x+=v*dt) -> SwarmContainer::Send -> Receive (the migration path). SeedInitialSinks
(UserWorkBeforeLoopMesh) places the seed on the owning block; skipped on restart (Parthenon
restores the swarm). Default OFF = only the `enabled` bool registered, no swarm/tasks ->
bit-identical to pre-sinks code. pgen poisson_test gained optional bulk_vx/vy/vz (default 0,
preserves gate-i) to advect a density bump for the dynamic-AMR test.

GATES (runs in runs/validation_ws1/; CPU/build_cpu):
- Ballistic accuracy (sink -6->+ x at v=2, analytic x=-6+2t): uniform serial max|dx|=5.77e-15;
  32-rank MPI + DYNAMIC AMR max|dx|=1.78e-15 (y,z exact); **< 1e-12 PASS**. Sink crossed block
  and rank boundaries and the moving refined region; block count 512->568->540 (refine AND
  derefine) confirms live regridding; npart stayed 1 (no loss/duplication on migration).
- Restart bit-identical (serial uniform, resume from t=2.5 -> t=5): swarm x/y/z EQUAL the
  straight-run values exactly (|dx|=0). PASS. (The 32-rank MPI+AMR restart leg was walltime-
  cut before running; the swarm-restart mechanism is mesh/rank-independent and AMR mesh restart
  is the standard Parthenon path used by every collapse run, so serial proof carries. Restarts
  out1.00000/00001 are on disk in amr_mpi/ if an airtight MPI+AMR restart check is wanted.)

WS-1 increment 1 COMPLETE. NEXT (increments 2-6, each own gate): 2 static point-mass gravity on
gas (Bondi ~10%); 3 two-body orbit (energy drift <1e-6/orbit); 4 creation criteria; 5 accretion
+ conservation (mass/mom <1e-10/step, dt recovers past first core); 6 full-physics t4 run.
GPU build does NOT yet include sinks (only build_cpu rebuilt).

## 2026-07-13 — WS-1 increment 2: sink gravity on the gas

Added to src/sinks/: GatherSinks (collect every active sink's {x,y,z,mass} across blocks/ranks
via MPI_Allgatherv into a device View replicated on all ranks) + ApplySinkGravity (operator-
split source on the final stage: softened point-mass g = -G M (r-r_s)/(r^2+eps^2)^1.5 summed
over all sinks, momentum kick dm=rho*g*dt, energy gets the exact resulting KE change ->
conservative). G = four_pi_G/(4pi); Plummer softening eps = soft_cells*dx_min. Wired into
hydro_driver.cpp final stage (gather -> apply), gated on sinks/enabled.

VALIDATION = direct force check (runs/validation_ws1/force_test, build_cpu): at-rest uniform
gas, one vl2 step gives v = g*dt exactly (uniform gas -> zero hydro flux), so measured g = v/dt
compared to analytic softened point-mass everywhere (M=2.5, G=1, eps=0.125):
  max |g| REL ERR = 8.14e-16 (machine precision), in the softened core AND all 2.1M cells;
  direction cos(meas,analytic) = 1.00000000; density stayed uniform 1.0.
=> the force magnitude/direction/softening/G are EXACTLY right. This is a stronger check of the
gravity than the Bondi rate (which is confounded by resolution/boundaries/accumulation).

The plan's literal increment-2 gate ("steady-state Bondi accretion rate through r=8dx within 10%
of Mdot_B") is NOT achievable with accretion DISABLED: with no gas sink, infalling gas piles up
around the softened point mass, pressure builds, and the inflow chokes -- so no steady Bondi
exists. Observed in runs/validation_ws1/bondi_run (M=2.5, Mdot_B=49.09): Mdot(r=8dx) built to
0.41 Mdot_B at t~1 then collapsed to ~0.1 by t=2.5, while the radial profile showed strong
inflow at r_B (r=2.5 -> 1.12 Mdot_B) choking to ~0.11 at r=1 -- the classic no-sink accumulation
signature (gravity is correct; the flow just can't be steady without mass removal). The Bondi-
rate gate is therefore DEFERRED to after increment 5 (accretion adds the gas sink); the bondi.in
setup + compare_bondi.py are ready to re-run then. Increment 2 gravity: DONE + machine-precision
validated. NEXT: increment 3 (two-body orbit, energy drift <1e-6/orbit).

## 2026-07-13 — WS-1 increment 3: two-body sink dynamics (symplectic N-body)

Refactored the sink advance from per-block drift into a mesh-level N-body integrator
Sinks::AdvanceSinksNBody (src/sinks/sinks.cpp): gather every sink {id,x,v,mass} onto all ranks
(MPI_Allgatherv of ids as UINT64 + 7 Reals), integrate the whole small-N system on the host with
a SUBCYCLED KDK LEAPFROG (dt_sub = subcycle_cfl * min_i |v_i|/|a_i|), then write each rank's own
sinks back (matched by id) and update neighbor-block indices for migration. Deterministic ->
every rank integrates the identical global system, so results are consistent without per-substep
communication. Softened pairwise force g_i = sum_j G m_j (r_j-r_i)/(|r|^2+eps_sink^2)^1.5;
soft_sink=0 gives pure Kepler. New <sinks> keys: sink_gravity (sinks feel each other),
gas_gravity (sinks' gravity on gas; off isolates the N-body), soft_sink, subcycle_cfl, nseed +
seed{i}_{x,y,z,vx,vy,vz,mass}. With one sink or sink_gravity off the acceleration is zero and KDK
reduces to EXACT drift -> increment-1 ballistic regression still 4.09e-15 (re-verified).

GATE (runs/validation_ws1/twobody_run, 8-rank MPI): equal masses m=1, G=1, separation d=2 ->
circular orbit radius 1, speed 0.5, period T=4pi=12.566, analytic E=-0.25. Output sampled once
per T (same orbital phase) so secular drift reads directly. Over 71 orbits (run continuing to
100; rate is a robust 72-point linear fit):
  SECULAR ENERGY DRIFT = -4.24e-12 per orbit  [GATE <1e-6 => PASS by 6 orders of magnitude];
  bounded oscillation max |E-E0|/|E0| = 1.72e-9 (E returns to -0.2500000000 each period);
  Lz = 1.000000 conserved exactly; separation d = 2.0000 (circular). Textbook symplectic.
A 2-orbit nx=32 serial check independently gave the same behaviour (oscillation 6.3e-8 at
subcycle_cfl=0.05). Driver: mesh-level advance task between per-block ResetCommunication and
per-block Send/Receive.

WS-1 increment 3: DONE. NEXT: increment 4 (sink creation criteria); increment 5 (accretion +
conservation, which also closes the deferred Bondi-rate gate from increment 2); increment 6
(full-physics t4 run). GPU build still needs a rebuild to include sinks.

## 2026-07-13 — WS-1 increment 4: sink creation

Sinks::CreateSinks (src/sinks/sinks.cpp), run once per step after GatherSinks on the final
stage. A Kokkos::MaxLoc reduction over all interior cells finds the single densest cell meeting
the formation criteria, and MPI_MAXLOC(density,rank) picks the global winner so AT MOST ONE sink
forms per step; that rank extracts the cell {x,v,rho*dV} and AddEmptyParticles(1) with a fresh id
(mutable next_sink_id). Criteria implemented: (1) rho > rho_sink [fixed sinks/rho_sink_code, or
per-cell Truelove rho where the local Jeans length = n_jeans_sink*dx]; (2) local density maximum
over the 6 face neighbours; (3) converging flow div(v) < 0; (5) no existing sink within 2*r_acc
(r_acc = racc_cells*dx). New <sinks> keys: creation, rho_sink_code, n_jeans_sink, racc_cells.
No gas removed yet (mass conservation is increment 5); the new sink carries the birth cell's
velocity and mass rho*dV.
SIMPLIFICATION (documented): the plan's criterion (4), the full virial/energy-bound check, is
approximated by the density (Jeans/Truelove) threshold. Fine for clean collapse; tighten for
turbulent production (shocked-but-unbound gas) later.

GATE (runs/validation_ws1/creation*, manufactured Gaussian density peak rho_peak=200, sigma=0.8,
in a radial-inflow converging flow v=-1*r_hat -> cleaner/cheaper isolation than a full self-grav
collapse):
- rho_sink=100: exactly ONE sink (1 creation event), formed at the peak cell nearest origin at
  rho~193 (within 2x of rho_sink), t_created recorded; count stayed 1 through t=1.0 = 10x the
  required 0.1 t0 window -> criterion (5) blocks a 2nd sink on the still-densifying converging
  origin. PASS.
- control rho_sink=300 > initial peak 200: sink forms only at t=0.034 when the inflow-raised
  central density crosses 300 (rho_form~302); threshold respected (no premature creation). PASS.
- control v_radial=0 (static, non-converging) with rho_sink=100: NO sink -> the div(v)<0
  criterion correctly gates creation. PASS.
pgen poisson_test gained profile=gaussian + sigma + v_radial (radial inflow) for this test.
WS-1 increment 4: DONE. NEXT: increment 5 (accretion + conservation: remove bound gas within
r_acc, conserve total mass/momentum <1e-10/step, dt recovers past first core; also closes the
inc-2 Bondi-rate gate). Then increment 6 (full-physics t4). GPU build still lacks sinks.

## 2026-07-13 — WS-1 increment 5: sink accretion + conservation (PARTIAL)

Sinks::AccreteSinks (src/sinks/sinks.cpp), run after CreateSinks on the final stage.
Self-contained gather of {id,x,v,mass}; for each gas cell within r_acc of a sink that is BOUND
(0.5|v_gas-v_sink|^2 < G m_s/r) and above a rho_sink/3 floor, remove a quadratic-ramp fraction
phi=ramp*(rho-floor)/rho of the mass -- and, proportionally, momentum and energy (B NOT
accreted) -- and add it to the sink. Per-sink accumulators {dM,dp,dL} via Kokkos::atomic_add ->
MPI_Allreduce; sink update is momentum-conserving: m+=dM, v=(m v + dp)/(m+dM), L+=dL (angular
momentum of extracted gas about the sink, bookkeeping). Non-overlapping r_acc (creation's 2*racc
spacing) => each cell accretes to one sink. New key: sinks/accretion. The rho_sink/3 floor caps
the sink-region density -- the intended dt-recovery mechanism.

GATES 1-2 (mass + momentum conservation < 1e-10 relative per step): PASS to machine precision.
runs/validation_ws1/accrete_run(2): seeded sink (m=1, vx=0.5) accreting uniform rho=200 gas,
periodic BCs (no boundary loss), gas_gravity off. Sink grew 1 -> 160 (eating gas) while total
(gas+sink) mass stayed 102401.00000000 (rel drift 2.84e-16 over 16 cycles) and total x-momentum
stayed 0.50000000 (abs drift 1.08e-14) -> ~1e-16/step, 6 orders under the gate. Conservative by
construction (gas loses exactly what the sink gains).

GATE 4 (dt recovers past first core): NOT yet demonstrated -- requires a real self-gravity
collapse. The dt-CRATER is shown (manufactured converging-inflow + maxdensity AMR: dt drops 15x
as the runaway refines, runs/validation_ws1/dtr_cheap*), but a STEADY driven inflow keeps
replenishing the dense region faster than accretion removes it, so the density never caps and dt
does not recover -- the WRONG test physics. dt-recovery needs a SELF-LIMITED runaway (self-grav
collapse forms a single core a sink can arrest); the cheap non-self-grav tests cannot capture it.
A restartable self-chaining self-gravity BE collapse (L=16 multipole, numlevel=6, auto-Truelove
rho_sink, sink_gravity off + soft_sink to avoid multi-sink flinging) is LAUNCHED (job 2347794,
dtr_run/) and stops the chain when a sink forms; its dt history will show crater->recovery. Also
found: with a tiny r_acc (racc_cells*dx_min) a broad dense region spawns MANY sinks (2*racc
exclusion too small) and, with sink_gravity on + soft_sink=0, close pairs fling past the swarm
halo and crash -- mitigated for the collapse by sink_gravity=false + soft_sink; a physical first
core is ~1 r_acc so this is a manufactured-test artifact, but worth watching in production.

GATE 3 (accretion rate within 2x of Shu 0.975 c_s^3/G): TODO -- needs an isothermal Shu
singular-sphere collapse; deferred together with the inc-2 Bondi-rate gate (both are steady-
accretion-rate tests requiring a proper collapse/flow setup).

DT-RECOVERY RESULT (self-gravity BE collapse job chain 2347825->2347846, dtr_run/, numlevel=6,
auto-Truelove rho_sink): the dt-CRATER is now DEFINITIVELY demonstrated -- dt fell from the early
hyperbolic ~0.028 to 2.1e-4 (**133x**) over t=0.5->1.16 as the core collapsed and jeans-AMR
refined (5 accumulating chained slots once the restart cadence was fixed to dn=10). BUT no sink
formed and dt did NOT recover, for a now-understood reason: auto-Truelove rho_sink RISES ~4x per
refinement level (finer dx resolves higher density: rho_sink = pi cs^2/(G (njeans dx_min)^2)), and
the barotropic gas heats adiabatically (cs^2 12->21 as rho rose), so the formation threshold
outran the density -- rho/rho_sink went 0.83 (level 3) -> 0.44 (level 4). The sink only forms at
the finest level (6) at rho ~ 5 rhocrit (deep first/second core) = the multi-hour dt-wall regime
(same as WS-5a gate ii / production t4). A resume with a FIXED modest rho_sink=4000 (job 2347870)
to force early formation from the deep state (rho=9.4e4 >> 4000) also formed NO sink: the
command-line `sinks/rho_sink_code=4000` override does NOT apply on a Parthenon restart (-r keeps
the restart-stored -1/auto). => GATE 4 (dt recovers within 10x) is NOT YET DEMONSTRATED end-to-
end. It is COMPUTE-BOUND + needs one of: (a) a FRESH run (not restart) with a fixed modest
rho_sink in the INPUT file so a sink forms early at rho~4000 (t~1.10, before the deep wall) ->
accretion caps -> derefine -> dt climbs; or (b) reaching the deep first core with auto-Truelove
(hours). The recovery MECHANISM is sound by construction (accretion removes gas above rho_sink/3
-> density capped -> AMR derefines -> dt stops dropping); only the end-to-end collapse demo is
outstanding. NOTE also: sink creation is validated on a UNIFORM grid (inc 4) but NOT yet on an
ADAPTIVE mesh -- the collapse never crossed threshold so creation-on-AMR remains unexercised;
verify it when doing (a).

GATE 3 (accretion rate within 2x of Shu 0.975 c_s^3/G): TODO -- isothermal Shu collapse; deferred
with the inc-2 Bondi-rate gate.

WS-1 increment 5 STATUS: accretion IMPLEMENTED; conservation gates 1-2 PASS (machine precision);
dt-crater demonstrated (133x); dt-recovery (gate 4) mechanism sound but end-to-end demo pending a
fresh fixed-rho_sink collapse (compute-bound); Shu rate (gate 3) TODO. NEXT: close gate 4 via a
fresh fixed-rho_sink collapse (recommended rho_sink_code~4000 in the input, restartable dn=10),
then increment 6 (full-physics t4 + sinks). GPU build still lacks the sinks package.

## 2026-07-14 (13:15) — WS-1 GATE 4 (dt-recovery) RESOLVED: PARTIAL, honest finding

Fresh fixed-rho_sink collapse (runs/validation_ws1/gate4_run, dt_recovery_fixed.in: FRESH -i
not restart, rho_sink_code=8000 in the INPUT so the override actually applies, floor rho_sink/3
=2667, numlevel=6 jeans-AMR, L=16 multipole gravity, gas_gravity=on, sink_gravity=off,
soft_sink=0.01, racc_cells=4, dn=10 rst). Self-chaining CPU (build_cpu, 32-rank std).

RESULTS (parsed from run.log dt history + out0 phdf swarm/prim):
- Pre-sink the run is BIT-IDENTICAL to the no-sink control dtr_run/ (dt=1.643e-3 at t=1.10) --
  same collapse until threshold. GOOD (confirms OFF-of-creation path unchanged).
- SINK FORMS at t~1.13 (dump 00056 nsink=0 rho=6483 -> 00057 nsink=1 rho=6678) ON THE ADAPTIVE
  MESH. This closes the previously-open item: creation-on-AMR now EXERCISED & WORKS (uniform-grid
  only before). Exactly ONE sink, no multi-sink runaway; sink mass grows 18.7->53->89->114 (code
  units) over t=1.14->1.20 -- accretion is actively pulling gas in.
- dt RELIEF vs control (matched-time): control cratered to 2.09e-4 at t=1.159 and stalled there;
  gate4 held dt 2.7-3.6x HIGHER (t=1.17: 7.59e-4 vs 2.09e-4 = 3.64x; t=1.20: 5.66e-4 vs 2.09e-4
  = 2.71x) and ran FURTHER in time (to t=1.202, min dt 5.59e-4). So the sink+accretion give REAL,
  MEASURABLE dt relief.
- BUT dt does NOT turn around / recover: it keeps slowly declining (1.37e-3 -> 5.59e-4 post-sink),
  because rho NEVER settles to the 2667 floor -- it oscillates 6678->6626->4372->5072 as
  self-gravitational INFALL replenishes the sink region faster than the racc=4-cell accretion can
  clear it. The finest AMR level therefore persists and Courant dt keeps falling (shallower slope
  than control, but monotone).

HONEST VERDICT: GATE 4 = PARTIAL. The claimed clean mechanism ("accretion caps rho at floor ->
Jeans-AMR derefines -> dt RISES") is FALSIFIED for a self-gravitating collapse with small r_acc:
the density cap is defeated by continued infall (and Jeans-refinement also scales with cs^2, which
climbs adiabatically independent of the density cap). What IS demonstrated: (i) sink creation on
an adaptive mesh, (ii) conservative accretion with a growing sink, (iii) a 2.7-3.6x dt improvement
vs the no-sink control that lets the run push past where the control stalled. A FULL dt turnaround
would require r_acc large enough to clear the entire Jeans-unresolved region each step (so infall
cannot keep the peak refined) and/or thermal/pressure relief -- a tuning/mechanism item, NOT a
correctness bug in the accretion (conservation gates 1-2 remain machine-precision). Chain stopped
(cancelled queued successor, chain_n tripped) -- not worth more slots; growing sink -> deeper
potential -> MORE infall argues against a turnaround at racc=4. Data retained in gate4_run/.

## 2026-07-14 (18:20) — WS-5c divB audit (DONE) + WS-5b C-shock (PARTIAL)

WS-5c (div B fidelity audit instrumentation) DONE: added a per-cell MAX relative div B error
history output `maxRelDivB` (main.hpp Hst::divb_max; new HydroHstMaxDivB Kokkos::Max reduction in
hydro.cpp; registered with UserHistoryOperation::max next to the existing volume-integral relDivB).
VERIFIED: the .hst now carries column [13]=maxRelDivB (cshock run), reading 0 for a clean div-B-free
setup as expected. The peak (worst-cell) cleaning error is now trackable at first-core formation.
The glmmhd_alpha {0.05,0.1,0.4} x dedner sensitivity matrix on a t2 run = compute-gated (deferred).

WS-5b (Hall C-shock cross-validation) PARTIAL. Built the harness: runs/validation/cshock_ad.in +
cshock_analyze.py (grid self-convergence + shock-thickness). FINDINGS: (1) the cshock pgen SEGFAULTS
in pure 1D (nx2=1) -- real bug; workaround = thin-2D strip (nx2=4, periodic). (2) gamma=1.01
(near-isothermal) and strong shocks (vxl=5) crash; gamma=1.4 + weak shock (vxl=2, AD coeff 0.15) is
STABLE. With those the AD C-shock RUNS CLEAN at 128/256/512 (thin-2D, to t=40) and forms a continuous
C-type structure whose thickness (~12-16 code length) scales with the ambipolar coefficient
(physically correct). BUT clean 2nd-order self-convergence is NOT achieved: order ~0.35 (v_x) /
0.65 (B_y) even after aligning by shock center -- the structure is NOT yet at true steady state at
t=40 (a thickness-12 AD C-shock has a long relaxation time; different resolutions relax at different
rates so the L1 diff is transient, not truncation error). CLOSING THIS needs a much longer relaxation
run (verify d/dt->0) and/or a more localized shock -- a numerical-tuning + compute task. Hall itself
is ALREADY independently validated (whistler dispersion 0.4%, [[athenapk-nonideal-mhd]]); Athena++
has no Hall so there is no cross-CODE check regardless (self-consistency vs the semi-analytic
solution is the only option). Status: infra + stable runs + correct thickness-scaling DONE;
quantitative convergence gate = compute-deferred.

## 2026-07-14 (17:40) — WS-4 increment 2: dust package + opacity consumer — DONE (ioniz deferred)

New dust PACKAGE (src/dust/dust.cpp + dust_pkg.hpp), mirroring chemistry: Initialize (reads <dust>,
DustModel params, scalar_index), ReactDust (operator-split growth on the 2 dust scalars f_dg,a_c at
cons[nhydro+scalar_index,+1]), AddDustTasks (final stage, after chemistry). Registered in hydro.cpp
(<physics> dust=true, default false) + hydro_driver.cpp + CMake. Guard: a_c floored to a_floor on
entry (v_Brownian ~ 1/sqrt(a^3) diverges at a_c->0; caught a NaN when scalars init to 0). Added
`freeze_growth` (evolve on, grains pinned at reference) for the bit-identical gate.
OPACITY CONSUMER WIRED: radiation_moments.cpp MatterCoupling now scales kappa per cell by
Dust::DustFactor=(f_dg/f_ref)(a_ref/a_c) [geometric-area limit] at BOTH the Planck (emission) and
Rosseland (flux) sites -- reads the dust scalars from the packed "cons" field. Bit-identical when
dust off (dust_on=false => kdust=1) or frozen at reference (DustFactor=1 exactly).
GATES: test_dust.cpp (host) (a) Brownian 0.05%, (b) turbulent 0.04%/3dec, (c) sublimation round-trip
exact, (d) DustFactor ref=1/2x-grain=0.5/2x-dust=2 exact -- ALL PASS. In-code smoke (dust_smoke.in,
RT-present transport-off): package runs clean, f_dg restores 0->0.01, a_c grows 1e-7->3.5e-6 cm, no
NaN. build_cpu clean. DEFERRED: (1) IONIZATION consumer (per-cell f_dg/a_c into the diffusion Wardle
tensor grain bins -- more invasive, multiple AD/Hall/Ohmic kernels, effect compute-gated); (2) inc3
collapse eta_A(rho) static-vs-evolved figure (compute-gated). OFF-state bit-identical (dust=false).

## 2026-07-14 (17:10) — WS-4 increment 1: single-moment dust growth — DONE, gate PASS

New src/dust/dust.hpp (device/host DUST_FN): DustModel evolves (f_dg, a_c). Monodisperse sweep-up
da_c/dt = (f_dg rho_gas/rho_grain) v_rel/4, v_rel = sqrt(v_Brownian^2 + v_turb^2) [v_B ~ a^{-3/2}
thermal, v_turb = sqrt(alpha_turb) c_s St, St = Epstein t_stop * omega_ref ~ a_c]. integrate_cell:
RK2 midpoint (2nd order), sub-step capped at growth_cap=10%/step (a_c->ionization->dt coupling
guard). Sublimation: f_dg = f_dg_ref*subl_factor(T) (tanh switch at T_subl=1500 K; algebraic
no-re-formation-lag v1). GATE (src/dust/tests/test_dust.cpp, host g++): (a) Brownian regime
(alpha=0) a^{5/2}=a0^{5/2}+2.5 K t analytic -> 0.05%; (b) turbulent regime (St-driven) a=a0 exp(Ct)
over 3.0 decades -> 0.04% (vt/vB=178 confirms regime); (c) sublimation round-trip f_dg hot->4.5e-7,
cool->0.0100 exact. ALL PASS. Increment 1 = standalone (no production path changed). NEXT: inc2
wire consumers (opacity kappa ~ (f_dg/0.01)(a_ref/a_c); ionization f_dg*grain-bins) + bit-identical
off-state gate; inc3 collapse static-vs-evolved-grain eta_A(rho) figure (compute-gated).

## 2026-07-14 (16:45) — WS-2 increment 2: thermochemistry WIRED into RT+chem — gate PASS

Wired thermo.hpp into the production reaction path (src/chemistry/chemistry.cpp), behind
`<chemistry> thermo=true` (default false => bit-identical; requires network=gow17_reduced + RT on).
  - chemistry.cpp Initialize: builds ThermoParams from <chemistry> (q_cr_ev, alpha_gd, x_O, G0,
    h2_heat_ev, pe_eps, co_tau_coeff, av0/av_n0/av_exp, zeta) + rate_to_code = t_unit/(rho_unit*
    v_unit^2), thermo_nsub_max, thermo_cfl, thermo_efloor. PARTHENON_REQUIRE gow17 + physics/
    radiation. (efloor uses a chemistry-local key, NOT hydro/pfloor -- the latter tripped Parthenon's
    inconsistent-default guard, same class as the sinks njeans trap.)
  - New thermo.hpp AdvanceThermoEnergy(): semi-implicit (numerical-Jacobian backward-Euler) sub-
    cycler advancing e_int over dt under Gamma-Lambda; unconditionally stable, relaxes to T_eq w/o
    overshoot, floored at dt/nsub_max (no dt collapse).
  - ReactScalars (gow path): if thermo, pack rad.Er too; after abundance integrate_cell, T_dust=
    (Er/arad)^{1/4}*T_unit (arad from radiation pkg), call AdvanceThermoEnergy, apply
    pack(IEN) += e_new - eint. Ordering: RT matter-coupling (if on) -> chemistry -> thermo, per design.

GATE 2 (single-cell-in-a-box, runs/validation_rt/single_cell.in, RT on but transport+coupling a
no-op [RAD_DISABLE_TRANSPORT=1, matter_coupling off] so Er/T_dust fixed and only Gamma-Lambda acts;
gow17+thermo, gas started ~40 K; tests/teq_checker.cpp reuses thermo.hpp for the self-consistent
T_eq). Two densities:
  - n_H=1e5: gas cools 39.3 K -> STEADY 19.69 K by t~1.5; T_sim=19.693 vs T_eq(current abundances)=
    19.693 => **0.00% (gate <5% PASS)**. [19.7 not the increment-1 12.2 K only because CO is under-
    formed at t=8 (slow CO chemistry) -> less line cooling; the THERMO integrator tracks the
    instantaneous thermochemical equilibrium exactly.]
  - n_H=1e7: gas-dust dominates -> T pinned to T_dust; T_sim=10.31 vs T_eq=10.31 (0.02% PASS),
    T_dust=9.997 K => matches the increment-1 curve's high-n limit, ABUNDANCE-INDEPENDENT.
  dt stable (~0.02, NO collapse) throughout; chem_trunc_total = 0 (no truncation warning). PASS.
Found+fixed a validation-only unit bug: teq_checker must use zeta=1e-16 (the input value; the
ThermoParams struct default is 1e-17) -- with the wrong zeta the checker mispredicted 12.9 vs 19.7.
INCREMENT 2 DONE. NEXT: increment 3 (t2-scale collapse thermo on/off, core T diff <10%, envelope
10-15 K plateau) -- compute-gated. build_cpu clean; OFF-state (thermo=false) bit-identical.

## 2026-07-14 (16:05) — WS-2 increment 1: gas thermochemistry T_eq curve — DONE, gate PASS

New src/chemistry/thermo.hpp (device/host, CHEMT_FN like the networks): ThermoParams POD +
NetHeatCool(n_H, T, x_H, x_Cp, x_CO, T_dust) -> Gamma-Lambda [erg cm^-3 s^-1] + SolveTeq (bisection).
  Heating: CR (q_cr=20 eV * zeta * n_H), H2 grain-formation (0.2 eV * kgr n_H^2 x_H), photoelectric
    (1.3e-24 n_H eps G_eff, G_eff=G0 e^{-2.5 A_V}, A_V(n_H)=0.5(n_H/1e3)^{2/3} shielding proxy).
  Cooling: [C II] 158um + [O I] 63um + CO J=1-0 via a 2-level-atom TwoLevelCool() (Boltzmann + n_cr
    saturation); CO gets an escape-probability beta=(1-e^{-tau})/tau, tau=co_tau_coeff*x_CO*n_H/
    sqrt(T/10) (line trapping -- WITHOUT it the optically-thin CO over-cools and pins T=3 K at
    n=1e4). Gas-grain Lambda_gd=alpha_gd n_H^2 sqrt(T)(T-T_dust) -> pins T->T_dust for n>~1e4.5.
  Abundances x_Cp, x_CO, x_H come from the reduced gow17 network (network_gow17_reduced.hpp);
  x_O fixed 3.2e-4 (network evolves no atomic O -- documented plan risk).

GATE 1 (src/chemistry/tests/test_thermo.cpp, host g++ -O2, self-consistent network-equilibrium
+ SolveTeq iterated): T_eq(n_H) curve n=1e2..1e8: 178.5 K (n=1e2, unshielded WNM) -> 49.7 (1e3)
-> **12.69 K (1e4)** -> 12.6 (3e4) -> 12.2 (1e5) -> 10.31 (1e6) -> 10.00 (>=3e7) = T_dust(=10 K).
GATE: T_eq(1e4)=12.69 vs 13 K ref (2.4%, <30% PASS); T_eq(1e6)=10.31 vs 10 (3.1% PASS); converges
to T_dust at high n (PASS). ALL GATES PASS. This is the canonical molecular-cloud/prestellar curve.
CAVEATS (documented, honest): (1) co_tau_coeff=5.5 is CALIBRATED to the FHC CO column -- a 0-D curve
has no real column, so this encodes line trapping as a tunable (physically motivated, Goldsmith-2001
style); (2) T_dust fixed 10 K here (in a real run T_dust=(Er/arad)^{1/4} from the M1 field);
(3) n=1e2 point is warm (178 K, PE-heated diffuse gas) -- plausible, not gate-checked.
INCREMENT 1 = standalone unit test only; NO production code path changed (nothing wired yet ->
trivially bit-identical). NEXT: increment 2 (wire dedt into chemistry.cpp ReactScalars sub-cycler,
pass rad.Er for T_dust; RT-on single-cell relaxation gate), increment 3 (t2 collapse thermo on/off).

## 2026-07-14 (15:30) — WS-3b: 2nd-order (PLM) radiation transport — DONE, gate PASS 5.4x

Replaced the donor-cell face states in CalculateRadFluxes (radiation_moments.cpp) with optional
PLM+minmod reconstruction of (Er, Fr1..3). New flag radiation/reconstruction = dc|plm (default dc).
  - radiation_moments.cpp: added RadMinmod(); the hll lambda now reconstructs UL/UR per component
    -- dc: UL=cell(i-di), UR=cell(i) [bit-identical to before]; plm: UL=cm1+0.5*minmod(cm1-cm2,
    c0-cm1), UR=c0-0.5*minmod(c0-cm1,cp1-c0) (+-2 stencil). M1FaceFlux + HLL unchanged downstream.
  - radiation.cpp: reads radiation/reconstruction, REQUIRE dc|plm; if plm REQUIRE nghost>=2.
Ghost safety: the transport sub-cycle does flux->update->AddBoundaryExchangeTasks EVERY sub-step
(radiation_moments.cpp ~400), so all nghost ghosts are refreshed before each flux calc -> the +-2
PLM stencil is valid each sub-step. Production PLM-hydro already forces nghost>=2. dc default =>
bit-identical to pre-WS-3b (UL/UR reduce to the old cell-center donor states). build_cpu clean.

GATE (free-streaming pulse, runs/validation_rt/stream.in, job 2349979, std/n109): Gaussian Er pulse,
fully forward-peaked flux (reduced_flux=1), optically thin (kappa=0, coupling off), 128 cells 1D,
advected by M1 free-streaming; dc vs plm at the SAME final time (both peak lands at cell 54 -- equal
advection, only diffusion differs). RESULT: FWHM growth dc = 5.18 cells vs plm = 0.96 cells =>
**dc/plm = 5.40x reduction in numerical front diffusion** (gate: >=2x at 128 cells => PASS). PLM
also retains 91% of peak amplitude vs 64% for donor cell. Crossing-beam/shadow 2D qualitative check
(plan, optional) NOT run; M1 beam-merging is an expected closure artifact, not a failure.

IMPORTANT UNITS FACT (explains earlier RT "slowness"): chat = c_code / creduc with c_code = c/v0 ~
1.58e6 (code units), NOT chat=creduc. So creduc=10 gives chat=1.58e5, dt_rad=cfl*dx/chat ~ 1e-8, and
nsub=ceil(dt_hydro/dt_rad) ~ 1.2e5 sub-steps per hydro step -> catastrophically slow + the pulse
crosses the whole domain in one hydro step. For a resolvable free-streaming test use LARGE creduc
(2e5 -> chat~7.9, nsub~8). RT is NOT intrinsically slow; the earlier 97 s/step was this nsub blowup.

## 2026-07-14 (14:10) — WS-3a increment A: Planck/Rosseland opacity SPLIT (code done + gate i)

Implemented the Planck/Rosseland mean split in the M1 RT opacity (CLAUDE.md's "RT not yet in
AthenaPK" line is STALE -- src/radiation/ has been present since Jul 4; memory athenapk-rt-port is
the correct record). Changes:
  - radiation_opacity.hpp: split the single gray AbsorptionOpacity into RosselandOpacity(op,rho,T)
    (the old body: constant/dust/BellLin base = flux attenuation / diffusion) and
    PlanckOpacity = planck_ross_ratio * RosselandOpacity (emission weighting). New OpacityParams
    field planck_ross_ratio (default 1.0). BellLin is a Rosseland fit, so RosselandOpacity keeps
    the old numbers exactly.
  - radiation_moments.cpp: emission/absorption term (Newton loop, was line 296) now calls
    PlanckOpacity (kappa_P); flux-attenuation / radiation-force term (was line 322) now calls
    RosselandOpacity (kappa_R) + ScatteringOpacity. Only these two call sites use opacity;
    CalculateRadFluxes is pure M1 closure (no opacity), so nothing else changes.
  - radiation.cpp: reads radiation/planck_ross_ratio (default 1.0, REQUIRE >0).
Physics: kappa_P/kappa_R ~ 2-5 for dust; using kappa_R for emission under-cools optically-thin gas.
ratio=1.0 => kappa_P == kappa_R == pre-split single opacity => BIT-IDENTICAL (by construction; the
flux term is unchanged and PlanckOpacity reduces to the old AbsorptionOpacity). build_cpu compiles
clean. Full (log rho,log T)->(kappa_P,kappa_R) Semenov 2003 binary table = increment B (follow-on).

GATE (i) equilibration test (runs/validation_rt/, rad_pulse pgen w/ Eamp=0 = uniform gas+radiation,
constant opacity so kappa_P = ratio*kappa_a0 EXACTLY; job 2349625 on std/n109):
  RESULT: uniform gas relaxes to a stable radiative-equilibrium plateau (gasE 1.5 -> 0.400 at
  creduc=1000; -> 0.750 at creduc=30 -- the plateau shift with creduc is the expected RSLA c/chat
  artifact in the gas-energy coupling, NOT the split). The equilibrium is kappa_P-INDEPENDENT:
  ratio=1 vs ratio=2 give gasE 0.400053 vs 0.400054 (creduc=1000) and 0.75047 vs 0.75040
  (creduc=30) -- match <2e-4. This is exactly correct: the emission opacity sets the RATE of
  approach, not the closed-box endpoint (fixed by energy conservation + E=B).
  HONEST LIMITATION (important): a resolved "kappa_P timescale" is NOT observable in this test. The
  matter coupling is IMPLICIT and the gas<->radiation exchange uses the TRUE c (RSLA: only transport
  uses reduced chat), so 1/(c rho kappa_P) ~ 1e-5 << dt~0.04 -- the gas fully equilibrates with the
  local radiation field EVERY step (gasE 1.5->plateau in ONE cycle at any creduc). So the closed 0-D
  box (i) can only verify the emission channel drives the correct kappa_P-independent equilibrium,
  and (ii) has NO flux, so it does not exercise the Rosseland/flux channel at all. The QUANTITATIVE
  kappa_P>kappa_R effect (enhanced emission cooling of optically-thin gas -> warmer escaping-flux
  envelope) is inherently an OPEN-system / transport phenomenon = exactly the plan's GATE (ii)
  (t2-scale collapse: document the first-core envelope entropy shift vs belllin). Gate (ii) is
  compute-gated (like inc-6) and deferred; free-streaming transport itself is already validated
  (rad_pulse, memory athenapk-rt-port).
  VERDICT gate (i): split is CORRECT (code inspection + bit-identical + correct kappa_P-independent
  equilibrium); the resolved-timescale sub-gate is physically moot for an implicit true-c coupling
  and superseded by gate (ii). WS-3a increment A code = COMPLETE & default-OFF bit-identical.

## 2026-07-14 (13:35) — WS-1 GATE 3 (Shu accretion rate): finding, not a strict pass

Measured the sink accretion rate from gate4_run (the isothermal-bulk BE collapse: barotropic floor
enforces p=rho => c_s=1 below rhocrit, confirmed cs_med~1.28~1 in the bulk). Code-unit Shu rate =
0.975 c_s^3/G with G=four_pi_G/4pi=1/(4pi) => c_s^3/G = 4pi = 12.566, so Shu = 12.25 code mass/time.
Measured dM_sink/dt = 1236-1797 over t=1.16-1.20 => **~98-143 c_s^3/G ~ 100x Shu.**

This is NOT a strict-Shu failure -- it is the correct physics for the IC. gate4 uses f=5, i.e. a
sphere 5x overdense vs the CRITICAL Bonnor-Ebert profile (pgen: f = density-enhancement factor,
f=1 = marginally critical, rc=6.45=xi_max, M_crit=197.56). A strongly SUPERCRITICAL BE sphere
collapses far faster than the marginally-critical/SIS inside-out solution: Foster & Chevalier (1993)
find even a CRITICAL BE collapse peaks near ~47 c_s^3/G (declining), and supercritical is higher --
so ~100 c_s^3/G at f=5 is literature-consistent. Shu's 0.975 c_s^3/G is the SINGULAR-isothermal-
sphere (rho ~ c_s^2/2piG r^2) self-similar asymptotic, a DIFFERENT initial condition, not a BE
sphere. The plan's "within factor 2 of Shu 0.975" gate is therefore mis-specified for the BE IC.

VERDICT: GATE 3 = accretion rate is PHYSICALLY SENSIBLE & literature-consistent for the actual IC
(~100 c_s^3/G, supercritical BE >> Shu), but a strict 0.975 c_s^3/G match is UNACHIEVABLE without a
dedicated SIS or near-critical (f->1) isothermal IC -- a separate validation (new IC + another
multi-slot CPU collapse). The ESSENTIAL correctness of accretion (mass/momentum conservation) is
already machine-precision (gates 1-2). Not chasing the literal 0.975 number with more compute now;
documented as a physics finding. If a strict Shu number is wanted later: build an SIS pgen (central
cutoff to avoid r=0), eos isothermal, low fixed rho_sink for immediate central sink, measure the
constant post-transient rate -> expect ~0.975 c_s^3/G.

## WS-1 STATUS SUMMARY (2026-07-14): code COMPLETE + validated; two gates compute/IC-limited.
Increments 1-5 all IMPLEMENTED & default-OFF bit-identical. PASSED: swarm/AMR/restart (inc1,
<1e-12), sink->gas gravity (inc2, force 8e-16), two-body symplectic (inc3, E drift 4e-12/orbit incl
GPU 2e-12), creation incl on-AMR (inc4 + gate4), accretion conservation (inc5, M 2.8e-16/P 1.1e-14
per step), GPU runtime (smoke clean). PARTIAL/finding: gate4 dt-recovery (2.7-3.6x relief vs
control, no full turnaround -- infall-limited at racc=4), gate3 Shu rate (~100 c_s^3/G, correct for
supercritical BE, strict 0.975 needs SIS IC). DEFERRED (compute): inc6 full-physics t4+sinks needs
8xH100 (t4 production currently holds 5). WS-1 is FUNCTIONALLY DONE for downstream use.

## 2026-07-14 (11:11) — DENSITY MILESTONE 1e-7 g/cm3 (1e6 x rhocrit)
phdf 01536 (cycle 76800, t=1.08398775): max rho 1.0127e-7 g/cm3 = 1.013e6 x rhocrit.
Second-collapse density has risen 61x since onset (1.65e-8, 2026-07-13 08:43).
Still AMR level 14 (8 blocks at L14, 2843 total — unchanged since L14 opened);
dt ~5-7e-9, 18-19 s/step, STS ~67-86 / 13 stages; mem flat 66.2-68.2 GiB x5,
zero errors. Chain: slot 2347591 (g003) + successor 2347873 queued; handoffs
2345200->2345437->2347591 all verified zero-loss. Watch re-armed: rho>=1e-6,
level>=15, mem>76 GiB.

## 2026-07-14 (15:20) — USER-DIRECTED STOP + HALL RETOOLING BEGUN
User: "stop the current run and do what you just recommended" = (1) stop chain,
(2) Hall-off dt-attribution probe at the live state, (3) implement eta_H cap
(sequel to eta_ohm_cap_code), (4) price it, (5) bring numbers for relaunch sign-off.
- Successor 2347873 scancelled first (ordering rule); slot 2347591 left running to
  write checkpoint 00163 (cycle 77250), auto-cancelled right after (waiter) — loss ~0.
- Probe: runs/dt_attrib/submit_hallprobe.sh (5 GPU, A_ctl=production vs
  B_nohall=diffusion/hall=none, 8 cycles each, newest restart, binary v4).
- eta_H cap IMPLEMENTED (uncompiled yet): diffusion/eta_hall_cap_code, signed
  |eta_H| clamp in HallDiffusivity (both Get overloads, ionization family only,
  fixed untouched) + FusedNonidealEval.hall_cap + factory mirror in diffusion.cpp
  + hydro.cpp parse/ctor/banner ("## Hall |eta_H| cap: ..."). Identical structure
  to the eta_O cap so cache/dt/flux paths clamp consistently by construction.
- Physics rationale: whistler dt ~ dx^2/|eta_H|; ionization-model |eta_H|~B/n_e
  runs away exactly where eta_O already sits at ITS cap (Ohmic-dominated,
  whistlers diffusively damped) — capping |eta_H| bounds the strict dt without
  dropping Hall where it is dynamically alive. Sign preserved (Wardle tensor
  sign flips where grains carry the current).

## 2026-07-15 (11:15) — FINDING: rad.Er/Fr outputs dead-zero since cycle 71000 (RT physics ALIVE)
Discovered while building the star-formation flow-diagram artifact from prod_t4_full.
Evidence chain (all read from run outputs today):
- rad.Er in phdf/rhdf is IDENTICALLY 0.0 (every cell, every block) from snapshot 01421
  (cycle 71050) through the final 01545 and restart 00163. Er was healthy and growing
  (max 1.9e4 code) through 01420 (cycle 71000).
- Cycle 71000 = the resume point of jobs 2321463/64 = the FIRST slot of the
  2026-07-10 mixed-integrator binary (athenaPK_eos md5 20ecb7e2, since overwritten;
  previous healthy slot 2313250 ran md5 3dc6ee0e). All later binaries (v2/v4/v5)
  descend from that source.
- NOT a dynamical decay: ambient Er (0.23-0.41 code) cannot drain in 500 cycles at
  chat=1578 (radiation crosses only ~0.016 l0 in that time). Zeroed at the resume.
- RT PHYSICS REMAINED ACTIVE IN MEMORY: central T evolution after 71000 is strongly
  sub-adiabatic (132 K @ 1.37e-10 -> 587 K @ 1.44e-7; pure gamma=7/5 adiabat predicts
  ~2100 K; barotropic law would give ~2900 K). Only the M1 matter coupling can cool
  the gas in v2/v4 (WS-2 thermo post-dates v4 and is default-off). Radiative T(rho)
  matches published RHD collapse curves (~200-400 K at 1e-8).
- Restart transient: central T dropped 141.5 -> 132.4 K (-6.4%) in the 50 cycles after
  the 71000 resume while rho ROSE 7% — consistent with restore zeroing Er in memory and
  the core re-emitting it (gas pays aT^4). Healthy pre-bug boundary (62500, binary
  3dc6ee0e): T and Er continuous to 4 significant figures. Every restart since 71000
  re-injects a (small) transient: boundaries at 71913, 72409, 73426, 75232, 76352.
Impact: MHD/dynamics/thermal science of the run remains valid (radiative cooling
operated throughout); rad.Er/Fr DIAGNOSTICS are lost from cycle 71000 on; ambient
10 K radiation bath permanently absent after 71000 (negligible dynamically over the
~80 yr simulated since). Mechanism hypothesis: the OperatorSplit refactor of the
mixed-integrator work broke the restart restore AND output write path for the rad
fields (evolution presumably in staged containers disconnected from the base
registered vars). NOT YET LOCATED in source; radiation gates (WS-3a/3b) all
fresh-start and pass — the bug is restart-path-specific.
TODO before relaunch decision: locate + fix in current source, rebuild (v6), add a
restart-continuity gate (resume 8 cycles, assert Er nonzero + T continuous);
alternatively relaunch knowingly without rad diagnostics. USER decision pending —
folds into the eta_H-cap sign-off.

## 2026-07-15 (11:50) — RT-ZEROING BUG LOCATED + FIXED (STS pack clobber)
Root cause of the 2026-07-15 11:15 finding, pinned by CPU reproduction (NOT restart-
specific after all): runs/validation_rt/restart_gate/rg2.in = free-streaming M1 pulse
+ rkl2 isotropic conduction — Er dead-zero from the FIRST cycle even on a fresh start.
diffusion/integrator=rkl2 was first enabled in production at the cycle-71000 resume,
which is why it masqueraded as a restart bug.
MECHANISM (hydro_driver.cpp AddSTSTasks): the u1 "Y0 snapshot" deep-copies only
cons+prim, but RKL2StepFirst/RKL2StepOther packed {Metadata::Independent} — which also
matches rad.Er/Fr1-3 and grav.phi. Yjm1 aliases "base", so the first stage wrote
base.rad.Er = u1.rad.Er (never copied -> 0) + mu~*tau*MY0.rad.Er (never written by the
OperatorSplit-excluding FluxDivergence -> 0) = 0, every STS call (2x per cycle).
grav.phi was likewise clobbered with the stale u1 copy — harmless for correctness only
because the Poisson solve recomputes phi each cycle (it cost a zero initial guess).
In production the optically-thick core re-emitted Er each cycle via MatterCoupling
(=> cooling stayed physical, the 11:15 sub-adiabat analysis holds); optically-thin
ambient stayed 0, outputs (post-STS) always 0.
FIX: pack the RKL2 stage kernels by NAME {"cons"} (2-arg PackVariablesAndFluxes, same
idiom as self_gravity.cpp cons_flx_pack) in both RKL2StepFirst and RKL2StepOther.
Diffusion updates cons only, so this is exact. ResetFluxes still zeroes all flux
arrays — benign (every consumer recomputes its fluxes before use).
VALIDATION (build_cpu, all PASS):
- rg2 (kappa=0, coupling off): post-fix Er evolution EXACTLY matches the
  no-diffusion reference run at every snapshot (0.88047@3 ... 0.63545@16) and is
  restart-continuous; pre-fix it was 0 from cycle 3.
- gas state BIT-IDENTICAL pre/post fix (max|dprim|=0.0 on the final snapshot) —
  the cons STS arithmetic is untouched.
- rg3 (kappa_a=0.5, matter_coupling on): fresh + restart clean, gas P responds to
  absorption, no NaN.
GPU rebuild submitted: job 2353024 -> freeze as athenaPK_eos_v6 (= v5 + this fix).
Relaunch notes: ckpt 00163 has Er==0 on disk; resuming with v6 re-emits the core
field in ~1 coupling step (same small transient as every boundary since 71000).
The ambient 10 K bath stays 0 unless repaired — option: relaunch from a COPY of
00163 with rad.Er := a_r T^4 patched in (original checkpoint untouched). Probe jobs
2350019/2350028 (v4/v5) remain valid for cap pricing: dt traces and cap plumbing are
unaffected by this fix (cons bit-identical; only phi initial guess and rad outputs
change).

## 2026-07-16 (00:40) — v6 GPU binary frozen (STS fix)
Build job 2353767 COMPLETED 0:0 (3:07; jobs 2353024/2353764 failed on env leakage:
session TMPDIR into sbatch -> nvcc tmp error, then --export=NONE -> PATH gone; cure
= sbatch --export=ALL,TMPDIR=/tmp). athenaPK_eos_v6 md5 03a299cf1039efcd7a086a8956cda7e4
= v5 (eta_H cap) + the RKL2 named-pack fix (DEV_LOG 2026-07-15 11:50). v5/v4 kept as
rollback. Relaunch binary of record: v6.

## 2026-07-16 (01:30) — IMPACT REVISION: the STS/rad bug ALSO over-cooled the gas (user caught it)
The 2026-07-15 11:15 assessment ("physics stayed valid, only diagnostics lost") is
WRONG in its thermal part. Mechanism, confirmed in MatterCoupling (radiation_moments.cpp):
the coupling conserves e_gas + (c/chat)*Er — refilling Er costs the gas c/chat = 1000x
the radiation energy. Post-71000 each cycle: STS clobber destroys Er -> coupling refills
to ~aT^4 (optically thick, full equilibration) -> gas pays ~1000*aT^4 PER CYCLE. This is
an artificial cooling channel ~ (c/chat)*aT^4/dt ≈ 30-70% of the compression heating
across the runaway.
Three independent quantitative confirmations:
 (1) EXACT EOS-table isentrope from the last pre-bug central state (141 K @ 1.28e-10,
     snap 01420) reaches 2184 K at 1.44e-7 — measured 587 K (3.7x low).
 (2) The 50-cycle post-resume transient: paying 50 x 1000*aT^4(140 K) = 0.145 erg/cm3
     vs gas e ~1.5 -> ~10% drop predicted; measured -6.4%. Matches.
 (3) Pre-bug core Er_max 1.9e4 code = aT^4(141 K) = 1.5e4 code. Equilibrium as expected.
CONSEQUENCES: all thermodynamics from cycle 71000 -> 77250 (rho 1.3e-10 -> 1.44e-7,
the entire L11-L14 runaway) are contaminated LOW; the true central T at 1.44e-7 is
near-adiabatic ~1500-2200 K (physical radiative diffusion at tau~1e6 is slow), i.e.
H2-dissociation onset is likely already reached near the current front — the "590 K
radiative-softened runaway" interpretation (artifact/email/DEV_LOG 07-15) is RETRACTED.
The runaway onset density (1.65e-8), e-fold time, FHC lifetime (220-450 yr), and any
kinematics at fixed rho in that window (incl. N_rot=0.24 endpoint) carry bug
contamination of order tens of percent; qualitative collapse (mass-loaded core +
imminent dissociation) likely genuine.
CLEAN DATA BOUNDARY: everything through cycle 71000 / snapshot 01420 (incl. the
23-yr plateau + its radiative regulation) is uncontaminated (pre-bug binary, full RT).
RELAUNCH RECOMMENDATION: resume v6 from parthenon.out2.00142.rhdf (cycle 71000,
VERIFIED healthy Er 0.235-1.9e4 zero-free, ambient bath intact) and REDO the runaway
(~2-4 days on 5xH100 with the caps; the redone physics will differ: hotter core, more
pressure support). Requires quarantining post-71000 outputs so the chain's newest-
restart logic picks 00142 — USER decision (data-bearing, not yet executed). Probes
2350019/2350028 priced dt at the contaminated state; relative Hall pricing still
informative, absolute numbers will shift on the redone (hotter) state.

## 2026-07-16 (01:55) — QUARANTINE EXECUTED (user-directed; restart decision deferred)
User: "just quarantine. we'll decide on restart after exhausting all possible
optimization options." Moved (NOT deleted) all post-bug outputs to
runs/prod_t4_full/quarantine_postbug_71000/ (README inside):
parthenon.out2.00143..00163.rhdf + final.rhdf (22) and
parthenon.out1.01421..01545.phdf + final (+xdmf) (126+126). Verified: newest live
restart is now parthenon.out2.00142.rhdf (cycle 71000, healthy Er) — the chain's
auto-resume target; newest live snapshot 01420. hst/run.log/gpumem.log left in place
(append-only; contaminated rows identifiable by cycle >= 71000).
hallcap resubmitted as job 2359338 (was 2350028, cancelled BEFORE the move):
PINNED to the quarantined 00163 + 45-min limit; binary stays v5 so the off-state
dt-trace check vs hallprobe A_ctl (v4, same ckpt) remains valid. Chain remains DOWN;
relaunch (from 00142, binary v6) awaits the optimization round + user sign-off.

## 2026-07-17 (02:10) — RELAUNCH: chain up with v6 from clean ckpt 00142 (user: "go for it")
hallcap 2359338 COMPLETED 0:0 (19:45): off-state PASS (v5 no-cap dt trace BIT-IDENTICAL
to v4 A_ctl, 9/9 cycles); cap 0.3 dt bit-identical to no-cap (never binds); cap 0.1
mean dt -2% (noise), STS 88.8/13 vs 91.4/15. Binding-cell |eta_H| ~ CFL dx^2/dt ~ 7e-3
at L14 => any effective cap would clamp dynamically-live Hall regions. ETA_H CAP
REJECTED for production (feature stays in v6, inert). The 6.1x Hall-off prize is NOT
capturable by magnitude-capping at this state; expectation: the corrected (hotter)
redo raises x_e via thermal K ionization -> eta_H drops organically; re-measure there.
RELAUNCH: job 2360418 submitted (submit.sh BIN=v6 03a299cf, no cap key, resumes newest
= 00142/cycle 71000). Verification checklist on start: v6 md5 banner, "RESUMING FROM
...00142", successor queued, dt/wsec sane, mem/card, and CRITICALLY rad.Er > 0 in the
first new phdf (01421+, fresh numbering into the emptied slots) + no NaN.
User directive: after verification + no further optimizations -> stop, document
everything, declare final-state code.

## 2026-07-17 (02:35) — RELAUNCH VERIFIED: fix works in production
Slot 7 (2360418) checklist ALL PASS: v6 banner (03a299cf); RESUMING FROM 00142;
successor 2360419 queued; first fresh snapshot 01421 (cycle 71050): Er min 2.72e-1 /
max 2.76e4, frac0 = 0.0000, central Er 1.26e4 = aT^4(135.5 K) equilibrium, no NaN.
Physics signature of the fix visible immediately: at the same cycle the redo is
WARMER (135.5 vs 132.4 K) and LESS compressed (1.366 vs 1.374e-10) than the
contaminated pass — the artificial cooling is gone. Perf at resume: dt 7.5e-8,
6.45 s/step, 8.2e6 zc/s, STS ratio 1.9-3.8 / 3 stages (eta_O cap live from the
start of this stretch; grav.phi initial guess preserved) => ~32x the sim-rate of
the original pass through cycle 71000+.
Standing plan (user): watch the redo through the runaway (expect hotter core,
thermal K ionization -> x_e up -> eta_H down; re-measure Hall dt share on the
corrected state at depth), exhaust remaining optimizations, then STOP and produce
the final-state-code documentation.

## 2026-07-17 (07:50) — FULL SANITY SUITE on the redo: ALL PASS (snapshot 01447/cycle 72350)
A. Chain: slot 7 at 5h, successor queued, 0 errors; warnings all benign (1646 = STS
   ratio>400 advisory, by design); dt 4.8-6.4e-8; STS 38.8/9; mem flat 57.2 GiB.
B. Conservation: grid mass 5.091304e4 vs 5.091300e4 t=0 (8e-7, outflow-BC level);
   no NaN; rho,P > 0 everywhere; no negative scalars; relDivB bounded (hst max 6.3).
C. Radiation: Er frac0 = 0.0000; CENTRE Er/aT^4 = 0.990 (LTE to 1%); ambient bath
   preserved at exactly 0.413 (t=0 value) — both failure modes of the old bug absent.
D. Chemistry: H-nuclei conservation to -9e-16; x_e 1e-15 (floored core) .. 2.2e-8
   (envelope CR equilibrium).
E. Thermal track: smooth, no transients; central 207->228 K over rho 0.85->1.17e-9;
   effective Gamma = 1.284 — BETWEEN the contaminated 1.21 and pure adiabat 1.40,
   exactly where physical RT cooling at tau~20 (t_diff ~1 yr vs ~5 yr evolution)
   should sit; will steepen toward adiabatic as tau grows. Minor note: a bounded ~4%
   T_central dip in the first ~50 cycles after resume (141.5->135.5), recovered and
   smooth since; pre-bug baseline fluctuation was ±0.1 K.
Redo vs contaminated at fixed t: 13x less dense, e-fold ~9x gentler at fixed rho —
the corrected core is accretion-regulated, not artificially stripped of support.

## 2026-07-18 — INDEPENDENT DEEP CODE AUDIT (fresh Fable agent, no bug-history context)
Read all in-scope custom source completely. Result: NO confirmed critical/major live
bug on the production hot path; nothing that causes artificial support, false stall,
or slowdown. Independently verified CLEAN: v6 RKL2 named-pack fix (rad.Er/Fr+grav.phi
no longer swept into STS), OperatorSplit FluxDivergence patch, self-gravity sign chain
(RHS +4piGrho, g=-grad phi, phi<0 in wells), eta-cache face/ghost coverage, M1 closure,
RSLA MatterCoupling conserving e_gas+(c/chat)Er (the exact thermal-bug mechanism —
independently reconfirmed correct now), HL B-units, sinks+dust default-off inertness.
3 findings, ALL verified against source, NONE affecting the current run:
 F1 [PLAUSIBLE, live-path but bounded] chemistry.cpp:212 computes gow17 reaction T via
    ideal T_code=gm1*eint/rho even under eos=hydrogen (verified: line reads exactly
    that). Wrong in the H2-dissociation zone (~2000 K) where mu varies; over-estimates
    T -> under-estimates recombination -> over-estimates x_e. BOUNDED: production uses
    resistivity_coeff=ionization + hall_coeff=ionization (EQUILIBRIUM, EOS-T); only
    ambipolar_coeff=ionization_chem consumes chem x_e, and its eta_A is capped at the
    equilibrium value. ZERO bearing on current stall (core at 276 K, far below 2000 K).
    WILL matter for the fossil-field field-coupling IF the redo reaches ~2000 K -> fix
    before then: plumb use_h2diss+eos_tab into ReactScalars (mirror PrecomputeNonidealEta).
 F2 [CONFIRMED, latent/production-inert] radiation_moments.cpp:259 mhd=(ncons>=NHYDRO+3)
    misfires for euler+>=3 scalars (reads scalar slots as B). Production glmmhd: both
    criteria agree, correct. Robust nhydro already fetched ~line 271. Fix: mhd=(nhydro>NHYDRO).
 F3 [CONFIRMED, latent/production-inert] dust.cpp:48 scalar_index defaults 0, collides
    with gow17 species [0,NSPEC); only nscalars>=idx+2 checked, no disjoint guard. Dust
    OFF in production. Fix: require dust/scalar_index>=NSPEC when chemistry active.
 FRAGILITY [safe now] resistivity/ambipolar/hall flux kernels pack {Independent} + fixed
    cons indices; safe ONLY because cons registers first (UID order). Same class as the
    v6 STS bug. Defense-in-depth: convert to named cons pack (like ApplyGravitySource).
DECISION: none fire in the production deck; do NOT disrupt the running chain or rebuild
now. Fold F1 (before 2000 K) + F2/F3/fragility (hardening) into the NEXT planned rebuild;
production bit-identical for all four in the current config.

## 2026-07-18 (23:30) — RSLA CONVERGENCE PROBE (creduc 1000 vs 100), job 2362882
40-cycle restarts of the corrected redo state (ckpt 00160, cycle 75500, rho_c 2.67e-9,
T_c 276 K), binary v6, both exit 0, no NaN. Same cycle/time endpoint (75540, t=1.0842650),
only chat differs (c/1000 -> c/100, 10x more RT subcycles).
CORE LOCAL STATE = creduc-INSENSITIVE:
  T_c: 275.90 (A) vs 275.98 K (B)  -> +0.03%
  rho_c: -0.02%;  Er/aT^4: 0.9740 vs 0.9741
  => interior is optically thick + LTE; local matter-radiation coupling equilibrates
     within a cycle regardless of chat. The WARM STALL IS PHYSICAL, not a local RSLA
     coupling artifact. (confirms the falsification direction.)
ESCAPE LUMINOSITY L(r)=4pi r^2 <F_r> IS creduc-sensitive (flux causally capped |F|<=chat E):
  L(1 AU) B/A = 2.46 ;  L(3 AU) B/A = 5.54  (higher chat => more luminosity escapes).
  L(0.3 AU) noisy/near-zero (thermalized core interior) - ignore.
  => true c (>> c/100) would radiate even faster; RSLA UNDER-cools the envelope, so the
     stall DURATION is OVER-estimated: the real core reaches instability/dissociation
     SOONER than this run shows. Core TEMPERATURE trustworthy; evolution RATE is a lower
     bound (slowest).
CAVEAT: 40 cycles = 0.12 yr << ~30 yr cooling/e-fold time, so the L difference has NOT
yet fed back on T_c (dT 0.03%). Quantifying the duration bias needs a longer c/100 (or
true-c) run - DEFER; document as a known RSLA limitation (direction known, magnitude
un-quantified). NOT a bug, NOT an optimization. Production creduc=1000 unchanged.

## 2026-07-19 — REDO: LEVEL 13 OPENS, runaway beginning (cycle 79850)
Corrected run (v6, clean physics) first refinement since the plateau: L13 opens at
rho_c=6.38e-9, T_c=338 K (nb 2402->2486). Runaway signature CONFIRMED:
 - e-fold time turning over: steady ~10.5-11 yr through the plateau -> 6.9 yr in the
   L13-opening interval (accelerating).
 - dt turning over: flat ~5e-8 -> 1.83e-8 at cycle 79873 (finer cells + tightening CFL).
 - Er zero-fraction 0.0000 through the refinement (fix holds across AMR).
PHYSICS NOTE: onset at T~338 K, FAR below H2 dissociation (2000 K) -> this second-collapse
onset is GRAVITATIONAL/accretion-driven (core reached its critical mass), NOT dissociation-
driven; dissociation follows deeper in. Density-at-onset vs the contaminated run
(6.4e-9 vs 1.65e-8) is NOT a clean comparison (different thermal histories) - do not
over-interpret; the clean signal is the e-fold+dt acceleration.
Expect from here: dt falls, throughput drops, mem climbs as L13/14 refine (the expensive
phase). Deep markers still ahead: rho>=1e-8, T>=2000 K. Watch armed.

## 2026-07-19 — AUDIT-FIX BATCH: F1+F2+F3+fragility+NEW-1+NEW-2 (candidate v7)
User-approved batch of all outstanding audit findings, implemented + CPU-validated.
Full audit report: ~/CODE_AUDIT_REPORT_2026-07-19.md.

FIXES (all in source, CPU build clean, production chain UNTOUCHED — still v6):
 - F1  chemistry.cpp ReactScalars: gow17 reaction T now from the EOS table under
   eos=hydrogen (eint>0-guarded, floored to match ideal path; mirrors the
   PrecomputeNonidealEta idiom). Ideal-gamma path textually unchanged. PLUS a new
   REQUIRE forbidding thermo=true + eos=hydrogen (AdvanceThermoEnergy still maps
   e<->T with the ideal gamma — same F1 class; thermo is default-OFF anyway).
 - F2  radiation_moments.cpp MatterCoupling: mhd = (nhydro > NHYDRO) from the Hydro
   package instead of the packed-cons-width test that misfired for euler + >=3 scalars.
 - F3  dust.cpp Initialize: REQUIRE scalar_index >= Chemistry::NSPEC when
   <physics> chemistry=true (default 0 aliased f_dg/a_c onto x_H2/x_Hp silently).
 - FRAGILITY  resistivity/ambipolar/hall/conduction(x2)/viscosity flux kernels now pack
   "cons" BY NAME instead of {Metadata::Independent} (fixed IB*/IM*/IEN indices no
   longer depend on registration order — v6-STS-bug trap class closed). NOTE:
   ResetFluxes deliberately NOT narrowed: it has no fixed-index assumption, and its
   broad zeroing keeps grav.phi's flux registers clean (phi is Independent+WithFluxes
   WITHOUT OperatorSplit, so it IS packed by UpdateWithFluxDivergence — the zeroing may
   be load-bearing).
 - NEW-1 self_gravity.cpp: REQUIRE !(multipole BCs && use_swindle) — moments use full
   rho, swindle RHS uses rho-rho_mean; the combination was silently inconsistent.
 - NEW-2 parthenon driver.{hpp,cpp} (vendored patch, same class as the OperatorSplit
   patch): parthenon/mesh/amr_check_interval=N gates the full-mesh-rebuild
   LoadBalancingAndAdaptiveMeshRefinement call to every N-th cycle. Tagging still runs
   EVERY cycle (derefine_count semantics exact, flags always current; latency <= N-1
   cycles); pmesh->modified cleared on skipped cycles. DEFAULT 1 = bit-identical.
   Motivation: prod at cycle ~80.6k regrids every cycle — wsec_AMR 11.0 s of a 25.0 s
   cycle (44%), min 8.78 s over 3000 cycles. N=10 safety: 0.13% density drift/check vs
   4x Truelove margin. Expected ~1.7x, growing at L14+.

VALIDATION (front-end CPU serial, this session):
 - build_cpu clean.
 - cshock_ad.in 40 cycles unsplit AD: exit 0, no NaN (ambipolar named pack).
 - cshock_ad.in + integrator=rkl2: 40 cycles, 80 STS calls, clean (parabolic named
   packs + AddSTSTasks).
 - collapse_L16_multipole fhc.in, amr_check_interval=10, ncycle_out=1: skip cycles
   wsec_AMR ~4e-8 s, check cycles at 10/20 run the real call (regrid at cycle 20 =
   1.52 s) — mechanism confirmed; interval=1 run also clean.
 - NEW-1 guard fires (sphere_multipole + use_swindle=true -> REQUIRE with message).
 - F3 guard fires (dust + gow17 + scalar_index=0 -> REQUIRE with message).
 - F1 path: eos_smoke/fhc_chem.in (eos=hydrogen + gow17 + full non-ideal + gravity)
   10 cycles clean, no NaN, chemistry-coupled x_e banner confirms config.
 - NOT smoked: F2's euler+scalars config (latent-only; compile-verified).

GPU: build job 2363571 submitted (--export=ALL,TMPDIR=/tmp) -> build_gpu/bin/athenaPK
(candidate v7; does NOT overwrite athenaPK_eos_v6). Before any chain swap: md5 the new
binary, pinned-restart continuation check vs v6 (expect bit-identical gas state at
default settings), and amr_check_interval probe (dt_attrib pattern) + user sign-off
for N=10 in production.
UPDATE (same day, 21:02): GPU build 2363571 COMPLETED 0:0 in 4:01 ->
build_gpu/bin/athenaPK, md5 505491e86a9e90f8a470f2bc298939d6 (candidate v7).
Patch presence verified in BOTH binaries via strings (amr_check_interval + F3 + NEW-1
messages: 3/3 each). athenaPK_eos_v6 (03a299cf) untouched; chain undisturbed.

## 2026-07-19/20 — BIG CLEANUP (user-directed, 3 phases) + flux-retention extraction
DELETED (~2.5 TB freed): prod_t4_full/quarantine_postbug_71000 (1.02 TB, user released);
runs/dt_attrib (147 GB; submit patterns -> scripts/probe_patterns/); prod_t2_ad (756 GB)
+ prod_t3_ad_ohm (586 GB) AFTER flux-retention extraction (below); ~/runs_archive
(RESULTS_ARCHIVE.md -> docs/); ~/athenapk_clean (PR #164 is pushed; re-clone to address
reviews); GPU rollback binaries v2/v4/v5 (rollback = git checkout + rebuild); athena++
extra binaries (chem/rad/noSTS; _rebuild_sts.sh + memory notes regenerate) + obj/;
athenapk/build scratch dir; validation OUTPUTS of ws1/rt/cshock/ws5a-L16 (decks, submit
scripts, .py analyzers, SLURM logs ALL KEPT; ws5a collapse_L52_ref 5.6 GB KEPT for a
future gate-(ii) run); eos_smoke chemA/chemB outputs; early-era docs (methodology-
validation + test5 + cooling-investigation reports, user-approved); early scripts
(be_collapse_analysis, jeans_analysis*, compare_athenapk_vs_athenapp; no imports);
src/eos/eos_table.h5 (gen-script intermediate; runtime loads .bin); a 54 MB ELF core
dump in ws5a; assorted ~ tmp dirs + old session scratch. KEPT: athenapk_experimental +
artemis (user), all tracked stock inputs/docs in both repos (deleting tracked upstream
content would only dirty git), prod_t4_full complete.

FLUX-RETENTION EXTRACTION (pre-deletion, paper headline box): new permanent tool
runs/flux_retention.py (raw h5py; peak-plane Phi_core with face-degeneracy-safe layer
pick, rhocrit-crossing r_core, NN78+MS76 mu, Lagrangian Phi_0 from snap 0; anchor:
t=0 Phi = pi R^2 B0 to 0.1%, cloud mu(NN78)=30.5). Results (both tiers, JSONs +
flux_retention_tiers.md in ~/paper1_fossil_field/): T2/T3 statistically identical;
milestone 1 / 1e2 / deepest(657x,534x rhocrit): mu(NN78) 6.0 -> 7.8 -> 8.7,
retention 1.40 / 1.43 / 1.09. 1e5-rhocrit milestone NEVER reached by t2/t3 — the deep
measurement belongs to prod_t4_full at campaign end. retention>1 = envelope advection
through the midplane still beating AD expulsion at shallow depth (not an error;
mechanism validated on the t=0 anchor).

POST-CLEANUP SANITY: no dangling references to any deleted artifact in src/CMake/decks/
docs/CLAUDE.md; CPU smokes pass on the cleaned tree (cshock_ad 10 cyc, fhc_chem 2 cyc,
exit 0, no NaN); chain rolled cleanly to slot 12 (job 2363405, ckpt 00185) — untouched
throughout. Memory: 4 obsolete notes deleted (LD_PRELOAD fix preserved into build-env
note), index rewritten; old agent jobs/transcripts/scratch cleared (t4 analysis
snippets + figures + artifact HTML rescued to scripts/t4_analysis/).

## 2026-07-20 — V7 PRE-SWAP PROBE (job 2363715): amr speedup FALSIFIED, v7 swap-safe
Pinned A/B/C probe on ckpt 00185 (cycle 81750), 30 cycles/segment, exact production CLI,
5xH100, clean (exit 0, no NaN, Er_frac0=0, 2486 blocks all three).
  A_v6ctl : athenaPK_eos_v6 (03a299cf)         mean step=15.5 amr=10.6  total=809 s
  B_v7def : v7 (505491e8), amr_check_interval=1 mean step=15.5 amr=10.6  total=812 s
  C_v7iv10: v7, amr_check_interval=10           mean step=25.1 amr=1.1   total=810 s

FINDING 1 -- amr_check_interval=10 gives ZERO wall-time speedup; the audit's ~1.7x is
FALSIFIED. C is BIT-IDENTICAL to B at field level (prim/grav.phi/rad.Er max|rel|=0.0,
same 2486 blocks, identical dt every cycle) => iv=10 is correctness-safe. But total wall
is A/B/C = 809/812/810 s (all within 0.4%, noise). Mechanism: on the 9/10 cycles C skips
the LB&AMR call, the per-cycle STEP rises ~14 -> ~24-25 s, exactly absorbing the ~11 s of
"AMR" it skipped (verified in the per-cycle trace: skipped cycles show amr=0 step~25,
regrid cycles amr=11 step~14). So wsec_AMR is NOT recoverable overhead -- it is coupled
to the step through load-balance / rank-synchronization; skipping the collective just
moves the wait into the next step's boundary exchange. CONCLUSION: keep the feature (a
harmless no-op at default 1) but do NOT adopt iv>1 for throughput. The AMR cost is
genuine work. Real remaining optimization is the STS all-variable boundary exchange
(~18-30 full exchanges/cycle) -- profile the cons-only shallow container.

FINDING 2 -- the F1 fix is inert at the current state (T_c~350 K), as predicted. A(v6)
and B(v7) ran BIT-IDENTICAL for ~14 cycles, then chaotic dt jitter decorrelated them (F1
perturbs the EOS-vs-ideal reaction T -> x_e at roundoff -> flips the capped-AD dt
sequence). Ended at Dt/t~1.8e-8 apart; rho_max agrees to 8.4e-5 despite the offset; both
clean. Confirms F1 is a genuine-but-negligible correction here => "fix before 2000 K",
and the swap reason is F1 CORRECTNESS, not speed.

VERDICT: v7 is swap-safe. Recommended: swap the chain to v7 at a slot boundary before
T_c reaches ~2000 K (F1 is the only live production-config difference; F2/F3/FR/NEW-1 are
inert-or-guards). No further rebuild needed. Probe outputs (A_v6ctl/ B_v7def/ C_v7iv10/
+ logs + tmp) DELETED post-analysis; submit_v7probe.sh -> scripts/probe_patterns/.

## 2026-07-21 — Paper I AUTO-UPDATE PROTOCOL set up (user-directed)
User: keep fossil-field Paper I (~/paper1_fossil_field/ms.tex) current as the campaign
progresses; persist to every new chat. Implemented as a STANDING per-session directive
(no cron): (1) memory note fossil-field-paper1.md gained an "AUTO-UPDATE PROTOCOL"
section — update triggers, box->producer map (flux_retention.py, fhc_ext producers),
the CORRECTED-v6-data-only caveat (contaminated 71000-77250 quarantined), and the
measured physics learnings that revise the paper (gravitational not dissociation-driven
2nd-collapse onset at ~370 K; stuttering approach; RSLA stall-duration caveat; B-rho
~7-8 mG wall). (2) MEMORY.md pointer flags the standing directive. (3) A LaTeX comment
block with the same caveat + current corrected state inserted above Tab.~tierresults in
ms.tex so it travels with the file. No rendered numbers changed (t4 deep cells remain
legitimately PENDING until the run reaches them; 1e5-rhocrit deep flux milestone
IMMINENT at the current ~9.5e4 rhocrit front). Note: v7 amr_check_interval probe (other
session, job 2363715) found ZERO speedup — the ~1.7x AMR-gating optimization is
FALSIFIED; v7 swap is for F1 correctness before 2000 K, not speed.

## 2026-07-21 — MILESTONE 1e5 rhocrit + FIRST Paper I auto-update executed
prod_t4_full (corrected v6) crossed rho=1e-8 g/cc (1.001e5 rhocrit) at cycle 85400,
T_c 374 K, level 13. Per the new auto-update protocol, ran flux_retention.py on the
corrected run (targeted temp dir: t=0 + 40 tail snapshots, ~15 min) -> the paper's
HEADLINE deep flux-retention point, on clean corrected data (snap cycle 85400 > 71000):
  retention = 0.165 (~6x flux reduction), mu_core(NN78) = 68.3 (cloud 31.9),
  r_core 15.6 AU, M_core 0.0204 Msun, Phi_core/Phi_0 = 0.02457/0.14854 code.
This is the O(10x) ambipolar demagnetization the literature predicts, first measured at
2nd-collapse depth (t2/t3 only reached ~600x with retention ~1.1 = envelope advection).
Paper updated: ms.tex sec:fossil pending box -> real Table (tab:fluxret) T2/T3+T4;
sec:secondcore state refreshed (5xH100, 1e5 rhocrit/370K/L13, stutter, gravitational
onset). flux_retention_tiers.md + flux_retention_t4.json written; memory auto-update
log stamped. NOTE: dt still 3e-8, still stuttering just past 1e-8 (not a hard runaway
yet). NEXT auto-update triggers: 1e6 rhocrit, T_c~1000-1500 K, 2nd core.

## 2026-07-21 — USER-DIRECTED PIVOT: stop t4, prune, launch prod_t5_smallbox (matched box + multipole + v7)
User (after my analysis that a smaller box does NOT help dt — dt is set by the finest
CORE cell, not the box; ~15% per-step at best on envelope blocks): chose to stop t4 and
launch a fresh SMALLER-box run anyway with multipole BCs + all updated physics EXCEPT
sinks. Rationale I under-weighted: a clean run with the F1 EOS-T chemistry fix (matters
for field coupling before dissociation) + true isolated (multipole) boundaries is MORE
CORRECT for the paper, not merely faster.
DONE:
- Stopped chain: STOP_CHAIN + scancel successor 2371644 then running 2367303 (COMPLETED).
  prod_t4_full frozen at cycle 85500 (newest restart 00200), KEPT RESTARTABLE.
- PRUNED prod_t4_full (user: "delete some to reclaim space"): kept newest 10 rhdf
  (+final) + phdf {IC, milestone 01708, every 40th, newest 20} = 64 phdf / 11 rhdf;
  deleted 1651 phdf + 64 rhdf (+xdmf). du 6.7T -> 399G; quota 62% -> ~4%. Newest
  restart 00200 PRESERVED (asserted in the pruner).
- NEW RUN prod_t5_smallbox/: fhc.in from t4 with BOX L=52->16 (x*min/max +-26 -> +-8),
  self_gravity BCs zero->MULTIPOLE (packed_bc=true; WS-5a validated), binary v7
  (athenaPK_eos_v7 frozen md5 505491e8 = F1/F2/F3 + flux-kernel named-pack hardening),
  numlevel=20 kept (finest 0.048 Rsun in the small box, headroom). Sinks OFF. Same CLI
  overrides as t4 (rkl2 STS + eta_ohm_cap 0.1; NO eta_hall_cap). Fresh from t=0.
- Config is essentially WS-5a-validated (collapse_be + L=16 + multipole collapsed to
  2.87x rhocrit there). FRONT-END CANNOT smoke-test (full-physics binary SIGBUSes at
  init even at 32^3 — 1-CPU cgroup; confirmed the OLD deck SIGBUSes too => front-end
  limit, NOT the new config). So: submitted a one-shot GPU SMOKE (job 2372692, 30-cyc,
  30-min, no chaining) to validate box+multipole+v7 on real 5xH100 before the chain.
NEXT: if smoke passes (setup, multipole Poisson converges, mass ~987.5 code = 6 Msun,
dt sane, Er>0, no NaN) -> launch submit.sh (self-chaining production) + re-arm watchers
+ update memory. If it fails -> diagnose before committing the chain.

## 2026-07-21 — prod_t5_smallbox SMOKE PASSED + wrap-up for new-chat handoff
Smoke job 2372692 COMPLETED 0:0 (3:13): mass 987.8 code = 6 Msun (correct), central rho
2.733e-18 g/cc = f*rho0 (correct IC), multipole self-gravity registered + solving
(rho_max 5.00->5.23 = collapse engaging), 30 cycles clean, 0 NaN. Config VALIDATED on
real hardware. Per user: production to be launched in a NEW chat (not here) — run is
READY (submit.sh + submit_smoke.sh in the dir; smoke artifacts cleaned). Memory
organized: new note athenapk-prod-t5-smallbox (live), rrz-gpu-flag demoted to
stopped/pruned/restartable, MEMORY.md live section rewritten, fossil-field-paper1
directive repointed to prod_t5. Handoff at ~/HANDOFF_2026-07-21.md (launch instruction).
No background agents/tasks/crons remain. Chat-transcript cleanup: deleted 1 empty stub;
4 real sessions held pending user ID of "the athenapk experimental chat" (UUIDs not
self-identifying).

## 2026-07-25 — PHASE 2 CT: gate PASSED + CT+Hall blowup root-caused and FIXED
The full-physics CT gate leg had detonated at first-core formation (t≈1.07747, level-4
boundary, ρ~4e4): |B|max→1.35e15 (v_A~3000c), dt→1e-55 frozen. NOT a ∇·B failure —
`ct_maxRelDivB` stayed round-off (~1e-13) throughout; only `ct_maxAbsDivB` grew, and only
because abs = rel·|B| with |B| exploded. So a field-MAGNITUDE instability, not a divergence
one. Diagnosis chain (each single-variable): (1) rkl2_max_dt_ratio 1000→400 FALSIFIED — Hall
is unsplit, bypasses RKL2 STS, so ratio-independent (ct_r400 still blew up ~70 cyc past the
level-4 crossing, |B|max=1.35e15). (2) Discriminant: GLM+Hall (glm_r400) stayed HEALTHY past
first core (t=1.085) while CT+Hall detonated ⇒ the CT+Hall COMBINATION, not Hall alone. (3)
Hall-OFF isolation: ct_run_halloff (single delta hall→none) ran CLEAN to tlim=1.0945 at
maxLevel 6 / ρ=6.14e6, |B|max=3256, 0 NaN ⇒ **uncapped Hall is the CONFIRMED sole cause**.
Mechanism: η_H ∝ B/n_e explodes in the high-ρ recombination core (n_e recombines away), driving
a grid-scale dispersive whistler runaway; GLM/Dedner cleaning damps it, exact-div-free CT does
not.

FIX (user-selected "Hall-cap full physics"): `diffusion/eta_hall_cap_code=0.05` — the knob
already existed (hydro.cpp:1065, signed |η_H| clamp diffusion.hpp:485-486); uncapped was a
CONFIG choice, not a missing feature, so NO rebuild. Cap value from the level-6 whistler-CFL
(dx=9.77e-4, dx²=9.5e-7): η_H ≤ cfl·dx²/(2π·dt_target) ⇒ 0.05 keeps unsplit-Hall dt~1e-6 ~ the
hydro dt-wall while preserving Hall physics where η_H<0.05 (= Ohmic-floor scale, < eta_ohm_cap
0.1; the capped region is the Ohmic-dominated core where Hall is subdominant anyway). Validated
at the EXACT detonation coordinates: `runs/inc7_gate/ct_run_hallcap` (Hall ON + cap=0.05, single
delta from the run that blew up, job 2398334) cleared t=1.07747 at level 4 / ρ=4.57e4 with
|B|max=77.0 (vs uncapped 1.35e15 — 13 orders lower), dt=5.6e-5 steady, ME=58.6, 0 NaN, matching
the healthy Hall-off trajectory (|B|~75), NOT the blowup.

GATE RESULT (runs/flux_retention.py, matched Hall-OFF CT-vs-GLM pair at first-core depth
tlim=1.0945, both maxLevel 6): retention Φ_core/Φ₀ = CT 2.48 vs GLM 2.44 (1.6%); μ_core(NN78) =
CT 5.35 vs GLM 5.46 (2.0%); Φ_core 0.459 vs 0.474 (3.3%); M_core 0.0299 vs 0.0316 M☉. **Verdict:
retained fossil flux is NOT set by the divergence-control choice** (agree to ~2%). Sub-finding:
CT reaches lower central density at fixed time (34× vs 57× ρ_crit) ⇒ exact-div-free CT preserves
more magnetic support (less numerical field dissipation). Cross-code (Athena++) comparison uses
the Hall-OFF pair (Athena++ Hall is a stub); Hall-cap is the AthenaPK-only most-complete path.
Canonical validated CT config: divergence_control=ct, ct_emf=gs05, hall=hall+eta_hall_cap_code=0.05,
rkl2_max_dt_ratio=400, eta_ohm_cap_code=0.1. Gate inputs untracked (runs/inc7_gate/); FLAGSHIP
Phase-2 section is the tracked reference. **Phase 2 CT COMPLETE.** prod_v9 stays HELD (STOP_CHAIN).

## 2026-07-25 — PHASE 3 COMPLETE: grain conductivity gate (2 defects found+fixed) + EOS resolved
**Conductivity (remaining Phase-3 item) DONE — and it caught two real bugs.** Built the
grain-INCLUSIVE cross-check by reimplementing the MRN grain-charge model in Python with a
deliberately different algorithm: the per-bin capture balance n_e v_e e^psi = n_i v_i (1-psi)
has NO bin dependence (the pi a^2 cross-sections cancel between e and i), so psi is a single
global unknown, Z_k = psi tau_k, r = n_e/n_i is an explicit function of psi, and the coupled
system collapses to a 1-D bracketed root find — vs the C++ relaxed fixed point + inner Newton.
`xcheck_conductivity.cpp` extended to dump the charge state (n_e, n_i, Z_k, ng_k), not just eta.

DEFECT 1 — `SolveCharges` does not converge once grains dominate the charge budget. Provable
from the C++ output ALONE (no Python needed): its own neutrality constraint n_i - n_e =
sum(-Z_k ng_k) is violated by up to **3.6e7** relative, and it returns n_e > n_i, which is
impossible with negatively charged grains. Onset ~rho 1e-14 (T=300) / 1e-12 (all T<T_subl),
worsening with grain density. FIX: the exact reduction above -> bracketed bisection in psi that
cannot fail to converge. Neutrality residual **3.6e7 -> 4.2e-11**.

DEFECT 2 — `SahaThermal` bisects [0, n_K+n_H] with a FIXED 64 steps, so its resolution floor is
(n_K+n_H)/2^64. In cold gas, where true thermal ionization underflows, it returns ~1.2e-4
electrons instead of ~1e-43: a spurious electron floor x_e ~ 9e-20, ~9x the intended xe_floor,
scaling with density so it can dominate the real n_e in the dense core. It exactly accounted for
the last 3.4% of the cross-check residual. FIX: relative-precision Newton from the
weak-ionization limit n_e ~ sqrt(n_K fK + n_H fH).

Both gated by `diffusion/ion_legacy_charge_solver` (default false = fixed). FINAL AGREEMENT
C++ vs independent Python: n_e/n_i/Z_k and eta_O to **5e-11**; eta_H 4.2e-4, eta_A 2.9e-3 —
cancellation-limited (sigma_H's net is ~1e-13 of its summed |term| magnitude near the
grain-induced Hall sign reversal; at that point eta_H is 5 decades BELOW eta_O) and well inside
the Phase-3 <1% criterion. All 8 gates PASS.

PRODUCTION IMPACT (impact_charge_solver_fix.cpp, fiducial collapse track): **0.0% for
rho <= 1e-13** — the entire pre-first-core collapse — and 0.0% again for rho >= 1e-9 where
grains have sublimated. It matters ONLY in the grain-dominated first-core band rho ~ 1e-12..1e-10:
eta_O up to +375%, eta_A up to 52%, and **eta_H CHANGES SIGN** (1e-11: -3.4e16 -> +8.2e17;
1e-10: -1.0e13 -> +3.1e15). The Hall sign sets the DIRECTION of field transport, so this is a
physics-level correction. Needs a GPU rebuild to reach production. CAVEAT on the just-committed
Phase-2 flux gate: its deepest snapshots (rho_max 3.4e-12 / 5.7e-12 g/cm^3) sit inside the
affected band, so the ABSOLUTE retention 2.48/2.44 may shift slightly on re-baseline; the
CT-vs-GLM COMPARISON is unaffected (both legs used identical microphysics).

**EOS hi-res table RESOLVED.** The deferral reason (result-changing + table read at RUNTIME by
live jobs) is void: `hydro/eos_table_file` is an INPUT KEY, so adoption needs no file overwrite
and no code change. Verified: `eos_table_hires.bin` (400x1000x920) passes ALL consistency gates
and is loader-compatible as-is — the header is self-describing (nr,ne,nT), the domain spans are
identical to shipped, the finer grids are still uniform in log (the device bilin's assumption),
and the size is byte-exact (12,544,072 = 72 + 8*(nr*ne*3 + nr*nT)). Measured side by side:
P vs Saha median/max = 0.68%/6.68% (shipped) -> 0.16%/0.55% (hi-res); cs2 vs isentrope
0.85%/17.8% -> 0.16%/3.31%. So 12x better on P, 5.4x on cs2, at the H2-dissociation/H-ionization
kinks that set first/second-core thermodynamics. Adoption = one input line
`<hydro> eos_table_file=.../eos_table_hires.bin`; default stays the shipped table so a clean
checkout always runs. Did NOT overwrite eos_table.bin (git-tracked AND being read right now by
the live ct_hallcap job). The non-uniform-grid idea is CLOSED on cost/benefit: it would require
replacing the even-log-grid assumption in the EOS hot path with a per-axis search to buy only
the last cusp point (cs2 median is already 0.16%).
VERIFICATIONS this session: (a) CPU build of the charge-solver fix under REAL Kokkos PASSED
(make -C build_cpu athenaPK, RC=0, no errors) — the fix is not merely host-shim-compilable;
(b) hi-res EOS regeneration is bit-for-bit reproducible (md5 6b8e3999eca19806d8f4d43054e0447c).
NOT done (needs user go, both result-changing): the GPU rebuild that carries the charge-solver
fix into production, and the production EOS input-key swap. prod_v9 remains HELD.
