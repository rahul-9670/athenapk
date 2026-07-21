# RESULTS_ARCHIVE.md — superseded-run results (archived 2026-06-20 before aggressive cleanup)

The heavy output data (phdf/rhdf/athdf/rst/xdmf, ~1 TB) of the runs below was deleted on
2026-06-20 to make room for the 8-run physics-tier matrix (`{athenapk,athena++}/runs/mhd[1-4]_*`).
Their **lightweight records are preserved** in the per-run subdirs here (input files, submit
scripts, `.log`/`.out`, `.hst` history, `units.json`). Headline science below; full per-cycle
detail lives in the preserved `.log`/`.hst`. See also memory `athenapk-vs-athena-comparison`.

## Headline results (GPU AthenaPK vs CPU Athena++, matched FHC collapse, L=16 box, root 256^3, njeans=16, numlevel=16)

- **Per-zone-cycle hardware (ideal pair, run_ideal_match / ideal_match):** 1 H100 ≈ 29–32 CPU
  cores; 2 H100 ≈ one 64-core std node. Wall-clock to same sim-time ~3.5× GPU. The solver is
  self-gravity-multigrid + AMR-regrid latency/comm-bound, so GPUs shed most of their FLOP edge.
- **To first-core density (AD pair, run2match / ad_run2match), the FAIR axis:** GPU reaches
  rho_crit in ~3.1 h / ~447 cycles vs CPU ~5.0 h / ~812 cycles → **~1.6× wall, ~1.8× fewer
  cycles.** (The old "36×" was a matched-*sim-time* artifact — at fixed t the CPU had already
  raced into the dt-collapse while the GPU was still pre-formation.)
- **Formation-time offset (~4%):** AthenaPK forms first core at t≈1.076, Athena++ at t≈1.033.
  NOT AMR granularity (ruled out by run2match_mb16) and NOT an AD bug (AD validated separately);
  attributed to GLM cell-centered B + divB cleaning (AthenaPK) vs CT face-centered B (Athena++)
  and AD-EMF discretization.
- **Meshblock size (2026-06-19/20):** 16³ blocks WRECK GPU throughput (~110 s/cyc vs ~25 s/cyc
  at 32³); match the *physics grid*, not the block decomposition — AthenaPK 32³, Athena++ 16³.
  Athena++@32³ test (ad_run2match_mb32): 2.3× faster than 16³ pre-collapse but over-refines in
  the core (net wash by formation); at equal 32³, AthenaPK stays 1.7–12× faster.

## Runs archived here
- `athenapk__run1`, `athenapk__run2` — original (unmatched L=52 box) ideal & AD production.
- `athenapk__run2match`, `athena++__ad_run2match` — matched-box AD production pair.
- `athenapk__run_ideal_match`, `athena++__ideal_match` — matched ideal pair.
- `athena++__test9` — Athena++ ideal baseline (mu=20, 2-stage).
- `athenapk__run2match_mb16`, `athena++__ad_run2match_mb32` — meshblock-granularity controls.
- `athenapk__run_gpu_nlb_tcp`, `__run_gpu_nlb_smc0`, `__run_gpu_tuned` — GPU transport/launch benchmarks.
