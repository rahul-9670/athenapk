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

