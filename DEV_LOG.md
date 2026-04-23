# Self-Gravity Port Development Log

## Current state (end of Apr 23-24, 2026)

### What's done and validated
- Self-gravity port from Artemis → AthenaPK (hydro + MHD)
- Jeans dispersion tests: hydro 1.245% error, MHD 1.272% error (cluster)
- BE collapse tests: hydro matches laptop-vs-cluster to 3 sig figs;
  MHD shows physical flux-freezing and magnetic support of collapse
- Supported on macOS/Clang (laptop) and Linux/GCC 13.3 + OpenMPI 5.0.7 (Hummel2)
- Parthenon v25.12 (commit d467a382b) — grid.type API fix applied

### Athena++ code comparison (MHD BE collapse)
Done at uniform grid (no AMR yet — AthenaPK doesn't support Jeans criterion).
  - 32^3: rho_max final agrees 97%, B_max 86%
  - 64^3: rho_max final agrees 98%, B_max 95%
  - Peak-collapse offset of ~0.1 time units (expected GLM vs CT diffusion)
  - Both codes under-resolve Jeans length at peak rho (uniform grid limitation)

### Known issues / TODOs
- **Neumann gravity BCs in AthenaPK broken**: gauge not fixed for all-Neumann case
  → dt crashes at cycle 1. Use `zero` (Dirichlet) for now.
- **Multipole gravity BCs not implemented** in AthenaPK (Athena++'s default for test5)
- **Jeans-length refinement not implemented** in AthenaPK
  → required for resolved collapse comparison vs Athena++ production runs

### Next development priorities
1. Implement Jeans refinement in AthenaPK (src/refinement/)
2. Re-run comparison with matched Jeans AMR on both codes
3. Consider multipole gravity BCs (needed for production-scale collapse in small boxes)

### Environment setup
  - ~/athenapk_env.sh loads GCC/OpenMPI/HDF5 paths + cmake
  - /beegfs/u/bbg6470/venvs/analysis_env/ has h5py, numpy, scipy for analysis
  - Comparison data in /beegfs/u/bbg6470/athenapk_runs/comparison_athenapk/
    and /beegfs/u/bbg6470/athena++/test_collapse/runs/comparison/
