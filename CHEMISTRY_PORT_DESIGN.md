# AthenaPK GPU chemistry network — port design

Status: DESIGN (2026-06-23). Companion to the completed M1 RT port and the NICIL-class
non-ideal MHD work. This is the "remaining port" (RT done; sinks dropped; chemistry greenfield).

## What we're porting
Athena++ has a full CPU chemistry framework (`src/chemistry/`): networks `gow17` (Gong+2017
H–C–O, ~13 species), `H2` (minimal, ~3 species), `kida` (KIDA database), solved per-cell with
**CVODE** (serial, stiff BDF) or an explicit `forward_euler` sub-cycler; plus `thermo`,
`shielding`, KIDA parsers. CVODE does not run on GPU, so the AthenaPK port needs a
**device-callable network + a device stiff-ODE integrator**.

## Key enabler (discovered)
AthenaPK ALREADY advects passive scalars: `<hydro> nscalars = N` registers `scalar_density_i`
(conserved) + `scalar_i` (primitive), packed into cons/prim and fluxed with the hydro
(hydro.cpp:923,952,1283). **Species = these scalars** (number density relative to n_H, or mass
fraction). No greenfield scalar/advection infrastructure needed — only the reaction source term.

## Architecture (mirror the M1 RT package)
1. **`src/chemistry/chemistry.{hpp,cpp}`** — a `chemistry` StateDescriptor, created in
   `Hydro::ProcessPackages` when `<physics> chemistry = true` (exactly like `radiation`).
   Reads `<chemistry>` block: `network` (`H2`|`gow17_reduced`), `nsub_max`, `cfl_cool`,
   CR rate `zeta`, unit scales (share the FHC calibration already in `ionization.hpp`).
   Asserts `<hydro> nscalars >= NSPECIES`.
2. **`src/chemistry/network_*.hpp`** — device-callable `KOKKOS_INLINE_FUNCTION` network:
   `RHS(y[NSPEC], n_H, T, zeta, G0, ydot[NSPEC])` returning d(abundance)/dt. Start with the
   **H2 network** (H, H2; grain formation kgr·n_H·y_H, CR + photodissociation destruction) —
   a direct port of `athena++/src/chemistry/network/H2.cpp` RHS. Extensible to a reduced
   gow17 (H2, H+, C+, CO, e−) reusing the ionization chemistry already in `ionization.hpp`.
3. **`src/chemistry/integrator.hpp`** — device `IntegrateCell`: explicit sub-cycled forward
   Euler with a chemical-time-step limiter `dt_chem = cfl·min(y/|ydot|)` (port of
   `forward_euler.cpp::GetChemTime`/`IntegrateOneSubstep`), capped at `nsub_max`. This is the
   GPU-tractable stand-in for CVODE; adequate for the moderately-stiff reduced network. A
   Rosenbrock-2 device solver is the upgrade path if stiffness demands it.
4. **`AddChemistryTasks(tc, pmesh, dt)`** — operator-split `par_for` over cells calling
   `IntegrateCell` on the `scalar_*` fields, run after RT in `hydro_driver.cpp` (reuse the
   OperatorSplit-flag exclusion already added for radiation so the scalars aren't double-updated).
5. **Coupling:** the network n_H comes from `prim_density`/(mu_n m_H); T from RT (if active) or
   the barotropic EOS; x_e from the network feeds back into the non-ideal coefficients (replace
   the standalone `ionization.hpp` equilibrium x_e with the evolved abundance when chemistry is on).

## Validation plan
- Unit: 1-zone relaxation to H2 equilibrium (compare to Athena++ H2 network steady state).
- Conservation: total H nuclei conserved to round-off under advection + reactions.
- Collapse smoke: tier-4 + `chemistry=true`, 1 H100, a few cycles, no NaN.

## Risk / scope notes
- **Flag-gated** (`<physics> chemistry=false` default) ⇒ zero impact on existing/queued runs;
  a new binary is backward-compatible.
- Explicit sub-cycling can blow up the substep count in cold dense gas (stiff) — `nsub_max`
  cap + abundance floors required (same lesson as the barotropic floors in `collapse_be`).
- Full gow17 (13 species, thermo + shielding) is the long-term target; the H2 + reduced-ion
  network is the first compilable, validatable milestone.

## Files to add / touch
- new: `src/chemistry/chemistry.{hpp,cpp}`, `network_h2.hpp`, `integrator.hpp`, `pgen/chem_relax.cpp` (test)
- edit: `src/hydro/hydro.cpp` (ProcessPackages branch), `src/hydro/hydro_driver.cpp`
  (AddChemistryTasks + OperatorSplit exclusion), `src/main.hpp` (enum if needed), `CMakeLists.txt`
