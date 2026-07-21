# Flagship Fossil-Field Model — Program Plan

Status: **launched 2026-07-21**. This document turns the independent audit's aspirational
"state-of-the-art roadmap" (Workstreams A–M in `INDEPENDENT_CORRECTNESS_AUDIT_2026-07-21.md`)
into a concrete, dependency-ordered, staged coding program for AthenaPK, scoped to the
science goal: **a converged, uncertainty-quantified prediction of the fossil magnetic flux
inherited by a protostar from the collapse of a magnetized Bonnor–Ebert sphere.**

## Honest framing

The full roadmap is a multi-year, multi-person research program; no single session or
person completes it. What is tractable — and what this document commits to — is:

1. a **staged decomposition** with explicit dependencies, so each phase is a shippable,
   independently-validatable unit rather than a monolith;
2. a **critical-path ordering** driven by what actually moves (and de-risks) the headline
   flux-retention number, not by feature completeness;
3. **hard validation/publication gates** between phases (a successful collapse is not
   validation);
4. a **begun Phase 1** (the unified microphysics/units foundation) — real code, this session.

Guiding rules inherited from the audit and the project working-style rules:
- Correctness before physics: no new model lands while a confirmed high-severity defect is open.
- One variable per experiment: never combine a numerics fix and a model change in one comparison.
- Every change gets a small analytic/unit test, a subsystem benchmark, and a coupled-collapse regression.
- Caps/floors/reduced-c/frozen-eta are **model changes**, never "free" optimizations.
- Pre-register quantitative acceptance criteria before looking at the preferred result.

## Where the code is today (post-audit)

- **Solver:** GLM-MHD (Dedner cleaning) + multigrid self-gravity + M1 grey radiation +
  gow17-reduced/H2 chemistry + Ohmic/AD/Hall non-ideal MHD + tabulated protostellar EOS,
  on Parthenon AMR / Kokkos (Hopper GPU).
- **Audit correctness fixes landed (branch `audit-fixes-2026-07-21`):** ionization-rate
  unification (#1), stage-consistent gravity (#2), non-silent radiation convergence (#6),
  cooling task-dependency race (#7), radiation coarse-fine reflux (A1), post-coupling
  radiation ghost refresh (A2), post-chemistry scalar ghost refresh before STS (A3), and
  the units.json B-field √(4π) source fix (#3). These are the substrate the flagship
  builds on; they must be GPU-validated and re-baselined before Phase 2 begins.

## Critical path for the flux-retention number

Ranked by how directly each unresolved item can change or invalidate the headline:

1. **Divergence control** (GLM vs constrained transport). Numerical ∇·B transport can
   masquerade as physical flux redistribution — this is the single largest methodological
   risk to a *flux* result. → Phase 2.
2. **Conductivity microphysics consistency** (chemistry/grains → η_O, η_H, η_A from one
   charge population). The non-ideal terms set how much flux is lost. → Phase 3.
3. **Thermodynamics** (free-energy EOS through H₂ dissociation / ionization) — sets the
   collapse's pressure support and thus the field amplification. → Phase 3.
4. **Radiation realism** (grey M1 → multigroup + hybrid) — sets the thermal structure that
   feeds ionization and EOS. → Phase 4.
5. **Turbulent + environmental ensemble and UQ** — one seed is not a prediction. → Phase 7.

Everything sits on a shared microphysics/units foundation → **Phase 1 first.**

---

## Phase 0 — Reproducible baseline + correctness (audit Phases 0–3)

**Deliverables:** the seven correctness fixes (done); provenance capture (source + submodule
revisions, `CMakeCache.txt`, compiler/CUDA/MPI/HDF5/Kokkos config, executable hash, fully
expanded runtime input, EOS/opacity/chemistry table hashes) written into every restart/output,
not just a sidecar; the small restartable test-configuration suite (advection, MHD waves, GLM
pulse, Ohmic/AD/Hall, Jeans, isolated sphere, radiation limits, EOS round-trips, reduced-res
collapse) emitting machine-readable norms.
**Gate:** all seven fixes GPU-validated; test suite green; a re-baselined collapse whose
flux-retention differs from the old scheme by a *quantified, understood* amount.

## Phase 1 — Unified microphysics & units foundation (audit Workstream A infra, D.2 prep) — **STARTED**

The substrate every later phase consumes; also the permanent root-fix for audit #3/#4.
**Deliverables:**
- **`PhysicalUnits`** — one immutable package built before all physics packages, storing
  base mass/length/time/temperature/density/energy-density/**magnetic**/opacity/diffusivity
  units, with `B_unit = sqrt(4π·ρ_unit)·v_unit` (Heaviside-Lorentz). Consolidates the
  scattered per-package calibrations (`chemistry` rho_unit/t_unit/mu_n/T_unit,
  `diffusion` ion_*_unit, `radiation` mu, `collapse_be` code_*_cgs). Written to every
  restart/output. Roots #3 (correct B everywhere) and #4 (single T₀; kills the mu=2.33
  vs 2.29 split).
- **`IonizationEnvironment`** — one struct owning cosmic-ray rate + attenuation,
  radionuclide floor, composition, grain population, thermal-ionization params; consumed by
  chemistry AND every non-ideal coefficient (generalizes the #1 fix from an init-time
  assertion to a shared object).
**Started this session:** #1 input-level unification + init consistency check; #3 pgen
source fix. **Next:** the `PhysicalUnits` header + package and migration of consumers;
consistency acceptance tests (1000-sample cgs round-trip, `B²/2 ↔ B_cgs²/8π`, identical
physical T across radiation/EOS/chemistry/ionization for identical code state).
**Risk:** low (centralization, not new numerics). **Unblocks:** clean comparisons in all
later phases.

## Phase 2 — Constrained transport (audit Workstream D.1)

**Deliverable:** staggered/face CT so discrete ∇·B = 0 to round-off, with AMR-compatible
divergence-preserving prolongation/restriction and reflux-curl; keep GLM as an optional
comparison path.
**Gate:** CT vs GLM flux-retention measured on the fiducial collapse — the result must not
be set by the divergence-control choice, or the difference is reported as systematic.
**Risk:** high (touches the core update, AMR operators, restart). Largest single effort.
**Unblocks:** a defensible *flux* claim.

## Phase 3 — Microphysics: EOS + conductivity (audit Workstreams A.1, A.2, D.2)

**Deliverables:** thermodynamically-consistent free-energy H/He(+metals) EOS (rot/vib H₂,
ortho/para, dissociation, ionization stages, degeneracy, Coulomb, radiation) with all
thermodynamic derivatives from one Helmholtz free energy (Maxwell-consistent); non-equilibrium
thermal chemistry with a robust implicit stiff solver (sparse Jacobian, positivity, exact
element/charge conservation); full Ohmic/Hall/Pedersen conductivity tensor from all charged
species + grain bins, sign changes retained (not clipped), coupled to the same chemistry/grain
model. Published microphysics tables with uncertainty bands.
**Gate:** EOS round-trips + entropy-along-adiabat to tolerance; conductivities vs an
independent solver <1% off sign changes; charge neutrality + element conservation to roundoff.
**Risk:** high. **Unblocks:** trustworthy non-ideal flux loss and collapse thermodynamics.

## Phase 4 — Radiation beyond grey M1 (audit Workstream C)

**Deliverables:** multigroup RHD (dust-IR → ionizing), consistent mixed-frame v/c terms,
scattering, IMEX/implicit transport for optically thick cells; hybrid moment + ray-trace/Monte-
Carlo for direct protostellar irradiation; asymptotic-preserving diffusion limit; local RSLA
validity check.
**Gate:** free-streaming/shadow/diffusion/radiative-shock benchmarks; grey-vs-multigroup and
M1-vs-hybrid differences quantified against the flux result. **Depends on:** Phase 3 opacities.

## Phase 5 — Grains & frequency-dependent opacity (audit Workstream B)

Evolve grain size + charge distribution (coagulation/fragmentation/sublimation/ice/charging);
composition- and size-dependent Planck/Rosseland/flux/energy means from one monochromatic
dataset. **Feeds:** Phases 3 (grain charge → conductivity) and 4 (opacity).

## Phase 6 — Multifluid where single-fluid fails (audit Workstream D.3)

Diagnose collision frequencies / drift / gyrofrequencies / coupling throughout the run; where
single-fluid breaks, evolve ions/electrons/charged grains separately or a generalized
multifluid Ohm law, with drag heating. **Gate:** single-fluid validity map; agreement with
multifluid where both valid.

## Phase 7 — Gravity/dynamics, ICs, ensembles, UQ, cross-code (audit E, F, K, J)

Stage-centered gravity with validated isolated/Green-function BCs (Phase 0's #2 is the first
step); torque-budget diagnostics; physics-based AMR (Jeans + pressure scale height + non-ideal
lengths + current sheets + disk scale height). Replace the single BE sphere with a controlled
ensemble (BE reference + cloud-extracted cores) sampling mass/rotation/magnetization/alignment/
turbulence/metallicity/CR/pressure. Separate uncertainty classes (discretization, AMR/boundary,
RSLA, microphysics, IC variance, subgrid, analysis); Latin-hypercube/Bayesian design; report
predictive *distributions* for flux retention, not one curve. Reproduce the fiducial collapse
with ≥1 independent RMHD non-ideal code; blinded analysis.
**Gate:** primary observables converge <5% (or credible extrapolation) between the two finest
resolutions; ensemble scatter quantified; two-code agreement within combined uncertainty.

## Phase 8 — Protostar & observational connection (audit G, L)

Continue through second-core accretion to a physically-defined protostellar surface (couple a
stellar-evolution/accretion-shock model, magnetically-launched outflows, irradiation); track
flux crossing the surface vs retained in disk vs expelled vs dissipated; distinguish inherited
large-scale flux from dynamo field. Forward-model polarized dust emission / Zeeman / RM /
outflow polarization vs observed ensembles.

## Staged publication sequence (audit)

1. Methods paper (verified CT-RMHD-non-ideal + microphysics + unit system + test suite).
2. Controlled BE-sphere paper (resolution/boundary/RSLA/microphysics convergence + turbulent ensemble to second core).
3. Environmental-ensemble paper (cloud-extracted cores + observational priors + UQ).
4. Protostellar-survival paper (inherited flux through accretion/convection/outflow/dynamo).

## Flagship readiness (final gate, audit)

Two independent codes agree within combined uncertainty on a shared benchmark; flux observables
spatially+temporally converged or credibly extrapolated; CT-vs-alternate divergence treatment
shows divergence control does not set the result; gravity BC/domain/Poisson-tolerance converged;
RSLA and grey/multigroup and M1/hybrid differences bounded; conductivity uncertainties
propagated; cap/floor activations disclosed and bounded; adequate turbulent/environmental
ensemble; field followed to a physical surface (not a sink); numerical/microphysical/sample
uncertainties reported separately; conclusions survive alternate analysis surfaces.

---

**Immediate next actions (this program's live edge):**
1. GPU-validate the seven audit fixes (job 2375916 build → smoke) and re-baseline the collapse.
2. Land the `PhysicalUnits` package (Phase 1) and migrate consumers off duplicated calibrations.
3. Scope Phase 2 (CT) design note before touching the core update.
