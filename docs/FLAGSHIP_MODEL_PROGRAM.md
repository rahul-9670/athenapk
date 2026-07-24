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

## Phase 1 — Unified microphysics & units foundation (audit Workstream A infra, D.2 prep) — **COMPLETE (2026-07-23)**

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
**Done:** `PhysicalUnits` (`src/units/physical_units.hpp`) + `IonizationEnvironment`
(`src/units/ionization_environment.hpp`) built and consumed by hydro/diffusion-ionization,
chemistry, radiation. Acceptance test `src/units/tests/test_physical_units.cpp` (host g++,
4/4: cgs round-trip 2e-16, `B²/2 ↔ B_cgs²/8π` bridge 6e-16, single T₀=10.0151 K, RSLA
c_code). Unit **provenance** stamped into every restart/output (13 `phys_units/*`
Restart-params in `hydro.cpp` Initialize; verified write + restart read-back).
**Units-authority resolved = Option A (2026-07-23):** the FHC base scales are DERIVED from
the BE IC (mass,temperature,f) via ONE shared function `DeriveBENormalization`
(`src/units/be_normalization.hpp`), used by both the collapse_be pgen AND `BuildPhysicalUnits`.
This eliminated the historical ~0.13% drift where the microphysics ran on rounded `<units>`
defaults (ρ 5.467e-19, l 2.81e16) while the dynamics used exact BE (ρ 5.4668e-19, l 2.8063e16);
`diffusivity_unit`/`opacity_unit` now match the dynamics exactly. **This is a one-time
re-baseline** of the AD/RT/chemistry coefficients (~0.13%) vs the prod_v9 lineage — reaches
production only on the next GPU rebuild; a `<units>` base-scale override for collapse_be is now
rejected, and a pgen guard asserts the two paths agree. **Risk:** low (centralization).
**Unblocks:** clean comparisons in all later phases.

## Phase 2 — Constrained transport (audit Workstream D.1) — **IMPLEMENTED, AT GATE (2026-07-23)**

**Deliverable:** staggered/face CT so discrete ∇·B = 0 to round-off, with AMR-compatible
divergence-preserving prolongation/restriction and reflux-curl; keep GLM as an optional
comparison path. **DONE** — branch `flagship-phase2-ct`, increments 1–4 (face Bf, edge EMF,
CT_UpdateBf curl, CT_ProjectBfToCC, RKL2 STS-curl for Ohmic/AD, unsplit Hall EMF). The
increment-4 Hall bug was a **task-graph race** (`CT_UpdateBf` missing the `emf` dependency,
dropping the last edge-EMF chain task on degenerate blocks) — fixed in 456b05a, whistler now
matches GLM. Full-production-physics smoke (chem+RT+non-ideal+AMR+gravity) PASSED under CT:
`ct_maxRelDivB ≈ 2e-14`, mass +0.007%, 0 NaN/floor over 50 cycles.
**Gate:** CT vs GLM flux-retention measured on the fiducial collapse — the result must not
be set by the divergence-control choice, or the difference is reported as systematic.
**IN FLIGHT** — matched gate pair `runs/inc7_gate/{glm_run,ct_run}` (= `fhc.in` + one-line
`divergence_control` toggle), full production physics to first-core formation (`tlim=1.0945`,
~cyc 804, ~2.5–3 h @ 4 GPU). Diagnostic = `runs/flux_retention.py`. NOTE: the *lean* ideal
gate is impossible — adaptive-AMR + multigrid-gravity is unstable without the full stack, so
the gate must run full production physics (GLM leg = re-baselined production, CT leg = new).
**Risk:** high (touches the core update, AMR operators, restart). Largest single effort. **[done]**
**Unblocks:** a defensible *flux* claim.

## Phase 3 — Microphysics: EOS + conductivity (audit Workstreams A.1, A.2, D.2) — **STARTED (2026-07-23)**

**Survey (evidence-first, 2026-07-23):** the core machinery already exists — a tabulated
**multi-Saha** EOS (`src/eos/eos_table.hpp` + `gen_eos_table.py`: H₂ dissociation, H & He/He⁺
ionization, H₂ rot/vib, inert He, in shared FHC units) and the full **Wardle conductivity
tensor** (`src/hydro/diffusion/ionization.hpp`: σ_O/σ_H/σ_P over e + i + charged grain bins,
signed η_H, two AD closures). Phase 3 is therefore a *consistency* phase, not greenfield.
**EOS consistency gate built + run** (`src/eos/tests/test_eos_consistency.py`, read-only over
the shipped `eos_table.bin`): charge neutrality **5e-15**, H/He element conservation exact,
and the fundamental free-energy identity `(∂u/∂v)|_T = T(∂P/∂T)|_ρ − P` holds to **2.5e-3**
(finite-difference-limited) — i.e. the underlying Saha EOS is **already thermodynamically
consistent**, so a full Helmholtz-free-energy rewrite is NOT needed for consistency.
**Actionable finding:** the shipped bilinear table is ~0.7% accurate typically but degrades to
**~5–7% right at the H₂-dissociation (~2200 K) and H-ionization (~7000–10000 K) γ-softening
kinks** — a *tabulation-resolution* limit at exactly the transitions that set first/second-core
thermodynamics.

**EOS kink-fidelity increment (2026-07-24):** resolution sweep (worst-case P error vs Saha):
180×220×200 (shipped) 6.85% → 240×520×480 2.72% → 300×760×700 1.26% → **400×1000×920 0.80%**
(median 0.12%, 12.5 MB, ~4.5 min to generate). Convergence is only ~O(h) at the near-cuspy
γ-drops, so uniform refinement reaches <1% but not much better without a *non-uniform* grid
(which needs a C++ loader change — the device loader assumes an even log grid — DEFERRED,
needs GPU rebuild). `gen_eos_table.py` now takes `build <nr> <ne> <nT> [out.h5]`; a plain
`build` still reproduces the shipped `eos_table.bin` **bit-for-bit** (md5 07e5423b…), so
production is untouched. Validated hi-res candidate = `src/eos/eos_table_hires.bin`
(gitignored, regenerable); passes ALL consistency gates. A new **Gate 5 (sound-speed /
entropy-along-adiabat consistency)** — the tabulated cs2 vs an independent Saha isentrope-slope
computation — shows cs2 is the more derivative-sensitive quantity: **17.8% worst at the
H-ionization cusp in the shipped table → 3.3% in the hi-res** (median 0.85% → 0.16%).
**DEFERRED / user-gated:** swapping production onto the hi-res table is result-changing AND the
table is loaded at RUNTIME (the currently-running GLM gate job is reading `eos_table.bin` live),
so the swap must be coordinated with the user (and a matched GLM/CT gate pair).

**Conductivity gate (2026-07-24):** `src/hydro/diffusion/tests/test_conductivity.py` — an
INDEPENDENT pure-Python reimplementation of the Pandey & Wardle (2008) tensor + η_{O,H,A} that
`ionization.hpp::Diffusivities` computes, evaluated on the same gas-phase CR↔recombination
balance. Validates: **σ_O == scalar conductivity Σ n_j Z_j² e²/(m_j ν_j) to 4e-16** (the B
cancels exactly — guards the β/ecB assembly), η_O = c²/(4πσ_O), **exact parity in B** (η_O,η_A
even; η_H odd — the Hall sign structure), σ_⊥² = σ_H²+σ_P², and the **canonical ambipolar→Ohmic
crossover** (η_A/η_O = 2.7e6 at 1e-17 g/cm³ → 2.7e-3 at 1e-8). ALL PASS.
**Direct C++-vs-Python numerical cross-check DONE** (`tests/xcheck_conductivity.cpp` +
`tests/host_shims/`): compiles the REAL `ionization.hpp` host-only and dumps η_{O,H,A}; the
shipped device tensor matches the independent Python reference to **η_O 3e-7, η_H 4e-7, η_A
2.4e-3** — all < the Phase-3 <1% criterion (η_A is cancellation-limited: `σ_P/σ_⊥² − η_O` with
η_A ≪ η_O at the worst point). Grains OFF (`f_dg=0`) so it matches the gas-phase reference.
**Remaining Phase-3 conductivity item:** the grain-INCLUSIVE cross-check needs `SolveCharges`'
MRN grain-charge model reimplemented in Python (the C++ model's extension beyond gas phase).

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

## Phase 4 — Radiation beyond grey M1 (audit Workstream C) — **BASELINE CHARACTERIZED (2026-07-24)**

**Survey:** the radiation package (`src/radiation/`) is **gray single-group M1** — Er/Fr moments,
Artemis-ported M1 closure (`radiation_closure.hpp`), explicit HLL transport (`radiation_moments.cpp`),
implicit matter coupling, gray Planck/Rosseland opacity (`radiation_opacity.hpp`, WS-3a split).
**No multigroup structure exists** — so unlike Phase 3, the core Phase-4 deliverable is genuinely
new, large, high-risk numerics (NOT to be greenfielded unilaterally in an autonomous run).
**Baseline gate DONE:** `src/radiation/tests/test_m1_closure.cpp` (+ `host_shims/`) compiles the
REAL closure host-only and validates every M1 identity — Eddington-factor limits χ_Edd(0)=1/3 /
χ_Edd(1)=1, Eddington-tensor trace ≡ E (2e-16), beam (P=n_i n_j) & diffusion (P=⅓δ) limits, causal
wave speeds (|λ|≤1; 1/√3 diffusion, 1 beam), flux clamping. ALL PASS.
**Opacity gate DONE** (`src/radiation/tests/test_opacity.cpp`): compiles the REAL
`radiation_opacity.hpp` host-only; validates Bell&Lin continuity on the collapse track, canonical
dust values (κ=0.02 at 10 K, 2.0 at 100 K; e-scatter floor 0.348), the T⁻²⁴ sublimation gap, and
the Planck/Rosseland split (ratio=1 → κ_P≡κ_R). ALL PASS. **FOUND A LATENT DEFECT (deferred):** at
ρ ≤ 1e-11 g/cm³ and T~4–9 kK the Bell&Lin transition temperatures go out of order and the regime
walk **skips regime 6 (Kramers), giving a κ discontinuity of up to ~4 decades** (regime 5→7). This
corner is OFF the cold-collapse track (cold at low ρ; high T only at ρ≥1e-10 second core) so it does
not affect the nominal FHC run, but it is reachable by radiation in hot low-density gas. The fix
(monotonic regime selection / min-of-regimes) is **result-changing to production opacity → DEFERRED,
user-gated**.
**DEFERRED (user-gated, large):** multigroup RHD + mixed-frame v/c + scattering + IMEX +
asymptotic-preserving diffusion + RSLA validity are the real Phase-4 build. The next tractable
increments are the transport benchmarks (free-streaming/shadow/diffusion/radiative-shock), which
need running the C++ M1 solver on test setups (CPU build).

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

**Immediate next actions (this program's live edge, updated 2026-07-23):**
1. Read the running Phase-2 gate result (`runs/inc7_gate/{glm_run,ct_run}`), run
   `flux_retention.py`, and record whether divergence control sets the flux number — closes Phase 2.
2. Land the `PhysicalUnits` package (Phase 1) and migrate consumers off duplicated calibrations —
   **this is now the next coding edge** once Phase 2's gate is read.
3. (Phase 0 tail) full provenance-into-restart capture + the machine-readable test suite.

See `~/HANDOFF_2026-07-23.md` for the live operational state (gate jobs, prod_v9 hold, binaries).
