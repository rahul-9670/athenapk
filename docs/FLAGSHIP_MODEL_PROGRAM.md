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

## Phase 2 — Constrained transport (audit Workstream D.1) — **VALIDATED, GATE PASSED (2026-07-25)**

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
**PASSED (2026-07-25).** Matched pair `runs/inc7_gate/{glm_run_halloff,ct_run_halloff}`
(single-variable `divergence_control` toggle, full production physics minus Hall — see below),
both run to first-core formation `tlim=1.0945` at maxLevel 6. `runs/flux_retention.py` at the
deepest (first-core) snapshot:

| metric | CT | GLM | Δ |
|---|---|---|---|
| ρ_max/ρ_crit | 33.6 | 57.2 | — |
| M_core (M☉) | 0.0299 | 0.0316 | 5.4% |
| Φ_core (code) | 0.459 | 0.474 | 3.3% |
| μ_core (NN78) | 5.35 | 5.46 | 2.0% |
| retention Φ_core/Φ₀ | 2.48 | 2.44 | **1.6%** |

**Verdict:** the retained fossil flux is NOT set by the divergence-control choice (CT vs GLM
agree to 1.6% on retention, 2.0% on μ_core). Sub-finding: at fixed time CT reaches a *lower*
central density (34× vs 57× ρ_crit) — exact-div-free CT preserves more magnetic support (less
numerical field dissipation than GLM/Dedner cleaning). Not perfectly ρ-matched (both stopped
at tlim), so μ_core is the more robust comparison than raw density.

**CT+Hall blowup — root-caused and FIXED (2026-07-24/25).** The first full-physics CT gate leg
detonated at first core (t≈1.07747, level-4 boundary, ρ~4e4): |B|→1.35e15, dt→1e-55. Isolation
(GLM+Hall stayed healthy; CT+Hall-OFF ran clean to tlim; ratio=400 falsified) pinned it on the
**unsplit, uncapped Hall term** — η_H ∝ B/n_e explodes in the high-ρ recombination core, driving
a grid-scale dispersive whistler runaway that GLM cleaning damps but exact-div-free CT does not.
It was never a ∇·B failure (`ct_maxRelDivB` stayed round-off ~1e-13 throughout). **Fix:**
`diffusion/eta_hall_cap_code=0.05` (signed |η_H| clamp, = Ohmic-floor scale, < eta_ohm_cap 0.1;
keeps the unsplit-Hall dt ~ the hydro dt-wall at lvl6 while preserving Hall physics where η_H<0.05).
Validated at the exact detonation coordinates: `runs/inc7_gate/ct_run_hallcap` (Hall ON + cap,
single-variable delta from the run that blew up) cleared t=1.07747 at level 4 / ρ=4.57e4 with
|B|max=77 (vs uncapped 1.35e15), dt=5.6e-5 steady, 0 NaN — matching the healthy Hall-off
trajectory, not the blowup. **Canonical validated full-physics CT config:** `divergence_control=ct`,
`ct_emf=gs05`, `hall=hall` + `eta_hall_cap_code=0.05`, `rkl2_max_dt_ratio=400`, `eta_ohm_cap_code=0.1`
(gate run inputs live untracked in `runs/inc7_gate/`; this doc is the tracked reference).
NOTE: the *lean* ideal gate is impossible — adaptive-AMR + multigrid-gravity is unstable without
the full stack. The cross-code (Athena++) comparison uses the Hall-OFF pair since Athena++'s Hall
is a non-functional stub; the Hall-cap config is the AthenaPK-only most-complete-physics path.
**Risk:** high (touches the core update, AMR operators, restart). Largest single effort. **[done]**
**Unblocks:** a defensible *flux* claim.

## Phase 3 — Microphysics: EOS + conductivity (audit Workstreams A.1, A.2, D.2) — **COMPLETE (2026-07-25)**

**Closing summary (2026-07-25).** Both remaining items are resolved and the Phase-3 gate
("EOS round-trips + entropy-along-adiabat to tolerance; conductivities vs an independent
solver <1% off sign changes; charge neutrality + element conservation to roundoff") is MET.

*Conductivity — grain-inclusive cross-check DONE, and it found two real defects.* The MRN
grain-charge model is now reimplemented independently in Python
(`tests/test_conductivity.py` gates 6–8) via a genuinely different algorithm: the per-bin
capture balance has no bin dependence (the πa² cross-sections cancel between electrons and
ions), so ψ is one global unknown, Z_k = ψ τ_k, and the system reduces to a 1-D bracketed
root find in r = n_e/n_i — versus the C++ relaxed fixed point. It immediately exposed two
defects, both provable from the C++ output alone:
1. **`SolveCharges` did not converge** once grains dominate (ρ ≳ 1e-12): its output violated
   its own neutrality constraint n_i − n_e = Σ(−Z_k ng_k) by up to **3.6e7**, and returned
   n_e > n_i (impossible with negative grains). Fixed by the exact reduction above →
   bracketed bisection that cannot fail; neutrality residual **3.6e7 → 4.2e-11**.
2. **`SahaThermal` had an absolute-resolution floor**: fixed 64-step bisection on
   [0, n_K+n_H] ⇒ it cannot resolve below (n_K+n_H)/2⁶⁴ and returned ~1e-4 electrons in cold
   gas where the true value underflows — a spurious floor x_e ≈ 9e-20 (~9× the intended
   `xe_floor`) that grows with density. Fixed by relative-precision Newton from the
   weak-ionization limit.
Both behind `diffusion/ion_legacy_charge_solver` (default false = fixed). Final agreement
C++ vs independent Python: charge state and η_O to **5e-11**; η_H 4.2e-4 and η_A 2.9e-3,
cancellation-limited (σ_H's net is ~1e-13 of its summed |term| magnitude near the
grain-induced Hall sign reversal) and well inside the <1% criterion. **Result-changing but
narrow:** identical (0.0%) for ρ ≤ 1e-13 (the entire pre-first-core collapse) and again for
ρ ≥ 1e-9 (grains sublimated); it matters only in the grain-dominated first-core band
ρ ~ 1e-12…1e-10, where η_O shifts by up to 375% and **η_H changes sign** — a physics-level
correction, since the Hall sign sets the direction of field transport. Reproduce with
`tests/impact_charge_solver_fix.cpp`. **Reaches production on the next GPU rebuild.**

*EOS — hi-res table RESOLVED (swap is an input key, no code change).* Verified 2026-07-25:
`eos_table_hires.bin` (400×1000×920) passes ALL consistency gates and is **loader-compatible
without any C++ change** — the binary header is self-describing (nr, ne, nT) and the hi-res
file spans the identical (log ρ, log esp, log T) domain on finer *uniform* log grids, with
byte-exact layout (12 544 072 B = 72 + 8·(nr·ne·3 + nr·nT)). Measured side by side:

| table | P vs Saha (median / max) | cs² vs isentrope (median / max) |
|---|---|---|
| shipped 180×220×200 | 0.68% / **6.68%** | 0.85% / **17.8%** |
| hi-res 400×1000×920 | 0.16% / **0.55%** | 0.16% / **3.31%** |

⇒ 12× better on P, 5.4× better on cs², both at the H₂-dissociation / H-ionization kinks that
set first/second-core thermodynamics. **Adoption = one input line**,
`<hydro> eos_table_file=…/src/eos/eos_table_hires.bin` (default stays the shipped table so a
clean checkout always runs; the hi-res file is gitignored but regeneration was VERIFIED
bit-for-bit — `gen_eos_table.py build 400 1000 920 <out>.h5` reproduces md5
`6b8e3999eca19806d8f4d43054e0447c`). Deliberately NOT overwriting
`eos_table.bin`: it is git-tracked and read at RUNTIME by live jobs. The remaining
*non-uniform* grid idea is CLOSED on cost/benefit — it would require replacing the device
bilinear's even-log-grid assumption with a per-axis search in the EOS hot path, to buy only
the last cusp point (cs² median is already 0.16%).

### Original survey and increments (2026-07-23/24)

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

## Phase 4 — Radiation beyond grey M1 (audit Workstream C) — **COMPLETE (2026-07-26): multigroup RHD built + validated + GPU + campaign**

**DONE (2026-07-26, full multigroup RHD):** per-group (Er_g,Fr_g) M1 transport + group-coupled
implicit matter coupling (single Newton-on-T, GPU-safe) + physical band-mean opacity + the
state-of-the-art **tabulated dust+gas opacity** (`gen_opacity_table.py` → `opacity_table.bin`;
self-consistent Planck mean anchored to Bell&Lin, multi-species sublimation). Gates all PASS:
n_group=1 bit-identical to gray (max|old-new|=0), 3-group free-stream sum==gray 1.7e-14,
multigroup==gray coupling 7e-11, frequency reddening exp(−τ_g), restart+AMR safe, Er conservation
EXACT (2.14e-16 incl. AMR reflux), production-config validated (256³, 5×H100), GPU compile clean.
A bounded 256³ run validated + STOPPED; a **tabulated multigroup campaign to first core is running**
(`runs/mg_prod_tab/`). Commits 6d378d6..6ee898e. See `runs/validation_rt/MG_MULTIGROUP_RESULT.md`.
Not done (deferred, needs 2nd core / hybrid): mixed-frame v/c higher-order terms, ray-trace/MC
hybrid for direct irradiation, formal AP-diffusion limit proof. Original baseline notes below.

### Original baseline (2026-07-24)

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
not affect the nominal FHC run, but it is reachable by radiation in hot low-density gas. **Fix IMPLEMENTED (default-off, 2026-07-24):** `radiation_opacity.hpp::BellLinKappaFixed` + the
`<radiation> bell_lin_fix_regime_skip` flag (default false → the plain walk, bit-identical). The
fix SKIPS regimes whose forward window has closed [cross(j,j+1) ≤ cross(i,j)] and bridges regime i
directly to the next active regime. VERIFIED (test_opacity.cpp Gate 6): continuous everywhere
(max |d ln κ/d ln T| = 24 at every density incl. the defect corner), gives the correct canonical
values (0.02 @ 10 K, 2.0 @ 100 K, 0.348), and is **exactly identical to the plain walk on the
collapse track** (rho ≥ 1e-10, rel diff 0). NB: an earlier note claimed κ = min-over-all-8-regimes
was the fix — that was FALSIFIED (it gives 2e-31 instead of 0.02 at 10 K; wrong by ~29 decades;
only slope-continuity had been checked, not values). Reaches production on the next GPU rebuild with
the flag enabled; default-off keeps every current run bit-identical.
**MULTIGROUP RHD — STARTED (scaffold increment 1, 2026-07-24):** `radiation_groups.hpp` adds the
frequency-group structure (`RadGroups`: n_group + log-spaced edges) and the group-integrated
Planck function (`PlanckCumFraction` = Clark-1965 tail series; `PlanckFraction(g,T)`), designed so
**n_group=1 reduces EXACTLY to gray** (fraction 1) → production bit-identical. VERIFIED
(`test_rad_groups.cpp`): Planck series vs an independent numerical integral 3e-10; F(0)=0/F(∞)=1/
monotonic; **Σ_g PlanckFraction = 1 exactly** (energy partition conserved); gray reduction exact.
A bug (series was the upper tail, needed 1−tail) was caught by the gate. This is a NON-transport
foundation only; the remaining large build (per-group Er_g/Fr_g moments + per-group M1 closure +
per-group Planck/Rosseland opacity means + IMEX group coupling + mixed-frame v/c + AP diffusion
+ RSLA validity) is the multi-increment effort ahead — best done with fresh context, NOT rushed.
Other next increments: the M1 TRANSPORT benchmarks (free-streaming/shadow/diffusion/radiative-shock;
need running the C++ M1 solver on test setups, CPU build).

**Deliverables:** multigroup RHD (dust-IR → ionizing), consistent mixed-frame v/c terms,
scattering, IMEX/implicit transport for optically thick cells; hybrid moment + ray-trace/Monte-
Carlo for direct protostellar irradiation; asymptotic-preserving diffusion limit; local RSLA
validity check.
**Gate:** free-streaming/shadow/diffusion/radiative-shock benchmarks; grey-vs-multigroup and
M1-vs-hybrid differences quantified against the flux result. **Depends on:** Phase 3 opacities.

## Phase 5 — Grains & frequency-dependent opacity (audit Workstream B) — **opacity half DELIVERED (2026-07-26)**

**Frequency-dependent opacity DONE:** the tabulated monochromatic dust+gas opacity
(`gen_opacity_table.py`) with composition/size-motivated Planck & Rosseland band means, multi-
species sublimation cascade (ice ~150 K, refractory ~1500 K), and a self-consistent Planck mean —
the "size/composition-dependent means from one monochromatic dataset" deliverable (Phase 4 consumes
it). WS-4 grain growth (coagulation/fragmentation/sublimation) already exists in `src/dust/`.
**Remaining:** fully coupling an evolving grain-size/charge distribution INTO the opacity table per
cell (currently the table is a fixed dust model; the DustFactor consumer scales it). Grain-charge →
conductivity feed already covered by the MRN Wardle tensor (Phase 3). **Feeds:** 3 (done), 4 (done).

## Phase 6 — Multifluid where single-fluid fails (audit Workstream D.3) — **GATE DELIVERED (2026-07-26): single-fluid ADEQUATE, multifluid NOT needed**

**Single-fluid validity map DONE** (`src/hydro/diffusion/singlefluid_validity.py`, result in
`docs/PHASE6_SINGLEFLUID_VALIDITY.md`): computed from the deepest FHC snapshot (ρ up to 3.6e-8
g/cm³, the 2nd-collapse threshold), the ion strong-coupling parameter χ_ion=ν_in/ω_ff = 6.6e6…2e12
(min ≫ 1) EVERYWHERE ⇒ ion inertia negligible ⇒ **single-fluid non-ideal MHD (the Wardle tensor the
code already evolves) is adequate through the whole collapse; full multifluid is NOT on the critical
path.** Self-validated: β_i maps ambipolar→Hall→Ohmic with density (canonical FHC sequence). A
first pass using ν_ni was FALSIFIED (falsely flagged 80% multifluid — contradicts established AD
theory; the correct criterion is ν_in). **Remaining (only if a non-fiducial regime demands it):**
the actual separate-fluid evolution — deprioritized since the map shows it is unnecessary here.

## Phase 7 — Gravity/dynamics, ICs, ensembles, UQ, cross-code (audit E, F, K, J) — **gravity-BC gate DONE (2026-07-24)**

**Gravity-BC gate:** `src/self_gravity/tests/test_multipole.py` turns the WS-5a multipole-BC
validation into a committed reproducible test — checks `MultipolePhi` (monopole + traceless
quadrupole exterior potential) against the analytic identities: monopole −GM/r (1e-16) + 1/r
scaling, quadrupole P₂ structure + 1/r³ scaling, superposition linearity (2e-16), and the
traceless-Q shell-average → 0. ALL PASS. (Formula mirror of `multipole.hpp`; the header pulls
heavy parthenon deps so a host-g++ shim isn't worth it.) The rest of Phase 7 (torque diagnostics,
physics-based AMR, IC ensembles, UQ classes, cross-code) is new capability → user direction.

**Physics-based AMR DONE (2026-07-26, opt-in):** `refinement/nonideal.cpp` adds
`type = jeans_nonideal` — Jeans PLUS current-sheet resolution (refine where the field-reversal
scale L_B=|B|/|curl B| is under-resolved by < `curr_nsheet` cells, so the thin current layers that
set the reconnection/flux-loss are always resolved). Default `type=jeans` untouched (production
bit-identical). Validated on a 32³ MHD collapse: jeans stays 64 blocks; jeans_nonideal grows
64→127→155→169 (increments DECELERATE 63→28→14 ⇒ converges to a stable refined state, targeted on
current sheets, NOT runaway — the OOM risk is controlled). Caveat for production numlevel=20: a
pervasively-turbulent field could refine aggressively; monitor block counts and tune `curr_nsheet`.
The remaining Phase-7 IC-ensemble / UQ / cross-code pieces are the science program below.

**Torque / magnetic-braking diagnostic DONE (2026-07-26):** `runs/magnetic_braking.py` —
gradient-free (AMR-safe) rotational-support v_phi/v_Kep + convention-free B_phi/B_pol (braking
activity) shell profiles. On the deep prod_v9 snapshot: envelope sub-Keplerian (v_phi/v_Kep~0.4 at
1000 AU, infall-dominated), toroidal field wound up in the core (|B_phi/B_pol|~0.3 => braking
active). (Two physics bugs caught by falsification: Eulerian-vs-Lagrangian AM normalization, and
an M_enc->0 small-r artifact — inner <100 AU flagged as needing a proper disk-finder.) The
remaining Phase-7 pieces (physics-based AMR, IC ensembles, UQ classes, cross-code) below.

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

**Convergence machinery DONE (2026-07-26):** `runs/convergence.py` (harness) + `runs/
convergence_ladder/` (parameterized njeans={4,8,16} launcher + workflow README). Compares the
primary flux observables (mu_core=M_core/Phi_core, M_core, Phi_core, peak|B|) at matched PHYSICAL
STATE (matched rho_max, not time — reusing flux_retention.measure_snapshot), reports pairwise %diff
between the two finest + resolution-reach ("how far we can push"). Harness logic unit-tested
(recovers an injected 3% offset exactly). **Running the ladder is the compute step** (3 chained
campaigns to a matched 1e-13 first-core epoch), gated on GPU availability + user go — deferred with
the ensemble/UQ (Phase-3 plan: validate + know limits first, ensemble last).

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
