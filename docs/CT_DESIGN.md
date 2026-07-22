# Phase 2 — Constrained Transport (CT) Design & Scoping

**Progress (implementation):**
- **Increment 1 (ideal, single level) — DONE.** `divergence_control = glm|ct` switch, face field
  `Bf`, Balsara–Spicer edge EMF, Stokes curl, face→CC projection, `ct_maxRelDivB` history
  diagnostic. GLM path bit-identical. Field-loop / OT / div-pulse div-free to round-off
  (~1e-13). Commit `1d62eff`.
- **Increment 2 (AMR reflux-curl) — DONE (2026-07-22).** The edge EMF is now assembled *before*
  the flux-correction round in `hydro_driver.cpp`, so the existing Load/Receive/Set flux-correction
  trio (which handles a Face variable's edge fluxes via `GetFluxCorrectionElements`) restricts the
  fine-block edge EMFs onto the coarse neighbour's shared edge — both sides curl the identical EMF,
  keeping ∇·B at round-off across coarse–fine boundaries. Validated on a static-AMR field loop whose
  C–F boundary (x1=0, 20 blocks) bisects the loop: `ct_maxRelDivB` max = **7.6e-13** (round-off).
  Falsification: the increment-1 ordering (EMF assembled *after* correction, no reflux) spikes to
  **5.5** on the identical mesh. Test = `inputs/field_loop_ct_amr.in`, suite = `runs/ct_tests/run_ct_tests.sh`.
  The Tóth–Roe div-free prolongation ops were already registered on `Bf` in increment 1.
- **Increment 3 (GS05 upwind EMF) — NEXT.** Replace the arithmetic edge-EMF average with the
  Gardiner & Stone (2005) upwinded reconstruction (field-loop dissipation gate; matches Athena++).

**Status:** design / scoping note (read-only investigation, no code changed).
**Scope:** replace GLM/Dedner hyperbolic divergence cleaning with staggered/face-centered
Constrained Transport so discrete ∇·B = 0 to round-off, keeping GLM as a selectable
comparison path. Corresponds to `FLAGSHIP_MODEL_PROGRAM.md` **Phase 2 / audit Workstream D.1**.
**Why it is the top numerical risk to the science:** the headline result is a *magnetic
flux* inherited by the protostar. GLM cleaning transports and dissipates ∇·B errors
advectively; that numerical transport is indistinguishable, in a flux diagnostic, from
physical flux redistribution. A defensible fossil-flux number requires that the divergence
control provably does not move flux — i.e. ∇·B = 0 by construction (CT), cross-checked
against GLM.

Provenance of the investigation this note is built on (read this session):
- Parthenon submodule `external/parthenon` @ `fe2262799` ("flux-correction comm fence fix,
  PR #1405").
- AthenaPK current MHD: `src/main.hpp` (field layout), `src/hydro/hydro.cpp`,
  `src/hydro/hydro_driver.cpp`, `src/hydro/rsolvers/glmmhd_hlld.hpp`,
  `src/hydro/glmmhd/dedner_source.cpp`, `src/hydro/diffusion/{diffusion,resistivity}.cpp`,
  `src/hydro/prolongation/custom_ops.hpp`, `src/pgen/collapse_be.cpp`.

---

## 0. Executive summary

- **The single largest effort driver — a face-field / CT AMR framework — already exists in
  Parthenon.** Face and Edge topological variables, staggered storage and packing,
  divergence-free prolongation (Tóth & Roe 2002), shared-face prolongation and restriction,
  a generalized Stokes-theorem (curl) update, face/edge ghost exchange, and EMF/flux
  correction across coarse–fine boundaries are all present, with a working reference in
  `external/parthenon/example/fine_advection` (its `do_CT_advection` path). This removes the
  highest-risk, highest-effort piece (AMR-consistent divergence-free operators) from
  AthenaPK's plate.
- **What AthenaPK must build is the MHD-specific glue:** a face-centered B field, the CT EMF
  built from the existing HLLD face fluxes, the curl update onto faces, a face→cell-center B
  projection so the rest of the code (Riemann, EOS, non-ideal, radiation, gravity, analysis)
  is unchanged, re-routing of the non-ideal (Ohm/AD/Hall) EMFs from cell-face fluxes to
  shared edge EMFs, IC/restart plumbing, and a GLM/CT switch.
- **Recommended flavor:** Gardiner & Stone (2005, 2008) upwind CT ("flux-CT" family), with
  simple arithmetic (Balsara–Spicer 1999) EMF averaging as the first-light fallback. Reasons:
  it reuses the existing HLLD face fluxes with minimal new solver work; it is the documented
  scheme in the Athena++ reference code, aiding the head-to-head cross-code comparison the
  project is built on; and it maps cleanly onto Parthenon's Stokes/Tóth–Roe machinery.
- **Effort/risk:** High — largest single item on the flagship path. It touches the core
  update, the induction-equation data model, AMR operators, and the restart format. The
  Parthenon support drops the estimate from "port a full staggered-MHD AMR infrastructure"
  to "wire an EMF and a projection into an existing staggered framework," but the coupling
  surface (every consumer of cell-centered B) is wide and must be handled by the projection
  discipline below. Rough order: **6–10 focused increments**, each independently validatable.
- **Caveat found in the framework:** Parthenon's Stokes curl path for *Face* variables
  carries an in-source `TODO(LFR): This is untested, need to test in parthenon-mhd downstream
  or add a test involving curls` (`example/fine_advection/stokes.hpp:177`). The scaffolding
  exists but the face-curl sign/geometry has not been battle-tested inside Parthenon; part of
  our validation burden is to be its first real MHD exercise.

---

## 1. Current state — how GLM-MHD is structured today

### 1.1 Field layout (cell-centered B)
`src/main.hpp`:
```
NHYDRO = 5, IB1 = 5, IB2 = 6, IB3 = 7, IPS = 8   // IPS = GLM scalar psi
```
B is three **cell-centered** components living inside the single `cons` variable (and its
`prim` twin), registered in `src/hydro/hydro.cpp` (~L1267):
```cpp
Metadata m({Metadata::Cell, Metadata::Independent, Metadata::FillGhost, Metadata::WithFluxes}, ...);
m.RegisterRefinementOps<refinement_ops::ProlongateCellMinModMultiD, ...>();
pkg->AddField("cons", m);
```
So today B is prolongated/restricted by the **generic limited cell-centered** operator
`ProlongateCellMinModMultiD` (`src/hydro/prolongation/custom_ops.hpp`) — which does **not**
preserve ∇·B across coarse–fine interfaces. That operator, however, is already written to be
topological-element-aware (`INCLUDE_X1/2/3` switch on `TE::CC/F1/F2/F3/E1/E2/E3`), so it can
serve face fields for the *shared*-face part if desired; internal faces need Tóth–Roe.

### 1.2 Hyperbolic update path
`src/hydro/hydro_driver.cpp` per stage: reconstruct → `calc_flux_fun` (Riemann) →
`StartReceive/LoadAndSend/Receive/SetFluxCorrections` → `Update::FluxDivergence` →
weighted-sum stage update → `AddSplitSourcesStrang` / Dedner source → boundary exchange →
`FillDerived` (cons→prim). Flux correction across coarse–fine boundaries for `cons` is
already wired (`parthenon::{LoadAndSend,Receive,Set}FluxCorrections`).

The GLM HLLD solver `src/hydro/rsolvers/glmmhd_hlld.hpp` returns, per x-face, the fluxes of
all conserved variables including the transverse magnetic components and the GLM pair:
```cpp
flxi[IB1] = psii;            // = 0.5(psiL+psiR) - 0.5 c_h (BxR-BxL)
flxi[IPS] = c_h^2 * bxi;     // bxi = 0.5(BxL+BxR) - 0.5/c_h (psiR-psiL)
... flxi[IB2] = ...by...; flxi[IB3] = ...bz...;   // transverse-B fluxes = ∓EMF components
```
**Key fact for CT:** the transverse-B face flux *is* an EMF component. In flux form the
induction equation stores, at an x-face, `F_x(B_y) = -E_z` and `F_x(B_z) = +E_y`. CT reuses
exactly these numbers; it only changes *where they are collocated* (cell face → cell edge)
and how they are differenced (flux-divergence → curl).

### 1.3 Dedner GLM source
`src/hydro/glmmhd/dedner_source.cpp` — parabolic ψ damping `psi *= exp(-alpha c_h dt/mindx)`
plus, in the `extended` variant, the non-conservative `-(∇·B) B` momentum and
`-B·∇ψ` energy source terms (Dedner+/Mignone–Tzeferacos). `c_h` is the global cleaning speed.
Under CT the entire ψ field, its flux, and this source term are **dropped** (CT path has no
ψ); this file becomes GLM-path-only.

### 1.4 Non-ideal EMFs (Ohm / AD / Hall)
`src/hydro/diffusion/` computes edge-current-based EMFs and **adds them to the same
cell-face induction fluxes**, e.g. `src/hydro/diffusion/resistivity.cpp`:
```cpp
cons.flux(X1DIR, IB2, k,j,i) += -eta * j3;   // -E_z at x-face
cons.flux(X1DIR, IB3, k,j,i) +=  eta * j2;   // +E_y at x-face
cons.flux(X1DIR, IEN, ...) += Poynting...
```
So AthenaPK already runs a **flux-CT-of-cell-centered-B** discretization: induction is a set
of face fluxes, ideal (HLLD) + non-ideal (Ohm/AD/Hall) summed at the same face, then
flux-differenced onto cell-centered B. This is the crucial structural asset — moving to true
staggered CT is a *re-collocation* of an EMF the code already assembles, not a new physics
derivation. AD works under RKL2 STS; Hall is unsplit-only and dispersive (see §5.4).

### 1.5 Initial condition
`src/pgen/collapse_be.cpp` sets a uniform vertical field cell-centered (`u(IB3)=B0z`,
`u(IB1)=u(IB2)=0`). Trivially divergence-free, so the face initialization is exact and the
IC is the easiest possible starting point for CT validation.

---

## 2. Parthenon CT support — what exists vs. what must be built

### 2.1 EXISTS in the framework (verified this session)

| Capability | Where | Notes |
|---|---|---|
| Face & Edge variable types | `src/interface/metadata.hpp` (`Metadata::Face`, `Metadata::Edge`); `variable.hpp:153` maps them to `{TE::F1,F2,F3}` / `{TE::E1,E2,E3}` | First-class staggered storage; a Face field carries 3 components (one per face orientation). |
| Topological packing/indexing | `pack(b, TE::F1, var(), k,j,i)` (see `example/fine_advection/advection_package.cpp:252+`); `PackVariablesAndFluxes` | Read/write staggered data with the same pack API as CC. |
| Divergence-free internal prolongation | `src/prolong_restrict/pr_ops.hpp:391` `ProlongateInternalTothAndRoe` | **Tóth & Roe (2002)** div-free face prolongation — the hard AMR piece, done. `OperationRequired` fires only for Face fields with `cel==CC`. |
| Shared-face prolongation + restriction | `pr_ops.hpp` `ProlongateSharedMinMod`, `RestrictAverage`, `ProlongateInternalAverage` | Complete the AMR operator set for faces. |
| Register the operator triple on a field | `m.RegisterRefinementOps<ProlongateSharedMinMod, RestrictAverage, ProlongateInternalTothAndRoe>()` (`advection_package.cpp:102`) | One line to make a Face field div-free-AMR-correct. |
| Generalized Stokes (curl) update | `example/fine_advection/stokes.hpp` `Stokes()/StokesComponent()` | Maps Face-variable ↔ Edge-flux (EMF), does the curl with correct submanifold sign handling. **Carries the untested-TODO for the Face path (see §0).** |
| Face/Edge ghost exchange | `Metadata::FillGhost` on a Face field; `bnd_info.cpp` handles `topo_comp=3` for Face/Edge | Halo exchange of staggered data works. |
| EMF / flux correction across coarse–fine | `src/bvals/comms/bnd_info.cpp:77` `GetFluxCorrectionElements` handles Face variables and their **edge** fluxes; `AddFluxCorrectionTasks` (`boundary_communication.cpp:470`) | This is the **reflux-curl / shared-EMF** step — a fine block's edge EMFs are averaged onto the coarse neighbor's shared edge so both integrate the *same* EMF ⇒ ∇·B stays zero at the C–F boundary. Already generic. |
| End-to-end reference | `example/fine_advection` `do_CT_advection` (`advection_driver.cpp:137`) | Divergence-free advection of a Face vector with edge "vector fluxes," Tóth–Roe AMR, and flux correction. This is our template. |

### 2.2 MUST BE BUILT in AthenaPK (framework does not provide)

1. A face-centered B field `Bf` (`Metadata::Face`) and the decision to make it (not `cons`'s
   IB1..3) the **primary/independent** magnetic variable on the CT path.
2. The **CT EMF**: combine the HLLD (and non-ideal) *face* fluxes of transverse B into
   *edge*-centered EMFs (the averaging is the physics choice — §3.2).
3. The **curl update** wiring Bf ← curl(EMF) via the Stokes machinery, inside the AthenaPK
   task list.
4. A **face→cell-center projection** `Bf → (IB1,IB2,IB3)` executed every substage before the
   Riemann solve / FillDerived, so every existing consumer of cell-centered B is untouched.
5. Re-routing the **non-ideal EMFs** to the shared edge EMF instead of the cell-face flux.
6. **IC**: initialize `Bf` on faces (from a vector potential in general; trivial for the
   uniform-B0z BE sphere), then project to CC.
7. **Restart/output**: `Bf` added to the restart set; phdf science output keeps CC B for
   analysis (derived), so notebooks are unaffected.
8. A **GLM ⟷ CT switch** (see §3.6) so both paths coexist for the comparison gate.

---

## 3. Design — staged CT plan

### 3.1 Data model
- Add `Bf` = `Metadata({Face, Independent, WithFluxes, FillGhost})` with
  `RegisterRefinementOps<ProlongateSharedMinMod, RestrictAverage, ProlongateInternalTothAndRoe>()`.
- **Keep** cell-centered `IB1..IB3` inside `prim`/`cons` as *derived* on the CT path (still
  `WithFluxes` is unnecessary for them on CT — the induction flux moves to Bf's edge flux —
  but keeping the slots avoids touching NHYDRO indexing and every kernel that reads B). The
  clean split: on CT, `IB1..3` are **outputs of a projection**, never independently evolved;
  their `cons` flux entries are ignored (or zeroed) on the induction side.
- Drop `IPS`/ψ on the CT path (compile/runtime-gated). Simplest first cut keeps the index
  slot allocated but unused to avoid re-numbering; a later cleanup can compact it.
- Edge EMF storage: either a `Metadata::Edge` derived field `emf` (3 components on
  E1/E2/E3), or use Bf's `flux` slots as the fine_advection example does (`CalculateVectorFluxes`
  writes the edge flux of the face field). Prefer the **Bf-flux** route to reuse
  `AddFluxCorrectionTasks` unchanged.

### 3.2 CT EMF construction (the physics choice)
Per stage, after the directional HLLD solves have populated the cell-face fluxes of
transverse B (already computed today), assemble the edge-centered EMF. Two documented options,
staged:

- **Stage A (first light) — simple arithmetic averaging (Balsara & Spicer 1999):**
  `E_z(edge) = ¼( F_x(B_y)|_{j} + F_x(B_y)|_{j+1} − F_y(B_x)|_{i} − F_y(B_x)|_{i+1} )` with
  signs from the flux↔EMF identities in §1.2. Trivial to implement from existing face fluxes;
  correct ∇·B = 0; known to be slightly under-dissipative on grid-aligned discontinuities.
- **Stage B (production) — Gardiner & Stone (2005) upwind CT:** replace the plain average
  with the directionally-upwinded reconstruction of the EMF derivatives (the GS05 "integration
  constant" chosen from the contact-mode upwind direction), which removes the field-loop
  advection artifacts and matches Athena++'s scheme. This is the recommended production flavor
  because the reference code (Athena++) uses the same family, tightening the cross-code
  comparison, and because it is the best-tested upwind-CT for this solver family.

Non-ideal contributions (Ohm/AD/Hall) are added into the **same edge EMF** before the curl
(see §3.5), so ideal and non-ideal flux transport share one divergence-free update.

### 3.3 Curl / Stokes update
Update `Bf` by the discrete curl of the edge EMF over each face, i.e.
`∂_t B_x = −(∂_y E_z − ∂_z E_y)` etc., implemented via the generalized-Stokes helper
(`stokes.hpp` pattern: Face variable, Edge flux, `IsSubmanifold` sign handling). Because Bf
is `WithFluxes` with edge fluxes, the existing `FluxDivergence`-analog for Face variables (the
`Stokes` update in the example) performs the curl and the same weighted-sum stage combination
used for `cons`. Validate the face-curl signs against the untested-TODO (§0) with the divergence
pulse and field-loop tests before trusting it.

### 3.4 Face → cell-center projection
After each Bf update and before the next Riemann solve / `FillDerived`, set
`IB1 = ½(Bf_x|_{i} + Bf_x|_{i+1})`, etc. (the same averaging the example uses to form `C_cc`,
`advection_package.cpp:252`). This single projection is what keeps the **entire rest of the
code unchanged**: EOS, fast-speed/`c_h`, HLLD wave speeds, non-ideal η evaluation
(`PrecomputeNonidealEta` reads `prim` B), radiation, self-gravity coupling, history/`.hst`
diagnostics, and the analysis notebooks all continue to read cell-centered B. This projection
discipline is the primary mechanism for containing the blast radius (§5).

### 3.5 Interaction with the non-ideal EMFs (Ohm / AD / Hall)
Today `resistivity/ambipolar/hall.cpp` add their EMFs to `cons.flux(dir, IBn)` at cell faces
(§1.4). On the CT path they must instead deposit onto the **shared edge EMF** used by the curl.
Two clean routes:
- Compute the non-ideal edge EMFs directly at edges (currents `J = ∇×B` are naturally
  edge-centered — the resistivity kernel already forms `j2,j3` from face-averaged B), and add
  to the edge EMF array. This is *more* natural than the current face-flux placement and
  should reduce averaging error.
- Or keep the existing face-flux kernels and fold their transverse-B face fluxes into the same
  arithmetic/upwind average as the ideal flux (minimal code change, slightly more averaging).
Constraints to preserve: **AD** must remain RKL2-STS-compatible — the STS registers
(`diffusion.cpp` `RKL2Step*`) currently super-time-step the parabolic induction flux; under CT
the STS must advance `Bf` via the same curl, so the STS `FluxDivergence` calls become
Stokes-curl calls on the Bf edge flux. **Hall** stays unsplit and dispersive (needs the Ohmic
floor); its EMF simply adds to the ideal edge EMF in the unsplit stage. The energy (Poynting)
flux stays cell-face on `IEN` and is unchanged.

### 3.6 Keeping GLM selectable
Add a `<hydro> divergence_control = glm | ct` (or reuse the fluid/solver selection). On `glm`:
current code path verbatim (ψ, Dedner source, cell-centered induction flux, cell prolongation).
On `ct`: register `Bf`, skip ψ/Dedner, route induction through the edge-EMF/curl path, run the
projection. The switch is set at package init (`hydro.cpp`) so the task list
(`hydro_driver.cpp`) branches once. This satisfies the Phase-2 gate requirement that CT-vs-GLM
be run on the *same* fiducial collapse.

### 3.7 Task-list changes (`hydro_driver.cpp`)
Per stage on the CT path, mirroring `example/fine_advection/advection_driver.cpp`:
1. reconstruct + directional HLLD (unchanged) → cell-face fluxes.
2. **assemble edge EMF** from face fluxes (+ non-ideal edge EMFs).
3. `AddFluxCorrectionTasks` on Bf (edge-EMF correction across C–F) — replaces/augments the
   current cons flux correction for the induction components.
4. **Stokes curl** update of Bf + weighted-sum stage combine.
5. cell-centered conserved update (`FluxDivergence` for hydro vars) as today, **excluding**
   the induction components.
6. `Bf → CC B` projection.
7. boundary exchange (Bf FillGhost + cons), `FillDerived` (cons→prim, now consuming projected B).

---

## 4. Validation suite

Each is a restartable, machine-readable-norm test (fits the Phase-0 test-configuration suite).
Pre-register tolerances before looking at the collapse.

1. **∇·B round-off check (every test):** assert max|∇·B|·Δx/|B| ≤ ~1e-12 (round-off) on CT vs
   O(1) truncation growth on GLM. The `divC/divD` diagnostic in the example
   (`advection_package.cpp:262`) is the template.
2. **Field-loop advection (Gardiner & Stone 2005 §):** advect a weak B loop diagonally; measure
   (a) ∇·B ≈ round-off, (b) magnetic energy decay rate, (c) loop-shape distortion. This is the
   decisive discriminator between arithmetic (Stage A) and upwind (Stage B) EMF averaging, and
   exercises the Bf ghost exchange and projection. Run **with AMR** (a refined patch over the
   loop) to exercise Tóth–Roe prolongation + edge-EMF correction — the highest-risk operators.
3. **Orszag–Tang vortex:** 2D nonlinear MHD; compare CT vs GLM vs Athena++ reference at fixed
   resolution (density/pressure structure, energy history). Confirms the CT update is not just
   div-clean but *accurate*.
4. **Divergence pulse (audit "GLM pulse" analog):** initialize a localized ∇·B ≠ 0 seed. GLM
   advects/damps it; CT must keep it at round-off from t=0 (CT cannot create ∇·B, and with a
   div-free IC it never appears). Also the direct check of the face-curl signs (§3.3, §0).
5. **Non-ideal CT regression:** damped Alfvén eigenmode (AD, 0.04% today) and Hall whistler
   (0.4% today) rerun on the CT path — proves the edge-EMF re-routing preserves the validated
   non-ideal physics and STS compatibility.
6. **Decisive flux-retention comparison (the Phase-2 gate):** run the fiducial magnetized BE
   collapse (`runs/prod_v8` config, reduced resolution first) with `divergence_control=glm` and
   `=ct`, matched in every other input. Compare the flux-retention curve / retained-flux number
   at matched central density decades. **Acceptance:** either the two agree within the
   pre-registered numerical tolerance (⇒ divergence control does not set the result, claim is
   defensible), or the difference is quantified and reported as a systematic uncertainty. This
   is the whole point of Phase 2.

---

## 5. Risk / effort / blast radius

### 5.1 Blast radius — files touched
- **New:** CT EMF assembly + curl (new files under `src/hydro/`, e.g. `ct/`), analogous to the
  glmmhd/ and rsolvers/ layout.
- **Modified core:** `src/hydro/hydro.cpp` (register `Bf`, GLM/CT switch, projection in
  FillDerived), `src/hydro/hydro_driver.cpp` (task-list branch), `src/main.hpp` (possibly a
  `MagneticDiscretization` enum; ψ slot gating).
- **Modified induction physics:** `src/hydro/diffusion/{resistivity,ambipolar,hall,diffusion}.cpp`
  (edge-EMF deposition + STS-curl), `src/hydro/glmmhd/dedner_source.cpp` (GLM-only gate).
- **AMR:** register Tóth–Roe on Bf (one line) — but the operator itself is framework-provided.
  `src/hydro/prolongation/custom_ops.hpp` unchanged (it already handles Face TEs if reused).
- **IC:** `src/pgen/collapse_be.cpp` (face-init of Bf; trivial for uniform B0z, general case
  needs a vector-potential init).
- **I/O / restart:** Bf added to restart; **restart format changes** ⇒ old GLM restarts are not
  CT-resumable (and vice-versa) — a fresh `t=0` CT baseline is required, not a mid-run swap.
- **Analysis:** notebooks unchanged *because of the projection* (they read CC B / `prim`).

### 5.2 Why high-risk
- It rewrites the induction half of the **core update** and changes the primary magnetic
  variable — a bug here silently corrupts every downstream flux number, i.e. the headline result.
- It touches **AMR operators** (div-free prolongation + edge-EMF correction); errors there
  appear only at coarse–fine boundaries under refinement, exactly where the collapse spends its
  resolution.
- **Restart-format change** breaks resume-compatibility with the live GLM production chain —
  must be a clean fresh-start branch, never an in-place swap of `runs/prod_v8`.
- The framework face-curl path is **not yet MHD-tested** (§0 TODO) — we are its first real user.
- GPU/Kokkos: the edge-EMF assembly and curl are new kernels with staggered index arithmetic;
  correctness must be established on-device (Hopper), not just on the CPU build.

### 5.3 Effort estimate
Roughly **6–10 increments**, each with its own gate:
(1) register Bf + projection + IC, GLM/CT switch, arithmetic EMF, ideal-only, single level →
field-loop + OT + div-pulse green. (2) AMR: Tóth–Roe + edge-EMF correction → field-loop-AMR
green. (3) upwind (GS05) EMF → field-loop dissipation gate. (4) non-ideal edge-EMF re-routing +
STS-curl → AD/Hall regressions green. (5) restart/output plumbing + provenance. (6) GPU
validation. (7) the CT-vs-GLM collapse gate. The framework support is what makes this *weeks of
increments* rather than *a multi-month staggered-AMR port*.

### 5.4 Dependencies
- **Builds on Phase 1** (`PhysicalUnits`): CT changes nothing about units, but the B-unit and
  the cross-code √(4π) convention must be single-sourced *before* the CT-vs-GLM and CT-vs-Athena++
  flux comparisons, or a unit discrepancy will be misread as a divergence-control systematic.
- **Must precede** Phase 3 conductivity work's flux-loss claims (a flux-loss number computed on
  a GLM ∇·B-transporting background is not trustworthy).
- **One-variable-per-experiment:** the CT-vs-GLM comparison must change *only* the divergence
  control — same resolution, EOS, non-ideal coefficients, radiation, gravity, IC seed.

---

## 6. Bottom line

Parthenon already carries the expensive, high-risk part — a face/edge staggered data model with
divergence-free AMR prolongation (Tóth–Roe), edge-EMF coarse–fine correction, staggered ghost
exchange, and a generalized-Stokes curl update, with a working `do_CT_advection` reference.
AthenaPK's job is the MHD-specific glue: a face-B field, an EMF built from the HLLD face fluxes it
*already* computes (arithmetic first, Gardiner–Stone upwind for production), a curl update, a
face→cell-center projection that keeps every existing B consumer untouched, and re-routing of the
already-face-based non-ideal EMFs onto shared edges. Recommended flavor: Gardiner & Stone upwind
CT, matching the Athena++ reference for the cross-code gate. High risk, wide coupling surface,
restart-format change — but staged into ~6–10 independently-gated increments, ending in the
decisive CT-vs-GLM flux-retention comparison on the fiducial collapse.
