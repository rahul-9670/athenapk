# AthenaPK Physics-Completion Plan — Execution Document

**Audience:** a future Claude (Opus/Fable) session executing one workstream at a time.
**Author basis:** all file paths, conventions, and existing-code facts below were
verified against the tree on 2026-07-07. Re-verify anything load-bearing before editing
(rule 1 of the working-style rules in `CLAUDE.md`).

This plan closes the five gaps identified in `PHYSICS_IMPLEMENTATION_REPORT.md` §"What's
genuinely missing":

| WS | gap | port source | size (est.) |
|----|-----|-------------|-------------|
| WS-1 | sink particles | `artemis/src/nbody/` + Parthenon swarms | L (4–6 increments) |
| WS-2 | chemistry thermal feedback | greenfield (gow17 paper) | M (3 increments) |
| WS-3 | radiation upgrades (means split, 2nd order, multigroup) | in-tree | M (3 increments, last optional) |
| WS-4 | dust evolution | minimal greenfield; `artemis/src/dust/` for the full option | M (2–3 increments) |
| WS-5 | numerical fidelity (multipole BCs, Hall cross-check, divB audit) | in-tree | S+S+S |

**Recommended execution order: WS-5a → WS-1 → WS-3a,b → WS-2 → WS-4 → WS-5b,c.**
Rationale: WS-5a (multipole gravity BCs) is small, independently testable, and shrinks
the box from L=52 to ~L=16 — every subsequent long run gets ~30× cheaper per level of
the base grid. WS-1 is the science unlock. WS-3a (Planck/Rosseland split) must precede
WS-2 (gas–dust coupling needs a trustworthy T_dust). WS-4 feeds back into opacity and
ionization, so it goes after both consumers exist. WS-5b/c are validation work that can
interleave anywhere.

---

## 0. Global ground rules (apply to every workstream)

These encode how the previous ports (self-gravity, RT, chemistry) were done successfully.

1. **Every new feature is behind an input flag, default OFF, and the OFF state must be
   bit-identical to current production.** Acceptance test for *every* increment: run
   `runs/prod_t3_ad_ohm`-style small config (or a 64³ shrink of it) with the feature
   flag off, before/after the code change, and diff the final `.phdf` (`h5diff` or the
   analysis venv). The matched-physics comparison with Athena++ must never silently
   shift. New physics is AthenaPK-only divergence — justified and documented per rule 6
   (Athena++ has no counterpart; the paired-artifact rule applies to the *matched
   tiers*, which keep flags off).
2. **Increment discipline** (rule 8): each increment = design note appended to
   `DEV_LOG.md` → implement → CPU build (`make -C build_cpu athenaPK -j`, front-end OK)
   → run the increment's validation on CPU → record the numbers in `DEV_LOG.md` → GPU
   build via `sbatch runs/submit_build_gpu.sh` → re-run one validation case on GPU and
   confirm agreement to round-off-accumulation level. Never stack two increments
   without a green gate between them.
3. **Units:** everything internal to kernels is code units (4πG=1, HL magnetic
   convention, p = ρT_code, T_unit = 10.015 K). Microphysics evaluated in cgs must
   convert on the way in and out, following the `ionization.hpp` pattern (explicit
   `*_unit` members, defaults = FHC calibration, overridable by input keys). State the
   unit of every new quantity at its point of use (rule 4).
4. **Task-list integration:** new per-step physics goes into `hydro_driver.cpp`
   following the existing pattern — self-gravity (`:583-590`), radiation (`:600`),
   chemistry (`:610`). Operator-split fields need the `OperatorSplit` metadata user
   flag (see `radiation.cpp:151-157`) so the hydro integrator's flux-divergence sweep
   leaves them alone.
5. **GPU correctness:** all device code is `KOKKOS_INLINE_FUNCTION` on POD structs
   captured by value (no host pointers in kernels) — the `IonizationModel` /
   `OpacityParams` / `EosTable` patterns are the templates to copy.
6. **Restarts:** any new evolved state must survive `parthenon.out2.*.rhdf` restart.
   Parthenon handles `Metadata::Independent` fields automatically; **swarms need
   explicit verification** (WS-1 increment 1 checks this first).
7. **Validation is quantitative** (rule 5): each increment below names its test and its
   numeric acceptance threshold. "Looks right" is not a gate.

---

## WS-1: Sink particles

### Goal
Continue the simulation past second-core formation: replace the unresolvable
protostellar interior with an accreting, moving point mass, so the run proceeds from
core formation (~1 t_ff) into the main accretion phase (~10⁴–10⁵ yr). Non-goal:
protostellar feedback (outflows, radiation from the sink photosphere) — leave hooks,
don't implement.

### Physics and equations
A sink of mass M, position **x**ₛ, velocity **v**ₛ, spin **L**ₛ obeys

```
dxₛ/dt = vₛ
dvₛ/dt = −∇Φ_gas(xₛ)  +  Σ_other sinks −G M_j (xₛ−x_j)/(|xₛ−x_j|² + ε²)^{3/2}  +  accretion recoil
dM/dt  = Σ_{cells in r_acc, bound} ΔM_cell/dt        (conservative extraction from the gas)
```

Gas feels the sink via direct softened gravity added as a momentum/energy source
(Artemis `nbody_gravity.hpp` pattern), **not** via the Poisson RHS — this avoids the
self-force problem and keeps the MG solver untouched. In code units G = 1/(4π).

**Creation criteria** (Federrath et al. 2010, all must hold in a candidate cell):
1. ρ > ρ_sink, with ρ_sink set by the Truelove limit at the finest level:
   λ_J = 2π c_s/√ρ (code units, cf. `refinement/jeans.cpp:23`) resolved by n_J cells ⇒
   `ρ_sink = 4π² c_s²(ρ_sink) / (n_J · dx_min)²` — iterate once with the EOS-table c_s.
2. Cell is on the finest AMR level and is a local minimum of Φ within r_acc.
3. ∇·v < 0 (converging flow).
4. Total energy (thermal + kinetic + magnetic − gravitational binding) of the control
   volume < 0.
5. No existing sink within 2 r_acc (else that sink accretes instead).

**Accretion:** cells within r_acc = `sink_racc_cells` (default 4) × dx_min that are
bound to the sink and inside the Bondi radius. Port the Artemis kernel
(`artemis/src/nbody/particle_base.hpp:193-240`): quadratic ramp in radius, extract
ΔM = f·(ρ − ρ_sink/3)·V_cell capped so the cell keeps ≥ ρ_sink/3; momentum extracted
proportionally; the *angular* momentum of extracted gas about the sink accumulates in
**L**ₛ (bookkeeping only). **Magnetic field is not accreted** — B stays on the grid
(standard practice; physically defensible here because the non-ideal terms have already
decoupled the flux at these densities; state this in the paper).

### Implementation design
- New package `src/sinks/` registering a Parthenon **swarm** (vendored Parthenon has
  full support: `external/parthenon/src/interface/swarm*.cpp`, including MPI migration
  in `swarm_comms.cpp`). Real fields: mass, x/y/z, vx/vy/vz, Lx/Ly/Lz, t_created;
  integer: id.
- Sink advance: KDK leapfrog, subcycled inside the hydro dt if
  |vₛ|dt > 0.5 dx_min. Gas Φ interpolated to xₛ by CIC from `grav.phi` (ghosts are
  filled — the solver already exchanges them).
- Gas source kernel: for each cell within `sink_gravity_rmax` (default: whole block set
  intersecting 32 r_acc; beyond that the sink's monopole is *added to the multipole
  boundary moments of WS-5a* — if WS-5a is not yet merged, apply direct force
  everywhere; N_sink is ≤ a few, cost is one fused kernel).
- Creation + accretion run once per step after the hydro update and Poisson solve,
  before output; both are global operations (candidate reduction + `MPI_Allreduce`).
- AMR: tag blocks containing a sink (or within 2 r_acc) `AmrTag::refine` to pin the
  finest level around it (add to `refinement/jeans.cpp` dispatch or a new criterion
  combinator).
- dt: add `min(dt, cfl · dx_min/|vₛ|max)` to the global constraint.

### Inputs (`<sinks>`)
`enabled` (false), `n_jeans_sink` (default = `refinement/njeans`), `racc_cells` (4),
`soft_cells` (2), `rho_sink_code` (−1 ⇒ auto from Truelove formula), `merge_dist_racc`
(2.0), `max_sinks` (64), `subcycle_cfl` (0.5).

### Increments and gates
1. **Swarm plumbing:** one inert particle advected through an AMR mesh across MPI
   ranks; survives restart. Gate: position error vs analytic ballistic path < 1e-12;
   restart bit-identical continuation.
2. **Static point-mass gravity on gas:** disable creation/accretion, place one sink,
   run a Bondi-flow setup (new small pgen or reuse `polytrope.cpp` box). Gate:
   steady-state Bondi accretion-rate diagnostic (measured flux through r = 8 dx) within
   10% of the analytic ṀB for γ=1.4 at the chosen resolution.
3. **Two-body dynamics:** two sinks, no gas forces. Gate: circular orbit energy drift
   < 1e-6 per orbit over 100 orbits (leapfrog symplectic check).
4. **Creation criteria:** collapse run (t3 config shrunk to numlevel=8, 128³): sink
   forms at ρ within 2× of `rho_sink`, at the density maximum, exactly once. Gate: no
   spurious second sink for ≥ 0.1 t0 after formation.
5. **Accretion + conservation:** continue increment 4. Gates: total (gas + sink) mass
   conserved to < 1e-10 relative per step; momentum likewise; accretion rate within a
   factor 2 of Shu (1977) 0.975 c_s³/G after the transient; **dt recovers** to within
   10× of the pre-first-core hyperbolic dt (this is the entire point — record the dt
   history).
6. **Full-physics run:** t4 config + sinks on 8×H100 (budget per the t4 memory note:
   `cpus-per-task=8`, `derefine_count=50`). Gate: runs ≥ 5000 yr past sink formation
   without dt collapse; mass budget closes.

### Risks
- Swarm ↔ AMR rebalancing interplay under heavy `derefine_count` churn — increment 1
  must test on an *adaptive* mesh, not uniform.
- GLM ψ and `dfloor` in strongly evacuated accretion cells: cap extraction (the ρ_sink/3
  floor) and watch `relDivB` in the hst output; acceptance: no growth vs no-sink run.
- The barotropic overwrite (non-RT runs) resets e_th in accretion-heated cells — fine
  (it does so everywhere), but document.

---

## WS-2: Chemistry thermal feedback (gas thermochemistry)

### Goal
Replace "barotropic or RT-only" gas temperature with actual heating/cooling where it
matters: the envelope (n = 10²–10⁷ cm⁻³). **Design constraint: only active when RT is
on** (RT owns e_th; the barotropic path would overwrite any Λ/Γ work — see
`collapse_be.cpp:580` and `hydro.cpp:566-570`).

### Physics (rates in erg cm⁻³ s⁻¹; from Gong, Ostriker & Wolfire 2017 and standard fits)
Heating:
```
Γ_CR  = q_CR · ζ · n_H            q_CR ≈ 20 eV per primary (density-dependent fit optional)
Γ_H2  = (0.2–4.48 eV branch) · kgr n_H² x_H   (formation heating; use 0.2 eV default — grain-borne)
Γ_PE  = 1.3e-24 n_H ε G_eff       with G_eff = G_0 e^{−2.5 A_V}; in the shielded core ≈ 0 (keep for the envelope)
```
Cooling:
```
Λ_C+   : 158 µm fine structure, 2-level atom with n_cr; needs x_C+ (already evolved!)
Λ_O    : 63 µm, fixed x_O = 3.2e-4 · depletion
Λ_CO   : rotational ladder — use the Whitworth & Jaffa (2018) or Goldsmith & Langer
         analytic fit Λ_CO(n, T, x_CO); needs x_CO (already evolved)
Λ_gd   = α_gd · n_H² · √T · (T − T_dust),   α_gd ≈ 3.2e-34 erg cm³ K^{-3/2}
         (gas–dust; DOMINANT for n ≳ 10^4.5 — this is what pins T ≈ T_dust in the core)
T_dust = (E_r/a_rad)^{1/4} from the M1 field  (code: (Er/arad_code)^{1/4} · T_unit)
```

### Implementation design
- Extend the chemistry package (`src/chemistry/`): a `ThermoParams` POD + device
  function `dedt(rho, T, x_i, T_dust) → Γ−Λ` in a new `thermo.hpp`, following the
  network-struct pattern. Integrate temperature *inside* the existing per-cell
  sub-cycler (`integrate_cell`) so stiff cooling and the abundance ODEs share the
  sub-step controller — add e (or T) as the 5th/6th integrated variable with the same
  semi-implicit P/L split (Λ terms that are ∝ T-dependent loss linearize the same way).
- The cell's IEN update happens once at the end of `ReactScalars` (subtract old e_int,
  add new), matching how MatterCoupling edits IEN. Ordering per step: RT matter
  coupling first (sets T toward radiative equilibrium with dust), thermochemistry
  second (adds the gas-phase Γ−Λ and gas–dust exchange on top).
- Λ_gd needs T_dust: pass `rad.Er` into the chemistry pack (extend the pack list in
  `chemistry.cpp:111`).

### Inputs (`<chemistry>`)
`thermo = true|false` (false), `q_cr_ev` (20), `alpha_gd` (3.2e-34), `x_O` (3.2e-4),
`G0` (1.7), `h2_heat_ev` (0.2), `dust_beta` … all with FHC defaults.

### Increments and gates
1. **Standalone equilibrium curve:** host-side unit test (pattern:
   `src/chemistry/tests/`) computing T_eq(n_H) for n = 10²–10⁸ cm⁻³ at equilibrium
   abundances. Gate: reproduces the canonical molecular-cloud curve — T_eq within
   ±30% of gow17 Fig. set values (~10–15 K at n=10⁴, converging to T_dust at n≳10⁵).
2. **Wired in, single-cell-in-a-box run:** uniform box, RT on, check relaxation to the
   increment-1 equilibrium. Gate: |T − T_eq|/T_eq < 5% after 10 cooling times; no dt
   collapse (the sub-cycler must absorb the stiffness — verify `chem_trunc_total`
   stays 0).
3. **Collapse A/B:** t2-scale run with `thermo=on` vs `off`. Gate: core (n>10⁷)
   temperature difference < 10% (dust coupling must dominate there — if not, the
   coupling constant is wrong); envelope shows the expected 10–15 K plateau.

### Risks
- Double-counting emission: Λ_gd moves energy gas→dust, and the dust re-emits via the
  *existing* RT dust opacity. Do **not** add a separate Λ_dust-continuum term — the M1
  field already is the dust radiation. Document the energy pathway explicitly.
- The gow17_reduced network lacks species some Λ fits want (no atomic O evolved): use
  fixed x_O, state the assumption.

---

## WS-3: Radiation upgrades

### 3a. Separate Planck and Rosseland means  (do first — WS-2 depends on it)
**Physics:** the matter coupling `S_E = ĉρκ(E−B)` should use the Planck mean κ_P
(emission-weighted); flux attenuation and the diffusion limit should use the Rosseland
mean κ_R. Bell & Lin is a Rosseland fit; using it for emission under-cools optically
thin gas where κ_P/κ_R can be ~2–5 for dust.

**Design:** extend `OpacityParams` with a second model channel; add `opacity_model =
table`: a binary (log ρ, log T) → (κ_P, κ_R) table loaded like `eos_table.bin`
(generator script `src/radiation/gen_opacity_table.py` from the Semenov et al. 2003
dust model, low-T range stitched to Bell & Lin above sublimation). `AbsorptionOpacity`
→ κ_P in `MatterCoupling`'s Newton loop; a new `RosselandOpacity` → the flux-attenuation
`a` (`radiation_moments.cpp:322`) — currently both call sites use the same function;
split them.
**Inputs:** `opacity_model=table`, `opacity_table_file`; keep `belllin` behavior
bit-identical when selected.
**Gate:** (i) equilibration test — uniform gas+radiation relaxes to T_eq with the
correct κ_P timescale (analytic, < 1%); (ii) t2-scale collapse: first-core central
entropy shifts by a *documented* amount vs `belllin` (expected: slightly warmer
envelope; record the number).

### 3b. Second-order transport
**Design:** replace the donor-cell states in `CalculateRadFluxes`
(`radiation_moments.cpp:95-116`) with PLM + minmod reconstruction of (E_r, F_r) — the
stencil grows to ±2, ghosts already 2+ for PLM hydro (verify `nghost` covers the rad
sub-cycle exchange). Flag `radiation/reconstruction = dc|plm` (default dc).
**Gate:** free-streaming pulse (`rad_pulse` pgen) — numerical diffusion of the front
reduced by ≥ 2× at 128 cells (measure FWHM growth); crossing-beam/shadow qualitative
check documented (M1 merging artifact is expected and *not* a failure).

### 3c. Multigroup (OPTIONAL — decide after 3a/3b science need)
N_g groups, each with (E_g, **F**_g), group Planck emission fractions from
precomputed ∫B_ν tables, per-group κ from the 3a table generator extended in ν.
Memory ×N_g on 4 fields; only justified if the FHC SED or the sublimation-gap
transient becomes a paper target. Increment plan deferred; do not start without an
explicit science requirement.

### 3d. creduc convergence study (no code)
Run matrix creduc ∈ {300, 1000, 3000} on a t2-scale config; gate: first-core properties
(M_fc, R_fc, T_c at fixed central density) vary < 5% between 1000 and 300 — else
production creduc must drop.

---

## WS-4: Dust evolution

### Recommended scope: single-moment grain model (NOT the full Artemis dust fluid)
Justification: at core densities the grains that matter (< µm) have Stokes numbers
≪ 1 — they co-move with the gas; drift is negligible. What *does* change the physics is
(i) grain **growth** raising a_char (lowers total cross-section → lowers κ and raises
η's), and (ii) **sublimation** (already switched at 1.5e3 K). The full Artemis
two-fluid dust + coagulation port (`artemis/src/dust/`, `dust/coagulation/`) is the
fallback if drift ever becomes a target — do not start there.

### Physics
Evolve two passive scalars (extend `nscalars` by 2): dust-to-gas ratio `f_dg` and
characteristic grain radius `a_c` (mass-weighted). Sources, operator-split like
chemistry:

```
da_c/dt = (f_dg ρ / ρ_grain) · v_rel / 4        (monodisperse sweep-up; Ormel-style)
v_rel   = √(v_Brownian² + v_turb²),  v_turb ≈ √(α_t) c_s · Re^{1/4} St   (standard closure)
f_dg    : advected; → 0 smoothly above T_subl (sublimation), restored below (no re-formation lag, v1)
```

Consumers:
- **Opacity:** κ_dust scale ∝ (f_dg/0.01) · (a_ref/a_c) in the geometric limit
  (κ ∝ total area ∝ 1/a at fixed mass); wire as multiplicative factors into
  `AbsorptionOpacity`/`RosselandOpacity`.
- **Ionization:** `IonizationModel` currently fills MRN bins host-side once
  (`SetupGrainBins`, `hydro.cpp:283-305`). Make f_dg a per-cell multiplier on `ng[k]`
  (one extra argument through `Diffusivities`; the bins' *shape* stays MRN with
  a_min/a_max scaled by a_c/a_ref — geometric rescale, no re-tabulation needed).

### Inputs (`<dust>`)
`evolve = false`, `alpha_turb` (1e-2), `a_ref_cm` (1e-5), `sublimation = true`, and the
scalar indices (auto-assigned after the chemistry block).

### Increments and gates
1. Single-cell growth test vs the analytic monodisperse solution (< 5% over 3 decades
   in a_c); sublimation switch round-trips without mass error.
2. Consumers wired: with `evolve=on` but growth rates zeroed, results bit-identical to
   `off` (f_dg=0.01, a_c=a_ref everywhere). Then growth on: t2-scale collapse — record
   Δκ(ρ) and Δη_A(ρ) profiles; gate is *documented sensitivity*, not a target value.
3. Full run + report figure: η_A(ρ) with static vs evolved grains.

### Risk
a_c feeding ionization couples dust to the dt estimator — cap a_c growth per step
(≤ 10%/step) to avoid dt thrash; the eta cache (`PrecomputeNonidealEta`) already
amortizes the cost.

---

## WS-5: Numerical fidelity

### 5a. Multipole gravity boundary conditions  (small, do FIRST overall)
**Physics:** replace Φ=0 Dirichlet walls with the exterior multipole expansion of the
interior mass, Cartesian moments through quadrupole about the center of mass:

```
Φ(r) = −(4πG/4π) [ M/|d| + (p·d)/|d|³ + ½ d·Q·d/|d|⁵ ],   d = r − r_com
M = Σρ dV,  p = Σρ(x−r_com)dV (≡0 about the COM),  Q_ij = Σρ(3x_i x_j − δ_ij x²)dV
```
(code units: the prefactor is `four_pi_G/(4π)`; with four_pi_G=1 that is 1/(4π)).

**Design:** compute {M, r_com, Q} in `FillPoissonRHS` (`self_gravity.cpp:214`) — it
already does a global reduction for the swindle; extend it (9 more scalars, one
Allreduce). Store in Mutable package params. New BC option `"multipole"` in the BC
enrollment (`self_gravity.cpp:106-133`) and in the packed `SetBoundary` fast path
(`poisson_equation.hpp:208`, add `type==3`): ghost value = expansion evaluated at the
ghost-cell center (FixedFace semantics: set the *ghost* so the face value equals the
expansion — mirror formula `ghost = 2Φ_mp(face) − interior`). WS-1 sinks add their
point masses to M/p/Q.
**Gates:** (i) uniform-density sphere in a box: Φ error vs analytic < 1% everywhere
(vs ~10–30% for zero-Dirichlet at L=2R); (ii) **the known failure becomes the test**:
the L=16 matched-box collapse (previously killed by BC clipping, cf.
`athenapk-mhd-collapse-bug` memory) must now track the L=52 run's central-density
evolution to < 5% up to 10⁴× rhocrit; (iii) BiCGSTAB iteration count does not grow.

### 5b. Hall term cross-validation (no counterpart in Athena++)
The tree already has `src/pgen/cshock.cpp` — verify what it initializes; extend/use it
for a **Hall C-shock**: steady oblique shock with AD+Hall admits a semi-analytic ODE
solution (Wardle-type). Deliverables: (i) a Python reference integrator in
`runs/validation/` (analysis venv); (ii) AthenaPK 1D runs at 3 resolutions.
**Gate:** downstream B rotation angle and shock thickness converge to the ODE solution
at 2nd order; error < 2% at the highest resolution. This plus the existing whistler
dispersion (0.4%) closes the "no cross-code check" hole.

### 5c. div B fidelity audit (explicitly NOT a CT rewrite)
Full constrained transport is out of scope: it replaces the entire cell-centered
GLM infrastructure (fluxes, prolongation, restarts) — a new code, not an increment,
and it would break bit-compatibility with all matched runs. Instead, *bound* the
cleaning error: (i) add a per-AMR-level max|∇·B|dx/|B| history output next to the
existing `relDivB` (`hydro.cpp` hst block, `:492-497`); (ii) sensitivity matrix over
`glmmhd_alpha ∈ {0.05, 0.1, 0.4}` and `dedner_plain` vs `extended` on a t2-scale run;
(iii) report section quantifying max relDivB at first-core formation.
**Gate:** documented numbers; if max relDivB at the core exceeds ~1e-2, escalate to a
design discussion (possible ψ-wave BC or higher c_h) — do not silently proceed.

---

## Execution mechanics for the driving session

- One workstream at a time; one increment per working session unless gates are trivially
  green. Append every design decision and every measured gate number to `DEV_LOG.md`
  (grep it first — several of these topics have history).
- Long runs go through SLURM (`sbatch`); never MPI on the front-end. GPU builds:
  `sbatch /beegfs/u/bbg6470/athenapk/runs/submit_build_gpu.sh`.
- Post-processing with `/beegfs/u/bbg6470/venvs/analysis_env/bin/python`; reuse
  `runs/fhc_ext.py` for profiles.
- Before each increment: re-read the touched files (they may have moved since this plan
  was written); if this plan contradicts the tree, **the tree wins and the mismatch is
  a finding** to log.
- Paper hooks: WS-1 gate 6, WS-3a gate (ii), WS-4 increment 3, and WS-5c produce
  numbers/figures that belong in `~/paper1_fossil_field/ms.tex` PENDING boxes — flag
  them when they land, don't edit the manuscript unprompted.
