# AthenaPK FHC Collapse — Physics & Implementation Report

**Scope.** This document describes every physical process implemented in this workspace's
AthenaPK fork (`/beegfs/u/bbg6470/athenapk`), the governing equations each module solves,
the discretization and code path that implements it, the runtime inputs that control it,
and the values used in production. All equations, constants, and defaults below were read
directly from the source files cited; production values come from
`runs/prod_t4_full/fhc.in` (the most complete "tier 4" configuration).

The target problem is the gravitational collapse of a magnetized, rotating, mildly
turbulent Bonnor–Ebert (BE) sphere through first-hydrostatic-core (FHC) and second-core
formation, with matched physics to a reference Athena++ setup.

---

## 1. The composite system being solved

Per timestep the code advances (in the shared code units of §2, Heaviside–Lorentz EM
convention, so no 4π appears in the MHD equations):

**Ideal GLM-MHD** (cell-centered, `src/hydro/glmmhd/`, Riemann fluxes in `src/hydro/rsolvers/`):

```
∂ρ/∂t   + ∇·(ρv)                                   = 0
∂(ρv)/∂t + ∇·(ρvv + (p + B²/2)I − BB)              = −ρ∇Φ            (gravity source)
∂E/∂t   + ∇·[(E + p + B²/2)v − B(v·B)]             = −(∇Φ)·(mass flux) + non-ideal Poynting + RT coupling
∂B/∂t   + ∇·(vB − Bv) + ∇ψ                         = −∇×E_nonideal
∂ψ/∂t   + c_h² ∇·B                                 = −(c_h/c_p)² ψ    (Dedner GLM cleaning)
```

with total energy `E = e_int + ½ρv² + ½B²` and, for the ideal-gas path, `p = (γ−1) e_int`.
Magnetic pressure is `B²/2` and `v_A = B/√ρ` (Heaviside–Lorentz — the #1 pitfall when
comparing to Athena++'s Gaussian convention).

**Poisson equation** for self-gravity (`src/self_gravity/`):

```
∇²Φ = 4πG (ρ − ρ_mean·swindle)          with 4πG = 1 in code units
```

**Non-ideal induction EMFs** (`src/hydro/diffusion/`), added as fluxes to B and E:

```
E_Ohmic     = η_O J
E_ambipolar = η_A (J − (J·b̂) b̂)        (current ⊥ B)
E_Hall      = η_H (J × B)/|B|  [+ η_floor J stabilizer]
J = ∇×B  (no 4π, HL units)
```

**M1 gray radiation transport** (`src/radiation/`), operator-split with reduced speed of
light ĉ = c/creduc:

```
∂E_r/∂t + ∇·F_r          = ĉ ρκ_a (a_R T⁴ − E_r)
∂F_r/∂t + ∇·(ĉ² P_r)     = −ĉ ρ(κ_a+κ_s) F_r
P_r = D(f) E_r  (M1 closure),   gas: ∂e/∂t = −c ρκ_a (a_R T⁴ − E_r), momentum kick from ΔF_r
```

**Chemical kinetics** (`src/chemistry/`), operator-split per cell on passive scalars
x_i = n_i/n_H (H2 network or reduced GOW17 H–C–O network; §8).

**Thermodynamic closure** — one of:
- barotropic pressure–density relation enforced every step (§4.1), or
- tabulated multi-Saha protostellar EOS (`eos=hydrogen`, §4.2) with RT owning the gas
  thermal energy.

---

## 2. Code units and normalization

Set in `src/pgen/collapse_be.cpp` from three physical inputs (`mass`, `temperature`, `f`).
The normalization (Tomida-style) fixes **4πG = 1**, **c_s(10 K) = 1**, **BE central
density = 1**:

```
v0   = 1.9e4 · √(T/10 K)  cm/s                      (isothermal sound speed)
m0   = M_total / (197.561 · f)                       (code mass unit; 197.561 = critical-BE mass)
rho0 = v0⁶ / (64 π³ G³ m0²)
t0   = 1/√(4πG rho0),      l0 = v0 · t0
B_unit = √(rho0) · v0                                (HL; written to units.json)
```

For the production IC (mass=6 M☉, T=10 K, f=5), `runs/prod_t4_full/units.json` records:

| quantity | value (cgs) |
|---|---|
| length l0 | 2.80630e16 cm (~1876 au) |
| time t0 | 1.47700e12 s (~46.8 kyr) |
| mass m0 | 1.20819e31 g |
| density rho0 | 5.46683e-19 g/cm³ |
| velocity v0 | 1.9e4 cm/s |
| B unit | 1.40482e-5 G |

Temperature: `p = ρ·T_code` with `T_unit = μ m_H v0²/k_B = 10.015 K` (μ=2.29, the value
*implied* by c_s=1.9e4 cm/s at 10 K). Note the deliberate two-μ convention: μ=2.29 for the
temperature calibration, μ_n=2.33 for number densities n = ρ/(μ_n m_H). Do not unify them
(comments in `radiation.cpp:61-64` and `ionization.hpp:74-77`).

---

## 3. Problem generator: Bonnor–Ebert sphere (`src/pgen/collapse_be.cpp`)

**Physics.** A critical BE sphere is the marginally stable isothermal self-gravitating
equilibrium; multiplying its density by an enhancement factor f > 1 makes it unstable and
it collapses on ~t_ff. The approximate BE profile (Tomida 2011 thesis) is used:

```
ρ(r) = f · (1 + r²/r_c²)^(−3/2),   r_c² = 26/3,   sphere radius r_BE = 6.45 (code)
t_ff = π √(3/(8f))  (code)
```

**Initial conditions per cell** (lines 194–256):
- Inside r < 6.45: BE profile density, optional m=2 azimuthal perturbation
  `ρ → ρ(1 + amp·(r/r_c)²cos 2φ)`, solid-body rotation Ω set by `omegatff = Ω·t_ff`.
- Outside: uniform ambient at ρ(r_BE) (pressure-continuous), at rest.
- Pressure p = ρ (isothermal, c_s = 1); the energy seed uses e_int = p/(γ−1), except for
  the tabulated EOS where cold H2 is monatomic-like so e_int = 1.5p.
- Uniform vertical field B_z = `B0z` (0.15 code = 7.47 µG recalled; mass-to-flux
  μ_Φ = 31.9 in the paper's convention — those two conversions are from memory,
  everything else here is from the file).
- GLM scalar ψ = 0.

**Initial turbulence** (lines 264–469, gated by `turb_mach > 0`): a deterministic random
Fourier mode sum, `nmodes` modes with |k| uniform in [k_min,k_max]·2π/L_box, amplitude
∝ k^(−α/2), mixture ζ of solenoidal vs compressive unit polarizations. The field is scaled
to RMS Mach `turb_mach` using the *analytic* ensemble RMS (√(½Σ t_m²)) and the *analytic*
region-mean velocity is subtracted (sphere window function W(x)=3(sin x − x cos x)/x³) so
the cloud gets no net momentum kick — both closed-form so every rank applies identical
factors (a mesh reduction is not per-block consistent here; this fixed the 0.34 c_s bulk
kick bug of the old IC).

**Radiation seeding** (lines 479–512): when RT is active, E_r = a_rad T⁴ (LTE with the
gas), F_r = 0, preventing a spurious t=0 cooling transient.

**Chemistry seeding**: with ≥5 scalars (gow17_reduced) a molecular-cloud state
x_H2 = 0.5, x_H+ = 1e-8, x_C+ = 1.6e-4 = x_e, x_CO ≈ 0; with 2 scalars (H2 network)
all-atomic hydrogen.

**Inputs** (`<problem/collapse_be>`): `mass` [M☉], `temperature` [K], `f` (enhancement,
>1 to collapse; production 5), `rhocrit` [g/cm³] (barotropic transition; production
1e-13), `amp` (m=2 amplitude, 0–~0.1; production 0), `omegatff` (production 0.02), `B0z`
(code B; production 0.15), and the turbulence block `turb_mach` (0.5), `turb_kmin/kmax`
(1/8), `turb_nmodes` (128), `turb_zeta` (0.5 = half solenoidal), `turb_alpha` (3.667,
i.e. Burgers-like E(k)∝k^-1.667·k²→ amplitude exponent), `turb_seed` (42), `turb_region`
("sphere" or "box").

---

## 4. Thermodynamics

### 4.1 Barotropic closure (default, no-RT runs)

`collapse_be::ApplyBarotropicCooling` (pgen lines 525–595), registered as an unsplit
source each stage. It **overwrites** the thermal energy every step (EOS enforcement, not a
cooling ODE), Masunaga–Inutsuka/Tomida form:

```
e_th = ρ/(γ−1) · √(1 + (ρ/ρ_crit)^(2(γ−1)))
```

→ isothermal (p = ρ c_s², T = 10 K) for ρ ≪ ρ_crit; adiabatic p ∝ ρ^γ for ρ ≫ ρ_crit.
This is the standard stand-in for radiative trapping at first-core formation
(ρ_crit = 1e-13 g/cm³ ↔ 182.9 code, γ = 1.4). The same kernel zeroes momentum outside
r_BE (fixed-velocity BC on the ambient), with a kinetic-energy correction when RT owns
the energy so the killed momentum isn't silently converted to heat.

### 4.2 Tabulated protostellar EOS (`eos = hydrogen`, second-core physics)

`src/eos/gen_eos_table.py` generates offline, `src/eos/eos_table.hpp` interpolates on
device. Full multi-Saha equilibrium at (ρ,T) over species {H2, HI, HII, e, HeI, HeII,
HeIII} (X = 0.716, Y = 0.284):

```
H2 ⇌ 2HI    (χ_d = 4.4781 eV, with nuclear-spin-weighted rot. partition function,
             θ_rot = 85.4 K, θ_vib = 5987 K)
HI ⇌ HII+e  (χ_H = 13.598 eV);  HeI ⇌ HeII+e (24.587 eV);  HeII ⇌ HeIII+e (54.418 eV)
e_int = 3/2 kT n_tot + n_H2 e_rotvib(T) + dissociation + ionization energies
p = n_tot kT;   c_s² = (∂p/∂ρ)_s by finite difference along the adiabat
```

This reproduces γ_eff ≈ 7/5 warm molecular, dips < 4/3 at H2 dissociation (~2000 K) and
H ionization (~1e4 K), → 5/3 ionized — the softening that triggers second collapse.

Runtime tables: forward grid (log ρ, log e_sp) → {P, c_s², log T[K]} and RT grid
(log ρ, log T) → e_sp; bilinear interpolation, pressure inverse by 30-step bisection.
Constraints enforced at parse time (`hydro.cpp:558-570`): requires `riemann=hlld` (HLLE's
Roe average bakes in constant γ) and `<physics> radiation=true` (else the barotropic
source would overwrite table energies). Inputs: `hydro/eos = adiabatic|hydrogen`,
`eos_table_file` (default `src/eos/eos_table.bin`), `gamma` (1.4; still used for c_h and
fallbacks).

**Floors** (`<hydro>`): `dfloor` (t4: 1e-6), `pfloor` (1e-8), optional `Tfloor/vceil/Tceil`.

### 4.3 Cooling: what "cooling" means in this code (three distinct channels)

There is **no optically-thin radiative-cooling ODE in the collapse problem**. The word
"cooling" appears in three unrelated places:

1. **Barotropic pseudo-cooling** (§4.1, non-RT runs). Physically: below ρ_crit the
   thermal relaxation time of dusty molecular gas is far shorter than the dynamical
   time, so the gas is pinned at 10 K — the code takes the infinite-cooling-rate limit
   by *overwriting* e_th with the barotropic value every step. It is an EOS
   enforcement, not a rate equation: compressive/shock heating below ρ_crit is
   discarded by construction, and above ρ_crit the √(1+(ρ/ρ_crit)^(2(γ−1))) form
   smoothly hands over to adiabatic p ∝ ρ^γ (radiative trapping mimic). Its only knobs
   are `rhocrit` and `gamma`.

2. **Radiative cooling/heating via the M1 matter coupling** (§8, RT runs). This is the
   *physical* cooling channel in t4: the gas exchanges energy with the radiation field
   at rate S_e = −c ρκ_a (a_R T⁴ − E_r) — emission cools, absorption heats — solved
   implicitly (Newton + Fleck factor) so the stiff optically-thin limit is stable. The
   envelope cools toward the E_r of the ambient field, the core heats when trapped.
   When RT is on, channel 1's overwrite is skipped entirely (`rad_owns_energy` in the
   pgen); only its outside-sphere momentum BC remains.

3. **`TabularCooling`** (`src/hydro/srcterms/tabular_cooling.cpp`, `<cooling>` block):
   upstream AthenaPK optically-thin tabulated Λ(T) cooling with adaptive RK12/RK45
   subcycling — built for the cluster/ICM problems, **compiled in but not used in any
   FHC run** (`cooling/enable_cooling = none` default; the collapse decks never set it).

The chemistry package deliberately does **no thermal chemistry** (no line cooling /
photoelectric heating): `network_gow17_reduced.hpp` states its purpose is only the
electron abundance, not the gow17 thermal network.

---

## 5. Magnetohydrodynamics core

`src/hydro/hydro.cpp` assembles the solver from input keys:

- `hydro/fluid = euler | glmmhd`. GLM-MHD is cell-centered with Dedner hyperbolic
  divergence cleaning: ψ advects ∇·B error at speed c_h (recomputed each step from the
  global min dx) and damps it with factor exp(−α c_h Δt/dx_min)
  (Mignone & Tzeferacos 2010 eq. 27; `glmmhd/dedner_source.cpp`).
  `glmmhd_source = dedner_plain | dedner_extended` — extended adds the non-conservative
  Powell-like −(∇·B)B momentum term and −B·∇ψ energy term for robustness (production
  uses `dedner_extended`, `glmmhd_alpha = 0.1`).
- `hydro/riemann = llf | hlle | hllc | hlld | none` (production: **HLLD**).
- `hydro/reconstruction = dc | plm | ppm | limo3 | weno3 | wenoz` (production: **PLM**;
  needs `nghost≥2`, wenoz/ppm need 3).
- `parthenon/time/integrator = rk1 | rk2 | rk3 | vl2` (production: **VL2**, predictor
  stage automatically dropped to donor-cell reconstruction).
- `hydro/first_order_flux_correct = true|false` — retry troubled cells with first-order
  fluxes.
- Timestep: CFL (`parthenon/time/cfl = 0.3`) on the fast magnetosonic speed, minimum over
  cells; separate diffusive and radiation constraints below.

Conserved variable layout (glmmhd): [ρ, ρv1..3, E, B1..3, ψ, scalars...]; `prim` mirrors
with [ρ, v, p, B, ψ, x_i].

---

## 6. Self-gravity (`src/self_gravity/`, ported from Artemis/Parthenon poisson_gmg)

**Equation.** ∇²Φ = 4πG(ρ − ρ̄) with `four_pi_G` from input (=1 in these units); ρ̄ is
subtracted only under the Jeans swindle (`use_swindle`, required for fully periodic
meshes, computed by a global mass/volume reduction each step).

**Discretization** (`poisson_equation.hpp`): conservative 7-point flux form. Face flux
F = −∇Φ by 2-point difference; A·Φ = div(F) with AMR flux correction at coarse–fine
boundaries on the leaf grid (disabled inside MG two-level composites, which have their
own coarse–fine handling). RHS assembled from `prim` density every step over the entire
domain including ghosts (`FillPoissonRHS`).

**Solver.** Runtime-selectable (`<self_gravity/solver_params> solver`):
- `BiCGSTAB` (default, production): Krylov outer loop preconditioned by geometric
  multigrid; tolerates a weak SRJ1 (one weighted-Jacobi sweep) smoother.
- `MG`: pure GMG V-cycles, no global inner products (latency-optimized), but needs
  SRJ2/SRJ3 smoothing on AMR grids (SRJ1 V-cycles diverge across fine–coarse boundaries).

Production settings: `max_iterations=200`, `residual_tolerance=1e-6`, `do_FAS=true`,
`max_coarsenings=20`, smoothers SRJ1, `preconditioner=Multigrid`.

**Boundary conditions** (`ix1_bc..ox3_bc`, per face): `zero` (Φ=0 Dirichlet via
FixedFace odd reflection), `neumann` (∂Φ/∂n=0, outflow copy), or `default`. **Multipole
BCs do not exist in this port and Neumann is unreliable — production uses `zero`
Dirichlet on all six faces**, which is why the box is padded to L=52 (≈4× sphere
diameter): the Φ=0 wall clips the potential in tight boxes (the L=16 failure was config,
not code). A packed whole-MeshData BC path (`packed_bc=true`, bit-identical to the
per-block path) removes the per-block kernel-launch bottleneck on GPU.

**Coupling to hydro** (`ApplyGravitySource`, called per stage with β·dt): momentum gets
ρ·g with g = −∇Φ by centered differences; energy gets the *flux-weighted work*
Σ_faces (mass flux)·ΔΦ·(½dt/dx) (Artemis style — uses the hydro mass fluxes still in
memory, which conserves total energy to the truncation of the flux divergence rather
than dotting g into cell-centered momentum).

**Inputs** (`<self_gravity>`): `units_override` (must be true), `four_pi_G` (1.0),
`use_swindle` (auto-true iff fully periodic), six `*_bc` keys, `packed_bc`, plus the
`solver_params` sub-block above. Requires `parthenon/mesh/multigrid = true`.

---

## 7. Non-ideal MHD (`src/hydro/diffusion/`)

All three terms are cell-centered flux/EMF additions computed in `CalcDiffFluxes`
(`diffusion.cpp`) after the Riemann fluxes; face currents J = ∇×B use 2-point normal
derivatives and 4-point averaged transverse derivatives; face B and η are arithmetic
face averages. Each EMF adds antisymmetrically to the two transverse B fluxes plus the
corresponding Poynting flux (E×B) to the energy flux — exactly the Athena++
field-diffusion structure, adapted to the GLM cell-centered grid.

### 7.1 The three terms

| term | EMF | diffusivity (fixed mode) | dt limit (3D) |
|---|---|---|---|
| Ohmic (`resistivity.cpp`) | η_O J | η_O = `ohm_diff_coeff_code` (const) | (1/6)·dx²/η_O |
| Ambipolar (`ambipolar.cpp`) | η_A (J − (J·b̂)b̂) | η_A = Q_A·B², Q_A = `ambipolar_coeff_code` | (1/6)·dx²/η_A |
| Hall (`hall.cpp`) | η_H (J×B)/\|B\| + η_floor J | η_H = Q_H·B/ρ, Q_H = `hall_coeff_code` | 1.0·dx²/\|η_H\| (whistler) + parabolic for η_floor |

All limits get the `diffusion/cfl` safety factor (default = hyperbolic CFL; t4: 0.3).
The fixed AD coefficient maps **1:1** to Athena++'s `eta_ad` (both η_A = coeff·B²).

**Integrators** (`diffusion/integrator`):
- `unsplit`: diffusive fluxes added to every stage's flux divergence; global dt takes
  min(hyperbolic, diffusive).
- `rkl2`: RKL2 super-time-stepping, Strang-split as two half-dt STS blocks around the
  hydro step (`hydro_driver.cpp:460,623`), with `rkl2_max_dt_ratio` capping the
  hyperbolic:parabolic ratio.
- **Hall is dispersive (whistler waves) ⇒ enforced `unsplit` only** (parse-time error
  otherwise), and warns unless Ohmic resistivity or `hall_ohmic_floor_code` > 0 is
  present (production floor: 0.05). Validation (recalled from DEV_LOG/memory): AD
  damped-Alfvén eigenmode to 0.04%, Hall whistler dispersion to 0.4%.

### 7.2 Coefficient modes

Each term independently selects `<term>_coeff = fixed | ionization | ionization_chem`:

- **fixed**: constants above (this is the Athena++-matched comparison mode; production
  AD-only tiers used Q_A = 0.1).
- **ionization**: self-consistent η(ρ,T,B) from the reduced ionization model (§7.3),
  evaluated per face (or once per cell with the eta cache).
- **ionization_chem**: same tensor machinery but the electron abundance comes from the
  time-dependent chemistry scalar (prim component `nhydro + xe_scalar_index`, default
  index 4 = the gow17_reduced electron). For AD this uses the single-fluid
  η_A = B²/(4π γ_AD ρ_i ρ) with ρ_i = x_e n_n m_ion, **capped at the equilibrium-model
  η_A** — the chemistry x_e collapses at the C→CO ionization minimum (~1e-13) and the
  uncapped single-fluid η_A blow-up was the dt-collapse that froze the RT+chem run.

Production t4 mixes modes: AD `ionization_chem`, Ohmic `ionization`, Hall `ionization`
(signed tensor Hall; no Athena++ counterpart — Athena++'s HallEMF is a commented-out stub).

**Performance plumbing:** `PrecomputeNonidealEta` fills a 3-component cell-centered cache
(η_O, η_H, η_A) once per cell per step (flag `diffusion/eta_cache`, ~1.7× speedup
recalled), and `EstimateNonidealTimestepIonizationFused` computes the Wardle tensor once
per cell for the combined dt reduction instead of three times (auto-enabled only when
every active term uses plain `ionization`).

### 7.3 Reduced ionization microphysics (`ionization.hpp`, NICIL-class)

Computes the charge state of dense molecular gas and from it the conductivity tensor:

**Charge carriers**: electrons, one representative ion (m_ion = 24.3 m_H, ~Mg+/HCO+),
and 5 MRN grain size bins (dn/da ∝ a^−3.5 over [a_min, a_max] = [5 nm, 250 nm],
f_dg = 0.01, ρ_grain = 3 g/cm³).

**Ionization balance** (`SolveCharges`): cosmic-ray pair production P = ζ n_n
(ζ = 1e-17 s⁻¹) balanced against dissociative recombination
α(T) = 2.4e-7 (T/300)^−0.69 cm³/s and capture onto grains; each grain bin's mean charge
⟨Z_k⟩ from OML electron/ion capture balance (Newton on the reduced potential ψ);
charge neutrality n_e = n_i − Σ(−Z_k)n_g,k; fixed-point iteration with 0.5 relaxation.
On top: thermal Saha ionization of potassium (x_K = 1e-7, χ = 4.34 eV — the ~10³ K
electrons that re-couple the field after grains sublimate at T_subl = 1500 K) **and** of
hydrogen (13.6 eV, dominant >5000 K into the second core), sharing n_e via bisection.

**Diffusivities** (Pandey & Wardle 2008 conductivity tensor):

```
σ_{O,H,P} = (ec/B) Σ_j n_j Z_j · {β_j, 1/(1+β_j²), β_j/(1+β_j²)},  β_j = Z_j e B/(m_j c ν_jn)
η_O = c²/(4π σ_O);   η_H = c²/(4π) σ_H/σ_⊥²  (SIGNED);   σ_⊥² = σ_H² + σ_P²
```

Ambipolar closure selectable (`ion_ad_closure`):
- `single_fluid` (default, Athena++-matched): η_A = B²/(4π γ_AD ρ_i ρ),
  γ_AD = 3.5e13 cm³ g⁻¹ s⁻¹;
- `tensor`: η_A = c²/(4π)·σ_P/σ_⊥² − η_O (grain-modified, ~3–5× lower, floored at 0).

Ohmic and Hall are *always* tensor-derived — the electron-only single-fluid Hall
over-predicts catastrophically at the ionization minimum (the bug the tensor fixed).

**The ionization story across the collapse.** (i) Diffuse envelope: CR ionization vs
gas-phase recombination gives x_e ≈ √(ζ/(α n_n))/n_n^0 ∼ n^(−1/2) — well coupled.
(ii) Dense core (n ≳ 10⁸ cm⁻³): free electrons freeze onto grains, the negatively
charged MRN grains (⟨Z_k⟩ < 0 from OML capture balance — electrons arrive faster than
ions) become significant charge carriers, x_e passes through the **ionization minimum**
where all three η peak and the field decouples (this is where the AD/Hall/Ohmic terms
matter most for the B-field structure of the FHC). (iii) T ≈ 1500 K: grains sublimate
(`T_subl`) — they leave both the charge balance *and* the opacity — briefly raising the
resistivity further. (iv) T ≳ 10³ K: Saha ionization of potassium (χ = 4.34 eV,
x_K = 1e-7) supplies electrons and starts re-coupling; K saturates at x_e ~ 1e-7.
(v) T ≳ 5000 K into the second core: hydrogen Saha (13.598 eV, same expression as the
EOS table so conductive and thermodynamic ionization agree) takes over and restores
near-ideal MHD. The two Saha donors share n_e via a bisection on charge neutrality
(`SahaThermal`).

**Chemistry-coupled variant.** With `ionization_chem`, the *free-electron density* is
taken from the time-dependent gow17_reduced x_e scalar instead of the CR-equilibrium
solve; grains are still charged, but re-equilibrated at that fixed n_e
(`SolveGrainsFixedNe`: n_i = n_e + Σ(−Z_k)n_g,k closes neutrality). This captures
recombination lag in fast-collapsing gas and removes the equilibrium solver's
non-convergence noise at the minimum; the AD cap against the equilibrium η_A (§7.2)
bounds the runaway where the chemistry x_e collapses to ~1e-13.

Note there are **three separate "ionization" computations** in the code, kept mutually
consistent by construction: this conductivity model (η's), the chemistry network's x_e
(non-equilibrium, feeds this model), and the EOS table's Saha ionization
(thermodynamics: p, c_s², γ_eff). The H-Saha expressions in (v) and in the EOS table are
deliberately identical.
Internally everything is evaluated in cgs (unit scales rho_unit = 5.467e-19,
T_unit = 10.015 K, B_unit = 4.98e-5 G [Gaussian G for the microphysics],
η_unit = 1.874e-21) and returned in code units. Under `eos=hydrogen` the input
temperature comes from the EOS table (T ≠ p/ρ when μ varies).

**Inputs** (`<diffusion>` `ion_*` keys, all defaulted to the FHC calibration):
`ion_zeta`, `ion_mu_n`, `ion_m_ion`, `ion_alpha0`, `ion_sigma_en` (1e-15 cm²),
`ion_sigv_in` (1.9e-9 cm³/s Langevin), `ion_gamma_AD`, `ion_f_dg`, `ion_rho_grain`,
`ion_T_subl`, `ion_a_min_cm`, `ion_a_max_cm`, `ion_mrn_p`, `ion_x_K`, `ion_chi_K`,
`ion_ad_closure`, plus the four unit conversions `ion_rho_unit_cgs / ion_T_unit_K /
ion_B_unit_G / ion_eta_unit`.

(Also present but unused in the collapse problem: isotropic/anisotropic **thermal
conduction** with fixed or Spitzer coefficient and Cowie–McKee saturation, and isotropic
**viscosity** with fixed coefficient — same `<diffusion>` block.)

---

## 8. Radiation transport: gray M1 (`src/radiation/`, ported from Artemis)

**Physics.** Two-moment gray radiation hydrodynamics: evolve radiation energy density
E_r and flux F_r with the M1 closure — the pressure tensor P_ij = D_ij(f)E_r
interpolates between Eddington (P = E/3, optically thick) and free-streaming
(P = n̂n̂E, |F| → ĉE) using the flux factor f = |F|/(ĉE). This captures both the
optically-thin envelope cooling and the radiative trapping that sets the first-core
entropy, replacing the barotropic law.

**Transport** (`radiation_moments.cpp`): donor-cell HLL Riemann fluxes for the 4-vector
(E_r, F_r) with M1 characteristic speeds; explicit FV update sub-cycled at the radiation
CFL (`cfl=0.4` · dx_min/ĉ) within each hydro dt (flux → update → ghost exchange per
sub-step); causal clamp |F_r| ≤ ĉE_r and floor E_r ≥ `efloor` after each update.
Reduced speed of light (RSLA) is mandatory: c/v0 ≈ 1.58e6; production `creduc = 1000`
(ĉ ≈ 1578 code).

**Matter coupling** (once per hydro dt, per cell, implicit): Newton iteration on
(E, B≡a_rad T⁴) for the linearized absorption–emission system with Fleck factor,

```
S_E = ĉ ρκ_a (E − B),   S_e = −(c/ĉ)·S_E,   convergence to inner_iteration_tol=1e-8
```

then explicit flux attenuation ΔF = −a/(1+a)·F with a = ĉ dt ρ(κ_a+κ_s), and the
**radiation force**: gas momentum gains −ΔF/(c·ĉ) (momentum-conserving sign of Artemis'
"Full" coupling; the original "Simple" port had the sign backwards — fixed in the 2026-07
audit). Radiation loses (ĉ/c)(ΔE_thermal + ΔE_kinetic), floored at `efloor`.
**When RT is active it owns the gas thermal energy** — the barotropic overwrite is
skipped and T comes either from the ideal law e = ρT/(γ−1) or from the EOS table
(with numerical C_v = de/dT).

**Opacity** (`radiation_opacity.hpp`, `opacity_model =`):
- `constant`: `kappa_a_code`, `kappa_s_code` (code units, test problems).
- `dust`: Bell & Lin (1994) ice-grain law only, κ = `dust_kappa0_cgs`·T² cm²/g
  (default k0 = 2e-4), **frozen** at its value at `dust_Tmax_K = 1500` K above that —
  adequate through first-core temperatures, wrong beyond.
- `belllin` (production): the full 8-regime Bell & Lin Rosseland-mean fit
  κ = k0 ρ^a T^b, regime chosen by walking up the density-dependent crossing
  temperatures of adjacent regimes:

  | regime | k0 | ρ-exp a | T-exp b |
  |---|---|---|---|
  | ice grains | 2e-4 | 0 | +2 |
  | ice evaporation | 2e16 | 0 | −7 |
  | metal/dust grains | 0.1 | 0 | +0.5 |
  | **dust sublimation gap** | 2e81 | 1 | **−24** |
  | molecular | 1e-8 | 2/3 | +3 |
  | H⁻ | 1e-36 | 1/3 | +10 |
  | Kramers (bf+ff) | 1.5e20 | 1 | −2.5 |
  | electron scattering | 0.348 | 0 | 0 |

**Why the dust law matters physically.** Below ~150 K the opacity is entirely
ice-mantled grains with κ ∝ T² (Rayleigh-limit grain emissivity). This is what makes
the collapsing core optically thick at ρ ≈ 1e-13 g/cm³ — radiative trapping ends the
isothermal phase and forms the first core (the barotropic ρ_crit is precisely a proxy
for this opacity physics; in RT runs it emerges instead of being imposed). The
sublimation gap at 1500–2000 K (κ plunging as T⁻²⁴) transiently *releases* radiation
just before H2 dissociation triggers second collapse, and dust destruction here is the
same `T_subl ≈ 1.5e3` K event that removes grains from the ionization charge balance
(§7.3) — the code applies it consistently in both places.

Gray absorption coefficient used by the coupling: α = ρ_code·κ_code with
κ_code = κ_cgs·(rho_unit·l_unit) — for the FHC units the conversion is
5.467e-19 × 2.806e16 = 1.534e-2 per (cm²/g). Scattering κ_s is a constant
(`kappa_s_code`, 0 in production); dust scattering is not modeled.

**Code-unit calibration** (no free knobs): a_rad_code = a_R T_unit⁴/e_unit,
c_code = c/v0, κ_code = κ_phys·(rho_unit·l_unit); overridable for tests.

**Inputs** (`<radiation>`): `closure` (M1|P1), `creduc`, `cfl`, `efloor`,
`tfloor_code`, `matter_coupling`, `mu`, `rho_unit_cgs`, `v_unit_cgs`,
`length_unit_cgs`, opacity keys above, `inner_iteration_max/tol`. Enabled by
`<physics> radiation = true`. The M1 fields carry an `OperatorSplit` metadata flag so the
hydro integrator's flux-divergence sweep leaves them alone.

---

## 9. Chemistry (`src/chemistry/`, greenfield GPU port)

Species ride on AthenaPK passive scalars (`hydro/nscalars`): conserved
scalar_density_i = ρ·x_i, so prim scalar = abundance x_i = n_i/n_H. Reactions are
operator-split at fixed ρ (and T), once per hydro dt. Two networks
(`chemistry/network =`):

**"H2"** (NSPEC=2: H, H2) — minimal formation/destruction:

```
dx_H2/dt = k_gr n_H x_H − k_cr x_H2       (k_gr = 3e-17 cm³/s grain formation,
dx_H/dt  = −2(...)+2(...)                  k_cr = 3·ζ, ζ = 2e-16 s⁻¹)
```

explicit sub-cycled Euler with per-species CFL (`cfl_cool=0.1`, `nsub_max=200`).

**"gow17_reduced"** (NSPEC=5: H2, H+, C+, CO, e⁻) — compact reduction of Gong, Ostriker
& Wolfire (2017) whose sole purpose is a time-dependent **electron abundance** for the
non-ideal diffusivities. Eight reactions (H2 grain formation / CR dissociation, H CR
ionization / recombination gas+grain, C ionization ζ_C = 3ζ / C+ radiative recombination
a_C = 4.7e-12(T/300)^−0.6, CO formation k_CO = 5e-16 n_H / CO CR destruction 10ζ), with
conservation constraints x_H + 2x_H2 + x_H+ = 1 and x_C+ + x_CO + x_C0 = x_Ctot
(=1.6e-4) enforced by joint renormalization, and x_e = x_H+ + x_C+ slaved to charge
neutrality (floor `xe_floor = 1e-15`). Integrated with the **semi-implicit
production/loss split** y_new = (y + dt·P)/(1 + dt·L) — unconditionally stable and
positivity-preserving, so stiffness no longer caps the sub-step (accuracy-limited,
`nsub_max=400` guaranteed to cover dt). T for the rates comes from the conserved energy
(ideal p = ρT). Cells that truncate are counted and warned once per run.

**Inputs** (`<chemistry>`): `network`, `zeta_cr_cgs`, `kgr_cgs`, `x_Ctot`, `xe_floor`,
`x_floor`, `nsub_max`, `cfl_cool`, unit overrides. Enabled by
`<physics> chemistry = true`; requires `nscalars ≥ NSPEC`. Coupling to the MHD is via
`<diffusion> *_coeff = ionization_chem` + `xe_scalar_index = 4`.

---

## 10. Adaptive mesh refinement

`<refinement> type = jeans` (`src/refinement/jeans.cpp`) — Truelove (1997) criterion.
With 4πG = 1 the Jeans length is

```
λ_J = 2π (c_s + v_A)/√ρ,    c_s = √(γp/ρ),  v_A = √(B²/ρ)  (HL; MHD-augmented, Athena++ convention)
```

Refine when λ_J/dx < `njeans` anywhere in the block; derefine when > 2.5·njeans
everywhere. Production: `njeans = 8`, `numlevel = 20` (finest dx = 3.9e-7 code
= 0.16 R☉, holds Truelove through second-core densities; comment in fhc.in),
`derefine_count = 50` to damp AMR flapping (the t1 OOM was churn-buffer growth).
The gravity solver participates in AMR via GMG restriction/prolongation ops registered
on Φ and rhs.

---

## 11. Orchestration: one timestep (t4 configuration)

From `hydro_driver.cpp` task graph:

1. (If `rkl2`) first half-dt STS block for the parabolic terms.
2. For each VL2 stage: reconstruct → Riemann fluxes (+ first-order flux correction) →
   `CalcDiffFluxes` adds unsplit non-ideal EMF/Poynting fluxes → flux divergence →
   unsplit sources (Dedner GLM, `ApplyBarotropicCooling` [skipped-energy variant when RT
   active], per-stage `SelfGravity::ApplyGravitySource` with the fresh Φ) → c2p.
3. `SelfGravity::SolvePoisson` (BiCGSTAB/MG) once per step on the updated density.
4. `Radiation::AddRadiationTasks`: sub-cycled M1 transport + one implicit matter
   coupling over the full dt.
5. `Chemistry::AddChemistryTasks`: per-cell network integration over dt.
6. (If `rkl2`) second half-dt STS block.
7. New dt = min(hyperbolic CFL, diffusive limits §7, `hydro/max_dt` if set); radiation
   sub-cycles internally so ĉ does not set the global dt.

---

## 12. Known caveats (as encoded in the source/comments)

- Gravity Neumann/multipole BCs: unsupported/broken — use `zero` Dirichlet with a padded
  box (§6).
- Hall: unsplit only; needs the Ohmic floor (~0.05) on the cell-centered grid.
- `ionization_chem` AD is capped at the equilibrium η_A to prevent the
  ionization-minimum dt collapse.
- The barotropic source *overwrites* e_th — any physical heating (shocks) below ρ_crit
  is discarded by construction in non-RT runs.
- Global dt still collapses near first-core formation (rhocrit=1e-13); sink particles
  are not implemented (dropped from the port); no flux-limited-diffusion mode — M1 only.
- Athena++ comparison: identical code-unit B numbers, but Athena++ is Gaussian
  (v_A = B/√(4πρ)) — the √4π is the standing analysis landmine.

## Source-file map

| physics | files |
|---|---|
| BE-sphere IC + barotropic source | `src/pgen/collapse_be.cpp` |
| MHD core, input parsing | `src/hydro/hydro.cpp`, `src/hydro/hydro_driver.cpp`, `src/hydro/glmmhd/`, `src/hydro/rsolvers/`, `src/recon/` |
| Self-gravity | `src/self_gravity/{self_gravity.cpp, poisson_equation.hpp, self_gravity_driver.cpp}` |
| Non-ideal MHD | `src/hydro/diffusion/{diffusion.cpp/.hpp, resistivity.cpp, ambipolar.cpp, hall.cpp}` |
| Ionization microphysics | `src/hydro/diffusion/ionization.hpp` |
| M1 radiation | `src/radiation/{radiation.cpp, radiation_moments.cpp, radiation_closure.hpp, radiation_opacity.hpp}` |
| Chemistry | `src/chemistry/{chemistry.cpp, network_h2.hpp, network_gow17_reduced.hpp}` |
| Tabulated EOS | `src/eos/{gen_eos_table.py, eos_table.hpp}`, `src/eos/eos_table.bin` |
| AMR criterion | `src/refinement/jeans.cpp` |
| Production deck | `runs/prod_t4_full/fhc.in`, `runs/prod_t4_full/units.json` |
