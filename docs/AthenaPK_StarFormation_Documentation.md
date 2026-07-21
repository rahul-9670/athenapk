# AthenaPK for Star Formation — A Complete Technical Documentation

*A self-contained, thesis-style account of the AthenaPK code: its architecture, its
numerical engine, and the implementation of every physical process relevant to the
collapse of molecular cloud cores into stars — hydrodynamics, magnetohydrodynamics,
self-gravity, driven turbulence, ambipolar diffusion, the Hall effect, Ohmic
resistivity, cooling, and adaptive refinement. Every algorithm is traced down to the
formula and the lines of code that implement it.*

---

## Table of Contents

**Part I — Foundations**
1. Prologue: the physics of star formation and why we need a code
2. What AthenaPK is: lineage, philosophy, and the Parthenon/Kokkos stack
3. Code architecture: packages, tasks, meshes, and data layout
4. The governing equations: ideal MHD in conservative form

**Part II — The numerical engine**
5. The finite-volume Godunov method
6. Reconstruction
7. Riemann solvers
8. The equation of state and wave speeds
9. Time integration: the van-Leer predictor–corrector
10. Divergence control: GLM / Dedner cleaning
11. Adaptive mesh refinement and the Jeans criterion

**Part III — The physics of star formation, as implemented**
12. Self-gravity: the Poisson problem
13. Driven turbulence
14. Thermodynamics: the barotropic EOS and the first hydrostatic core
15. The non-ideal MHD framework
16. Ohmic resistivity
17. Ambipolar diffusion
18. The Hall effect
19. Optically-thin cooling
20. Super-time-stepping (RKL2)

**Part IV — Synthesis**
21. Anatomy of a single timestep
22. A collapse simulation end-to-end
23. Units
24. Performance portability: the GPU story

---

# Part I — Foundations

## 1. Prologue: the physics of star formation and why we need a code

Stars are born inside **giant molecular clouds** — cold (≈10 K), dense (n ≈ 10²–10⁴
cm⁻³), turbulent, magnetized reservoirs of mostly molecular hydrogen. Left to itself, a
parcel of such gas is held up against its own gravity by a combination of thermal
pressure, supersonic turbulence, and magnetic fields. When a region becomes massive
enough that gravity overwhelms these supports — quantified by the **Jeans criterion** —
it begins to collapse.

The collapse is not a simple free-fall. It is a multi-physics drama unfolding over many
orders of magnitude in density and length scale:

- **Isothermal phase.** While the gas is optically thin, the compressional heating
  radiates away almost instantly. The gas stays at ≈10 K and collapses nearly
  isothermally, the density running away as ρ → ∞ in finite time at the centre.
- **First hydrostatic core (FHSC).** Once the central density reaches ρ ≈ 10⁻¹³ g cm⁻³,
  the core becomes optically thick to its own dust emission. Heat can no longer escape;
  the gas heats up, pressure rises faster than gravity, and a tiny (~AU-scale),
  pressure-supported **first core** forms and bounces. This is the birth of the
  protostar's hydrostatic seed.
- **Magnetic mediation.** The cloud's magnetic field is dragged inward with the
  collapsing, partially-ionized gas. In ideal MHD, magnetic flux is "frozen" to the
  fluid, and the contracting field would brake rotation so efficiently that no
  rotationally-supported disc could ever form — the **magnetic braking catastrophe**. The
  resolution lies in **non-ideal MHD**: the gas is only ~1 part in 10⁷–10⁸ ionized, so
  the neutrals and the (field-tied) ions are imperfectly coupled. **Ambipolar diffusion**
  lets neutrals slip across field lines; the **Hall effect** introduces a handedness that
  can spin gas up or down depending on field/rotation alignment; **Ohmic resistivity**
  dissipates field in the densest, best-shielded regions. Together they regulate how much
  flux is lost, whether discs form, and how outflows are launched.

To follow this story quantitatively one must solve, simultaneously and over a dynamic
range of >10⁵ in spatial scale, the equations of **compressible magnetohydrodynamics**
coupled to **self-gravity**, with a realistic **thermodynamic** treatment and the three
**non-ideal MHD** effects, all while resolving the ever-shrinking **Jeans length**. That
is precisely what AthenaPK is built to do.

This document explains how.

## 2. What AthenaPK is

**AthenaPK** is a performance-portable, block-structured adaptive-mesh-refinement (AMR)
astrophysical MHD code. Its physics and numerical methods are inherited from the widely
used **Athena++** code (Stone et al. 2020), but its execution model is entirely
rebuilt on two modern libraries:

- **Kokkos** — a C++ performance-portability layer. The same kernel source compiles and
  runs efficiently on multicore CPUs, NVIDIA GPUs, AMD GPUs, etc. AthenaPK kernels are
  written once as `par_for` / `par_reduce` loops with `KOKKOS_LAMBDA` bodies; Kokkos maps
  them onto whatever backend the binary was built for.
- **Parthenon** — a block-AMR framework (a descendant of Athena++'s mesh, K-Athena, and
  Parthenon-proper) that provides the mesh, the AMR machinery, MPI communication, I/O,
  and a **task-based execution model**.

The practical consequence, and the reason AthenaPK exists, is the GPU: a single
simulation can run on thousands of GPUs with the identical source that runs on a laptop
CPU. The trade-off — important for this project — is that GPUs reward large,
flop-dense, regular work and punish latency-bound, communication-heavy work such as the
iterative Poisson solve of self-gravity (see §24).

The directory layout of the source (`src/`) mirrors the physics:

```
src/
  main.cpp, main.hpp        global enums and variable indices
  hydro/                    the core MHD package
    hydro.cpp/.hpp          package setup, flux driver, timestep
    hydro_driver.cpp        the per-cycle task collection (the "conductor")
    rsolvers/               Riemann solvers (hlle, hllc, hlld, llf, ...)
    glmmhd/                 GLM divergence-cleaning source
    diffusion/              non-ideal MHD: resistivity, ambipolar, hall, viscosity, conduction
    srcterms/               gravitational field, tabular cooling
    prolongation/           custom AMR prolongation ops for MHD
  recon/                    reconstruction (dc, plm, ppm, weno3, wenoz, limo3)
  eos/                      adiabatic hydro & GLM-MHD equations of state
  self_gravity/             Poisson solver (multigrid-preconditioned BiCGSTAB)
  refinement/               AMR criteria (gradient, jeans, ...)
  utils/few_modes_ft.*      few-mode Fourier turbulence driver
  pgen/                     problem generators (collapse_be, turbulence, sod, ...)
```

## 3. Code architecture: packages, tasks, meshes, and data layout

### 3.1 The mesh hierarchy

Parthenon decomposes the domain into a forest of **MeshBlocks** — small logically-
Cartesian sub-grids (e.g. 32³ cells), each padded by **ghost zones** (`Globals::nghost`,
typically 2) that hold copies of neighbouring data for stencils. AMR proceeds by
**refining** a block into 2ⁿᵈⁱᵐ children at half the cell size, organised as an octree.
A "level" is one factor-of-two in resolution; a star-formation run uses `numlevel` ≈ 16
levels, a dynamic range of 2¹⁵ ≈ 32000 between coarsest and finest cells.

For computation, blocks are gathered into **packs**: a `MeshData<Real>` is a pack of many
blocks presented to a kernel as a single multi-dimensional array indexed
`(b, var, k, j, i)` where `b` is the block within the pack. This amortizes kernel-launch
overhead — essential on GPUs.

### 3.2 Data layout: primitive and conserved variables

Two state vectors are carried for every cell, defined in `src/main.hpp`:

```cpp
enum ConsIndex { IDN = 0, IM1 = 1, IM2 = 2, IM3 = 3, IEN = 4, NHYDRO = 5,
                 IB1 = 5, IB2 = 6, IB3 = 7, IPS = 8 };   // conserved
enum         { IV1 = 1, IV2 = 2, IV3 = 3, IPR = 4 };     // primitive velocity/pressure
```

- **Conserved** `cons` = (ρ, ρv₁, ρv₂, ρv₃, E, B₁, B₂, B₃, ψ). These are what the
  finite-volume scheme actually evolves, because their volume integrals obey exact
  conservation laws. `E` is the total energy density; `ψ` is the GLM scalar (§10).
- **Primitive** `prim` = (ρ, v₁, v₂, v₃, P, B₁, B₂, B₃, ψ). These are what the physics
  (reconstruction, Riemann fluxes, wave speeds) is naturally expressed in.

For pure hydrodynamics `NHYDRO = 5`; GLM-MHD adds the three field components plus ψ for
nine variables (`NGLMMHD = 9`). The `Fluid` enum selects the model:

```cpp
enum class Fluid { undefined, euler, glmmhd };
```

The conversion cons → prim ("FillDerived") is done by the EOS package every time the
conserved variables change.

### 3.3 The task-based execution model

A timestep is not a monolithic function. Instead, `HydroDriver::MakeTaskCollection`
(`src/hydro/hydro_driver.cpp`) builds a **directed acyclic graph of tasks** —
`reconstruct → Riemann flux → flux-correct at coarse/fine boundaries → update → source
terms → boundary exchange → ...` — and Parthenon executes it, overlapping computation
with MPI communication wherever the dependency graph allows. Tasks are grouped into
**regions**; within a region, one task list runs per block-pack, and independent lists
may run concurrently. We dissect this graph in §21.

## 4. The governing equations: ideal MHD in conservative form

AthenaPK evolves the equations of compressible, inviscid, ideal MHD in
conservation form. In Heaviside–Lorentz units (the magnetic-pressure factor is B²/2, with
**no** 4π — see §23), with gravitational potential Φ:

**Mass:**
$$\partial_t \rho + \nabla\cdot(\rho\mathbf{v}) = 0$$

**Momentum:**
$$\partial_t (\rho\mathbf{v}) + \nabla\cdot\!\left[\rho\mathbf{v}\mathbf{v} - \mathbf{B}\mathbf{B} + \left(P + \tfrac{1}{2}B^2\right)\mathsf{I}\right] = -\rho\nabla\Phi$$

**Total energy** (E = P/(γ−1) + ½ρv² + ½B²):
$$\partial_t E + \nabla\cdot\!\left[\left(E + P + \tfrac12 B^2\right)\mathbf{v} - (\mathbf{v}\cdot\mathbf{B})\mathbf{B}\right] = -\rho\mathbf{v}\cdot\nabla\Phi + (\text{non-ideal, cooling})$$

**Induction:**
$$\partial_t \mathbf{B} - \nabla\times(\mathbf{v}\times\mathbf{B}) = -\nabla\times\mathbf{\mathcal E}_{\rm non\text{-}ideal}$$

with the divergence constraint ∇·**B** = 0. The non-ideal EMF collects the three effects
documented in §16–18:
$$\mathbf{\mathcal E}_{\rm non\text{-}ideal} = \underbrace{\eta_O\mathbf{J}}_{\text{Ohmic}} + \underbrace{\eta_A\,\mathbf{J}_\perp}_{\text{ambipolar}} + \underbrace{\eta_H\,\frac{\mathbf{J}\times\mathbf{B}}{|\mathbf B|}}_{\text{Hall}},\qquad \mathbf{J}=\nabla\times\mathbf{B}.$$

The Poisson equation closes the gravity coupling:
$$\nabla^2\Phi = 4\pi G\rho.$$

The remainder of this document is the story of how each term on these right-hand sides is
discretized and solved.

---

# Part II — The numerical engine

## 5. The finite-volume Godunov method

AthenaPK is a **Godunov finite-volume** code. The domain is tiled by cells; the code
stores **volume averages** of the conserved variables. Integrating the conservation law
∂ₜU + ∇·F(U) = S over a cell and a timestep gives the exact update

$$\bar U_i^{n+1} = \bar U_i^{n} - \frac{\Delta t}{\Delta x}\!\left(F_{i+1/2} - F_{i-1/2}\right) + \Delta t\,\bar S_i,$$

so the entire problem reduces to estimating the **interface fluxes** F at cell faces. The
Godunov recipe for each face is three steps:

1. **Reconstruct** left/right states at the face from cell averages (§6).
2. Solve a **Riemann problem** between them to get a single upwind flux (§7).
3. **Update** the conserved variables by the flux divergence, then add **source terms**
   (gravity, GLM, non-ideal, cooling) (§9–§20).

In code, the per-face flux for one direction is computed in `Hydro::CalculateFluxes`
(`src/hydro/hydro.cpp`), which loops over a block-pack, calls the templated
`Reconstruct<recon, DIR>` into scratch arrays `wl`, `wr`, and then the templated
`Riemann<fluid, solver>::Solve` to write `cons.flux(dir, var, k, j, i)`. The flux
divergence and update are Parthenon library calls
(`Update::UpdateWithFluxDivergence`), wired into the task graph (§21).

## 6. Reconstruction

The accuracy and robustness of a Godunov scheme is set largely by **reconstruction** —
how the smooth sub-cell profile is rebuilt from discrete averages. AthenaPK offers a
menu (`enum class Reconstruction { dc, plm, ppm, wenoz, weno3, limo3 }`):

- **dc** — donor cell (1st order, diffusive, robust).
- **plm** — piecewise-linear with a slope limiter (2nd order); the workhorse for collapse.
- **ppm** — piecewise-parabolic (3rd order).
- **weno3 / wenoz / limo3** — high-order (essentially) non-oscillatory schemes.

The default for star-formation runs is **PLM**, which balances accuracy and cost and is
total-variation-diminishing (no spurious oscillations at shocks). Its core is a
**van-Leer harmonic-mean limiter** (`src/recon/plm_simple.hpp`):

```cpp
KOKKOS_INLINE_FUNCTION
void PLM(const Real &q_im1, const Real &q_i, const Real &q_ip1, Real &ql_ip1, Real &qr_i) {
  Real dql = (q_i - q_im1);          // left slope
  Real dqr = (q_ip1 - q_i);          // right slope
  Real dq2 = dql * dqr;
  Real dqm = 0.0;
  if (dq2 > 0.0) {                   // only if slopes agree in sign (monotone)
    dqm = dq2 / (dql + dqr);         // harmonic mean -> van Leer limiter
  }
  ql_ip1 = q_i + dqm;                // right-face value of cell i
  qr_i   = q_i - dqm;                // left-face value of cell i
}
```

The logic is the heart of TVD reconstruction: if the two one-sided slopes have opposite
signs the cell is an extremum, the limited slope is set to **zero** (donor cell there),
killing oscillations; otherwise the limited slope is the harmonic mean
`2·dql·dqr/(dql+dqr)/2`, which is small whenever either one-sided slope is small. The
wrapper `Reconstruct<plm, XNDIR>` applies this along whichever axis is being swept
(`X1DIR/X2DIR/X3DIR`), so the same routine serves all three dimensions. Reconstruction is
done on the **primitive** variables (smoother, physically bounded) rather than conserved.

## 7. Riemann solvers

Given the reconstructed left/right primitive states at a face, a **Riemann solver**
returns the single, upwinded interface flux. AthenaPK's MHD solvers are
`{ hlle, hlld, llf }`; the default and most accurate is **HLLD**.

The HLL family approximates the fan of waves emanating from the discontinuity by a small
number of constant states separated by the fastest signal speeds:

- **LLF / Rusanov** — one wave speed, maximally diffusive.
- **HLLE** — two waves (fastest left/right), no contact or Alfvén resolution.
- **HLLD** — *five* waves for MHD: fast-left, Alfvén-left, contact/entropy, Alfvén-right,
  fast-right, resolving the rotational discontinuities and contact that carry the
  field and density structure crucial to MHD turbulence and discs.

The implementation (`src/hydro/rsolvers/glmmhd_hlld.hpp`, after Miyoshi & Kusano 2005) is
a textbook five-state solver. Its skeleton:

```cpp
// Step 1: load L/R primitive states; decouple Bx and psi (the GLM pair, Mignone & Tzeferacos 2010 eq.24)
Real bxi  = 0.5*(wli[IB1]+wri[IB1]) - 0.5/c_h*(wri[IPS]-wli[IPS]);
Real psii = 0.5*(wli[IPS]+wri[IPS]) - 0.5*c_h*(wri[IB1]-wli[IB1]);
flxi[IB1] = psii;  flxi[IPS] = SQR(c_h)*bxi;

// Step 2: fast magnetosonic speeds -> outer wave speeds spd[0], spd[4]
const auto cfl = eos.FastMagnetosonicSpeed(wli[IDN], wli[IPR], wli[IB1], wli[IB2], wli[IB3]);
const auto cfr = eos.FastMagnetosonicSpeed(wri[IDN], wri[IPR], wri[IB1], wri[IB2], wri[IB3]);
spd[0] = min(wli[IV1]-cfl, wri[IV1]-cfr);
spd[4] = max(wli[IV1]+cfl, wri[IV1]+cfr);

// Step 3: L/R physical fluxes
// Step 4: contact speed S_M (spd[2]) and Alfven speeds (spd[1], spd[3])
spd[2] = (sdr*ur.mx - sdl*ul.mx + (ptl-ptr)) / (sdr*ur.d - sdl*ul.d);  // Miyoshi-Kusano eq.38
spd[1] = spd[2] - |bxi|/sqrt(ulst.d);    spd[3] = spd[2] + |bxi|/sqrt(urst.d);

// Step 5: the * and ** intermediate states
// Step 6: pick the flux for the region containing x/t = 0 (the face)
if      (spd[0] >= 0) flxi = fl;                 // supersonic from left
else if (spd[4] <= 0) flxi = fr;                 // supersonic from right
else if (spd[1] >= 0) flxi = fl + ulst;          // left *  state
else if (spd[3] <= 0) flxi = fr + urst;          // right * state
else if (spd[2] >= 0) flxi = fl + ulst + uldst;  // left ** state
else                  flxi = fr + urst + urdst;  // right ** state
```

The two GLM lines in Step 1 are the only modification beyond classic HLLD: the normal
field component `Bx` and the cleaning scalar `ψ` form a linear 2×2 system that advects at
the cleaning speed `c_h` and is solved exactly and separately from the seven MHD waves
(see §10).

## 8. The equation of state and wave speeds

AthenaPK uses an **adiabatic (ideal-gas) EOS**, P = (γ−1)e, with γ set per run (γ=5/3 for
atomic gas, γ=7/5=1.4 for the molecular collapse runs). The EOS object provides the
thermodynamic closure and the characteristic speeds the solvers need
(`src/eos/adiabatic_glmmhd.hpp`):

```cpp
Real SoundSpeed(const Real prim[NHYDRO]) const {
  return std::sqrt(gamma_ * prim[IPR] / prim[IDN]);          // c_s = sqrt(gamma P / rho)
}
Real FastMagnetosonicSpeed(Real d, Real p, Real bx, Real by, Real bz) const {
  Real asq = gamma_ * p;                                     // gamma P
  Real ct2 = by*by + bz*bz;                                  // transverse B^2
  Real qsq = bx*bx + ct2 + asq;
  Real tmp = bx*bx + ct2 - asq;
  return std::sqrt(0.5*(qsq + std::sqrt(tmp*tmp + 4.0*asq*ct2)) / d);
}
```

The fast magnetosonic speed
$c_f = \sqrt{\tfrac{1}{2}\big[(a^2+B^2/\rho) + \sqrt{(a^2+B^2/\rho)^2 - 4a^2 B_x^2/\rho}\big]}$
(with a²=γP/ρ) is the fastest MHD signal and therefore sets both the Riemann wave fan and
the hydrodynamic CFL timestep (§21). The EOS package also performs the
conserved→primitive inversion and applies density/pressure floors (`dfloor`, `pfloor`)
that keep the collapse robust as ρ spans many decades.

## 9. Time integration: the van-Leer predictor–corrector

The default integrator is **VL2** — the van-Leer two-stage, second-order predictor–
corrector used in Athena++. Schematically:

- **Stage 1 (predictor):** compute fluxes from Uⁿ, advance a half step to U^{n+1/2}.
- **Stage 2 (corrector):** recompute fluxes from U^{n+1/2}, advance the full step from Uⁿ.

In the task graph this appears as two passes (`stage = 1, 2`) through the same
`MakeTaskCollection`, each calling the registered flux function and then a Parthenon
update with the stage's Butcher-like coefficients:

```cpp
const auto flux_str = (stage == 1) ? "flux_first_stage" : "flux_other_stage";
FluxFun_t *calc_flux_fun = hydro_pkg->Param<FluxFun_t *>(flux_str);
auto calc_flux = tl.AddTask(none, calc_flux_fun, mu0);
...
auto update = tl.AddTask(set_flx, parthenon::Update::UpdateWithFluxDivergence<MeshData<Real>>,
                         mu0.get(), mu1.get(),
                         integrator->gam0[stage-1], integrator->gam1[stage-1],
                         integrator->beta[stage-1] * integrator->dt);
```

`mu0` ("base") holds the working state Uⁿ→U^{n+1}, `mu1` ("u1") holds the frozen Uⁿ that
the corrector blends back in via the `gam0/gam1` coefficients. Source terms are added
after the flux update on each stage (unsplit) or once per step (Strang-split); see §21.

## 10. Divergence control: GLM / Dedner cleaning

A profound difficulty of numerical MHD is maintaining **∇·B = 0**. Athena++ enforces it
exactly by constrained transport on a staggered (face) field. AthenaPK instead stores
**cell-centred** B and uses the **GLM (generalized Lagrange multiplier) / Dedner**
scheme: the constraint is *advected and damped* rather than enforced to round-off.

An auxiliary scalar field ψ (index `IPS`) is introduced and coupled to B so that
divergence errors propagate away at a finite "cleaning speed" `c_h` and are simultaneously
damped:

$$\partial_t \mathbf{B} + \nabla\psi = \dots,\qquad \partial_t\psi + c_h^2\,\nabla\cdot\mathbf{B} = -\frac{c_h^2}{c_p^2}\,\psi.$$

Two pieces implement this in AthenaPK:

**(i) The hyperbolic part** is built into the Riemann solver — the `bxi`/`psii`
decoupling lines in §7 advect the (Bx, ψ) pair at ±c_h.

**(ii) The cleaning speed** is set each cycle to the maximum the CFL allows
(`src/hydro/hydro.cpp`):

```cpp
hydro_pkg->UpdateParam("c_h", cfl_hyp * mindx / dt_hyp);   // fastest stable cleaning speed
```

i.e. `c_h = CFL · Δx_min / Δt_hyp`, the grid speed of light for the chosen timestep.

**(iii) The source/damping** is the Dedner source (`src/hydro/glmmhd/dedner_source.cpp`).
The parabolic damping of ψ is applied as an exponential factor (Mignone & Tzeferacos
2010), and the optional "extended" (non-conservative but more stable) Powell-like terms
feed divergence errors into the momentum and energy equations:

```cpp
const auto coeff = std::exp(-alpha * c_h * beta_dt / mindx);   // psi damping factor
...
if (extended) {
  const Real divB = 0.5*((prim(IB1,k,j,i+1)-prim(IB1,k,j,i-1))/dx1
                       + (prim(IB2,k,j+1,i)-prim(IB2,k,j-1,i))/dx2
                       + (prim(IB3,k+1,j,i)-prim(IB3,k-1,j,i))/dx3);
  cons(IM1,k,j,i) -= beta_dt * divB * prim(IB1,k,j,i);   // Powell momentum source
  cons(IM2,k,j,i) -= beta_dt * divB * prim(IB2,k,j,i);
  cons(IM3,k,j,i) -= beta_dt * divB * prim(IB3,k,j,i);
  cons(IEN,k,j,i) -= 0.5*beta_dt*( prim(IB1,k,j,i)*(prim(IPS,k,j,i+1)-prim(IPS,k,j,i-1))/dx1 + ... );
}
cons_pack(b, IPS, k, j, i) *= coeff;                     // damp psi
```

GLM is the right choice for AthenaPK because cell-centred B with GLM is far simpler to
implement on a performance-portable AMR framework than staggered constrained transport,
and the residual ∇·B is controlled to acceptably small values for collapse problems. (It
does mean that, unlike CT, ∇·B is not machine-zero — a point that matters for the Hall
term, §18.)

## 11. Adaptive mesh refinement and the Jeans criterion

Collapse spans >10⁵ in length scale, so uniform grids are hopeless; AMR is essential. The
governing rule for *where* to refine in self-gravitating collapse is the **Truelove
criterion** (1997): the local **Jeans length** must be resolved by at least some number
of cells, or artificial fragmentation results. AthenaPK implements this directly
(`src/refinement/jeans.cpp`).

The Jeans length in code units with 4πG = 1 is λ_J = 2π·v/√ρ, with v the relevant signal
speed: c_s for hydro, or (c_s + v_A) for MHD (the Athena++ convention; v_A = √(B²/ρ) in
Heaviside–Lorentz units). The number of cells per Jeans length is n_J = λ_J/Δx. The
kernel finds the *minimum* n_J on a block and tags it:

```cpp
const Real fac = 2.0 * M_PI / dx;
pmb->par_reduce("jeans refinement", kb.s,kb.e, jb.s,jb.e, ib.s,ib.e,
  KOKKOS_LAMBDA(int k, int j, int i, Real &lnjmin) {
    const Real rho = w(IDN,k,j,i);
    const Real cs  = std::sqrt(gam * w(IPR,k,j,i) / rho);
    Real v = cs;
    if (mhd) {
      const Real bsq = w(IB1,..)^2 + w(IB2,..)^2 + w(IB3,..)^2;
      v += std::sqrt(bsq / rho);                  // + Alfven speed (HL units, no 4pi)
    }
    const Real nj = v / std::sqrt(rho);
    lnjmin = std::min(lnjmin, nj);
  }, Kokkos::Min<Real>(njmin));
njmin *= fac;                                      // = 2 pi (c_s+v_A) / (dx sqrt(rho))

if (njmin < njeans)        return AmrTag::refine;       // under-resolved -> refine
if (njmin > 2.5 * njeans)  return AmrTag::derefine;     // over-resolved  -> coarsen
return AmrTag::same;
```

With `njeans` ≈ 16 the finest level keeps the local Jeans length covered by ≥16 cells
throughout the collapse, suppressing numerical fragmentation. *(A subtle units pitfall
lives here: B is stored in Heaviside–Lorentz form, so v_A = √(B²/ρ) with no 4π. Mixing
the Gaussian Alfvén speed B/√(4πρ) into a diagnostic underestimates n_J by √(4π) ≈ 3.5 —
a real trap when post-processing.)*

---

# Part III — The physics of star formation, as implemented

## 12. Self-gravity: the Poisson problem

Gravity is what drives the whole story, and it is the hardest term to compute well,
because it is **elliptic**: the potential everywhere depends on the density everywhere.
AthenaPK solves ∇²Φ = 4πG ρ each step with an **iterative linear solver on the AMR mesh**
(`src/self_gravity/`).

**The operator.** Discretizing the Laplacian in conservative flux form, the face flux of
Φ is −∇Φ and A·Φ is its divergence (`poisson_equation.hpp`):

```cpp
// flux = -grad(phi) at each face:
pack.flux(b, X1DIR, var_t(), k,j,i) =
    (pack(b,te,var_t(),k,j,i-1) - pack(b,te,var_t(),k,j,i)) / coords.Dxc<X1DIR>(k,j,i);
// A.phi = divergence of fluxes (FluxMultiplyMatrix):
out += (pack.flux(b,X1DIR,..,k,j,i)*Al - pack.flux(b,X1DIR,..,k,j,i+1)*Ar) / Vol;  // + y,z
```

This yields the standard 7-point Laplacian on a uniform grid, but written in
face-area/volume form so that it is **conservative across coarse–fine boundaries** (flux
correction is applied there). `SetDiagonal` supplies the matching diagonal for the
smoother.

**The right-hand side** is assembled every step from the gas density
(`FillPoissonRHS`), including the optional **Jeans swindle** (subtracting the mean density
so a periodic box has a well-posed potential):

```cpp
rhs = four_pi_G * (rho - grav_mean_rho);     // grav_mean_rho = 0 unless use_swindle
```

**The solver.** The linear system A·Φ = RHS is solved by a **BiCGSTAB Krylov iteration
preconditioned by geometric multigrid (GMG)** — the modern, scalable choice for elliptic
problems on AMR:

```cpp
using PoissEq        = PoissonEquation<grav::phi>;
using preconditioner_t = parthenon::solvers::MGSolver<PoissEq, prolongator_t>;
using SolverT          = parthenon::solvers::BiCGSTABSolver<PoissEq, preconditioner_t>;
```

Multigrid attacks all spatial wavelengths of the error by cycling between fine and coarse
levels (smoothing high-frequency error on fine grids, coarse grids handling the
long-wavelength error); BiCGSTAB wraps it to handle the non-symmetry introduced by the
AMR boundaries. The solve runs to a relative residual tolerance (e.g. 1e-10) each step.

**Applying gravity.** Once Φ is known, the source term updates momentum and energy
(`ApplyGravitySource`). Momentum gets ρ**g** = −ρ∇Φ (centred difference); energy gets the
**flux-weighted gravitational work** (Athena++/Artemis style), using the *hydro mass
flux* across each face so that the energy and momentum updates are mutually consistent:

```cpp
cons(IM1,k,j,i) += rho * hdtodx1 * (dpl1 + dpr1);          // rho * g * dt
...
Real de = hdtodx1 * (cons_flx.flux(X1DIR,IDN,k,j,i)   * dpl1 +
                     cons_flx.flux(X1DIR,IDN,k,j,i+1) * dpr1);   // mass-flux-weighted work
cons(IEN,k,j,i) += de;
```

Gravity is applied **once per full step, on the final integration stage** (§21), after
the MHD update, because the Poisson solve is by far the most expensive single operation
and need not be repeated per RK stage.

*Practical note baked into the code:* multi-rank GPU runs of the GMG-preconditioned solve
require `CUDA_LAUNCH_BLOCKING=1` to avoid an asynchronous race in the V-cycle ghost
exchange that otherwise sends Φ → NaN. This is documented in the solver source and is the
reason production GPU runs set that environment variable.

## 13. Driven turbulence

Molecular clouds are supersonically turbulent, and that turbulence both supports clouds
and seeds the density structure that collapses. AthenaPK can **drive** turbulence with a
stochastic, large-scale forcing implemented as a **few-mode Fourier acceleration field**
(`src/utils/few_modes_ft.cpp`), rather than a full FFT — far cheaper when only a handful
of large-scale modes are excited.

**Spectral content.** A set of wavevectors with |k| in a band around `k_peak` is chosen
(`MakeRandomModes`); their target amplitude follows a parabolic spectrum peaked at
k_peak:

```cpp
ampl = SQR(k_mag/k_peak) * (2.0 - SQR(k_mag/k_peak));   // peaks at k = k_peak, zero at 0 and sqrt2*k_peak
```

**Temporal evolution — an Ornstein–Uhlenbeck process.** Each mode's complex amplitude is
evolved as a correlated random walk with correlation time `t_corr`, so the forcing is
smooth in time (not white noise):

```cpp
const auto c_drift = std::exp(-dt / t_corr_);          // memory of previous state
const auto c_diff  = std::sqrt(1.0 - c_drift*c_drift);  // fresh random kick
var_hat(n,m) = var_hat(n,m)*c_drift + var_hat_new(n,m)*c_diff;
```

**Solenoidal vs compressive split.** A Helmholtz projection with weight `sol_weight`
controls how much of the forcing is divergence-free (solenoidal, sol_weight=1) versus
compressive (sol_weight=0):

```cpp
var_hat_new(0,m) = var_hat_new(0,m)*sol_weight + (1 - 2*sol_weight)*dot.real()*kx;  // + y,z
```

This matters physically: compressive forcing makes denser, more clustered structure than
solenoidal forcing at the same Mach number.

**Back to real space.** An explicit complex-to-real inverse transform over the few modes
builds the acceleration field on the grid; precomputed per-axis phase factors
(`phases_i/j/k`, computed in `SetPhases`, level-aware so AMR blocks get the right global
phase) make this a cheap sum:

```cpp
var_pack(b,n,k,j,i) += 2.*(var_hat(n,m).real()*phase.real() - var_hat(n,m).imag()*phase.imag());
```

The field is then applied as a momentum (and energy) source. In the collapse problem the
turbulence is typically used to set up the **initial velocity field** (a one-time
perturbation with a prescribed power spectrum and Mach number — e.g. mach 0.5, slope
−11/3) rather than continuously driven, seeding the angular momentum and substructure of
the collapsing core.

## 14. Thermodynamics: the barotropic EOS and the first hydrostatic core

A faithful collapse needs to capture the isothermal→adiabatic transition at the first
core (§1). The full solution is radiation transport, which AthenaPK does **not** have
(only optically-thin cooling, §19). Instead, the collapse problem generator
(`src/pgen/collapse_be.cpp`) uses the standard **barotropic-EOS workaround** — a
prescribed P(ρ) that mimics the thermal behaviour up to and through first-core formation.

The thermal energy follows the Masunaga–Inutsuka / Tomida form:

$$e_{\rm th} = \frac{\rho}{\gamma-1}\sqrt{1 + \left(\frac{\rho}{\rho_{\rm crit}}\right)^{2(\gamma-1)}}.$$

For ρ ≪ ρ_crit the square root → 1 and e_th = ρ/(γ−1), i.e. an **isothermal** gas with
c_s = 1 in code units (the 10 K cloud). For ρ ≫ ρ_crit the second term dominates and
e_th ∝ ρ^γ, i.e. an **adiabatic** gas that heats on compression — exactly the stiffening
that halts collapse and forms the pressure-supported first core. The transition density
ρ_crit ≈ 10⁻¹³ g cm⁻³ is the observed onset of optical thickness.

This is implemented as a **source term**, `collapse_be::ApplyBarotropicCooling`, wired
into the task graph immediately after the unsplit MHD sources and gated on the problem
having registered `collapse_be_rhocrit`:

```cpp
if (hydro_pkg->AllParams().hasKey("collapse_be_rhocrit")) {
  after_cooling = tl.AddTask(source_unsplit, collapse_be::ApplyBarotropicCooling, mu0.get(), tm,
                             integrator->beta[stage-1] * integrator->dt);
}
```

The initial condition is a **Bonnor–Ebert (BE) sphere** — the classic
pressure-truncated isothermal equilibrium — with profile ρ(r) ∝ (1 + r²/r_c²)^{−3/2},
mass and temperature set by the user (e.g. 6 M_⊙, 10 K), optionally seeded with the
turbulent velocity field of §13 and threaded by a uniform vertical magnetic field whose
strength is set by the mass-to-flux ratio μ. This is the setup that, with ambipolar
diffusion added (§17), forms the production runs documented separately.

## 15. The non-ideal MHD framework

In a 10 K molecular core the ionization fraction is minuscule (xᵢ ~ 10⁻⁷–10⁻⁸). The
magnetic field is tied to the *ions*, but the bulk mass is *neutral*; the imperfect
collisional coupling between them means ideal MHD fails. The corrections are the three
non-ideal terms, each a contribution to the induction-equation EMF:

| Effect | Diffusivity | EMF | Character | Dominant regime |
|---|---|---|---|---|
| Ohmic | η_O = const | η_O **J** | dissipative, isotropic | densest, best-shielded core |
| Ambipolar (AD) | η_A = Q_A B² | η_A **J**_⊥ | dissipative, anisotropic | intermediate densities, the disc-forming regime |
| Hall | η_H = Q_H B/ρ | η_H (**J**×**B**)/\|B\| | dispersive, handed | intermediate densities |

**Architecture.** AthenaPK adds non-ideal terms not via a staggered EMF (it has no
staggered field) but as **contributions to the cell-centred face fluxes** of the
conserved B components, plus the matching Poynting term in the energy flux. The mapping,
shared by all three effects, is read straight off the Ohmic template
(`OhmicDiffFluxIsoFixed`) and generalized: for an EMF vector **E** = (E₁,E₂,E₃) and
face-averaged **B** at a given face,

- **X1 face:** flux(IB2) += −E₃, flux(IB3) += +E₂, flux(IEN) += E₂B₃ − E₃B₂
- **X2 face:** flux(IB1) += +E₃, flux(IB3) += −E₁, flux(IEN) += E₃B₁ − E₁B₃
- **X3 face:** flux(IB1) += −E₂, flux(IB2) += +E₁, flux(IEN) += E₁B₂ − E₂B₁

The energy term is the Poynting flux **S** = **E**×**B**. For Ohmic (**E**=η**J**) this
reproduces the standard resistive flux and dissipates; for AD it dissipates (**E**_A·**J**
> 0); for Hall (**E**_H ⊥ **J**) it does **no net work** — exactly as physics demands.

All non-ideal fluxes are dispatched from `CalcDiffFluxes` (`diffusion.cpp`):

```cpp
const auto &resistivity = hydro_pkg->Param<Resistivity>("resistivity");
if (resistivity != Resistivity::none) OhmicDiffFluxIsoFixed(md);
const auto &ambipolar = hydro_pkg->Param<Ambipolar>("ambipolar");
if (ambipolar != Ambipolar::none)     AmbipolarDiffFluxIsoFixed(md);
const auto &hall = hydro_pkg->Param<Hall>("hall");
if (hall != Hall::none)               HallDiffFluxIsoFixed(md);
```

These are **parabolic** (AD, Ohmic) or **dispersive** (Hall) operators, which impose a
diffusive timestep limit Δt_diff ∝ Δx²/η that can be catastrophically small on refined
grids. AthenaPK offers two integration paths (`enum class DiffInt { none, unsplit, rkl2 }`):
**unsplit** (the diffusive term limits the global Δt directly) and **rkl2**
super-time-stepping (§20), which lets a single hyperbolic step take many cheap diffusive
sub-steps. Ambipolar and Ohmic are parabolic and may use either; Hall is dispersive and
**must** use unsplit (RKL2 is invalid for it).

## 16. Ohmic resistivity

The simplest non-ideal term: a constant scalar resistivity that dissipates field
isotropically, **E** = η_O **J**, η_O = const. It dominates only in the densest, most
shielded gas where even the ions decouple from the field. The diffusivity is a constant
(`resistivity.cpp`):

```cpp
Real OhmicDiffusivity::Get(const Real pres, const Real rho) const { return coeff_; }
```

`OhmicDiffFluxIsoFixed` computes the face-centred current **J** = ∇×**B** by finite
differences of the cell-centred field, forms **E** = η_O**J**, and applies the §15
flux↔EMF mapping. It serves as the structural template for the AD and Hall kernels. Its
timestep is the standard parabolic limit Δt = cfl_diff·fac·Δx²/η_O with
fac = 1/2, 1/4, 1/6 in 1D/2D/3D.

## 17. Ambipolar diffusion

**The physics.** Ambipolar diffusion (AD) is the dominant non-ideal effect over the
density range where discs form (n ~ 10⁸–10¹¹ cm⁻³). The field is frozen to the ions; the
neutrals drift across it under the ion–neutral drag. The result is a diffusion of the
field *relative to the neutral gas* that acts only **perpendicular** to **B** (drift along
field lines does nothing). AD is what relieves the magnetic braking catastrophe — by
letting flux slip out of the collapsing core, it removes the excess magnetic support and
torque and allows rotationally-supported discs to form. Its diffusivity grows steeply
with field strength, η_A = B²/(γ_AD ρ_i ρ) ∝ B², which is why it bites hardest exactly
where the field has been amplified by compression.

**The implementation** (`src/hydro/diffusion/ambipolar.cpp`) follows the
constant-coefficient model of Athena++:
$$\eta_A = Q_A\,B^2,\qquad \mathbf{E}_A = \eta_A\,\mathbf{J}_\perp = \eta_A\big(\mathbf{J} - (\mathbf{J}\cdot\hat{\mathbf b})\hat{\mathbf b}\big).$$

The diffusivity and the perpendicular-current EMF are:

```cpp
Real AmbipolarDiffusivity::Get(const Real bmag, const Real rho) const {
  return coeff_ * SQR(bmag);                                 // eta_A = Q_A * B^2
}
void PerpCurrentEMF(Real eta, Real j1,Real j2,Real j3, Real b1,Real b2,Real b3,
                    Real &e1,Real &e2,Real &e3) {
  const Real bsq   = b1*b1 + b2*b2 + b3*b3 + TINY_NUMBER;
  const Real jdotb = j1*b1 + j2*b2 + j3*b3;
  e1 = eta * (j1 - jdotb*b1/bsq);                            // J_perp = J - (J.bhat) bhat
  e2 = eta * (j2 - jdotb*b2/bsq);
  e3 = eta * (j3 - jdotb*b3/bsq);
}
```

Each of the three directional kernels (X1/X2/X3) computes the **full** face-centred
current — all three components j₁,j₂,j₃ via finite differences of the cell-centred field
(the transverse derivatives use 4-point averages to land on the face) — because the
perpendicular projection mixes all components. It then face-averages **B** and ρ, forms
η_A and **E**_A, and applies the mapping. The X1 kernel, abbreviated:

```cpp
// j2 = d3B1 - d1B3, j3 = d1B2 - d2B1, j1 = d2B3 - d3B2   (face-centred curl B at i-1/2)
const Real b1 = 0.5*(prim(IB1,k,j,i-1)+prim(IB1,k,j,i));    // face-averaged B
const Real b2 = 0.5*(prim(IB2,k,j,i-1)+prim(IB2,k,j,i));
const Real b3 = 0.5*(prim(IB3,k,j,i-1)+prim(IB3,k,j,i));
const Real bmag = std::sqrt(SQR(b1)+SQR(b2)+SQR(b3));
const Real eta  = ad_diff.Get(bmag, rho);
Real e1,e2,e3; PerpCurrentEMF(eta, j1,j2,j3, b1,b2,b3, e1,e2,e3);
cons.flux(X1DIR, IB2, k,j,i) += -e3;
cons.flux(X1DIR, IB3, k,j,i) +=  e2;
cons.flux(X1DIR, IEN, k,j,i) += e2*b3 - e3*b2;              // Poynting flux
```

**The timestep.** AD is parabolic, so the bare stable step is Δx²/η_A with the
dimensional factor:

```cpp
Real fac = 0.5; if (ndim==2) fac = 0.25; else if (ndim==3) fac = 1.0/6.0;
...
const auto eta = ad_diff.Get(bmag, prim(IDN,k,j,i));
min_dt = fmin(min_dt, SQR(coords.Dxc<1>(k,j,i)) / (eta + TINY_NUMBER));   // + y,z
return cfl_diff * fac * min_dt_ad;
```

Because η_A ∝ B² and Δt ∝ Δx²/η_A, this step **collapses** as the field amplifies and the
mesh refines toward the core — the central practical challenge of AD collapse runs, which
is why they are run with RKL2 super-time-stepping (§20).

**Validation.** AD is verified quantitatively in AthenaPK via a damped-Alfvén
eigenmode (`diffusion` pgen, `iprob=50`): the measured By decay rate matches the exact
root of s² + η_A k² s + k² v_A² = 0 to 0.04% at 128 cells with clean second-order
convergence, and the η_A ∝ B² signature is confirmed. A 1-D oblique **C-shock**
(`cshock` pgen) reproduces the continuous ion–neutral shock structure AD is famous for.

## 18. The Hall effect

**The physics.** The Hall effect is the strangest of the three. It is **dispersive, not
dissipative**: it does no net work, but it makes the induction equation depend on the
*sign* of **B** relative to the rotation/current — it introduces a handedness. Physically
it arises because, at intermediate densities, the ions begin to decouple from the field
while the electrons remain tied to it, so the current-carrying species drifts. Its
dramatic astrophysical consequence is **bimodal disc formation**: depending on whether the
magnetic field and the rotation axis are aligned or anti-aligned, the Hall term either
spins up or spins down the inner core, producing large discs in one case and almost none
in the other from otherwise identical initial conditions. Its diffusivity scales as
η_H = Q_H B/ρ.

**The implementation** (`src/hydro/diffusion/hall.cpp`). The diffusivity and EMF:
$$\eta_H = Q_H\,\frac{B}{\rho},\qquad \mathbf{E}_H = \eta_H\,\frac{\mathbf{J}\times\mathbf{B}}{|\mathbf B|} + \eta_{\rm floor}\,\mathbf{J}.$$

```cpp
Real HallDiffusivity::Get(const Real bmag, const Real rho) const {
  return coeff_ * bmag / (rho + TINY_NUMBER);                // eta_H = Q_H * B / rho
}
void HallEMF(Real eta_h, Real eta_floor, Real bmag, Real j1,j2,j3, Real b1,b2,b3,
             Real &e1,Real &e2,Real &e3) {
  const Real jxb1 = j2*b3 - j3*b2;                           // (J x B)
  const Real inv_b = 1.0/(bmag + TINY_NUMBER);
  e1 = eta_h * jxb1 * inv_b + eta_floor * j1;                // Hall EMF + optional Ohmic floor
  // e2, e3 analogously
}
```

Two implementation points are notable:

1. **Beyond Athena++.** Stock Athena++ *does not actually evolve* the Hall term — its
   `HallEMF` is commented out and only Ohmic+AD EMFs are applied. The Hall EMF in
   AthenaPK is a genuinely working implementation, going beyond the parent code.
2. **The Ohmic floor.** On a **cell-centred** grid (no constrained-transport
   divergence-free field), the Hall term is numerically unstable without a little
   dissipation. A small `hall_ohmic_floor_code` adds an η_floor**J** Ohmic component that
   stabilizes it; runs should always use the floor or a companion Ohmic resistivity.
3. **Integration path.** Being dispersive, Hall is **incompatible with RKL2**; AthenaPK
   guards against this and forces the unsplit path. Its timestep is the whistler limit
   Δt ∝ Δx²/η_H — note the Δx² (not Δx) scaling means Hall, like AD, becomes very
   expensive at high resolution; that is physical, not a bug.

**Validation.** Hall is verified via a circularly-polarized **whistler** eigenmode
(`diffusion` pgen, `iprob=60`): the measured phase advance matches the combined
Hall+Ohmic-floor dispersion relation to 0.4% on the whistler branch, and the two
helicities split cleanly (fast whistler vs slow ion-cyclotron) — the unambiguous Hall R/L
signature absent from AD and Ohmic.

## 19. Optically-thin cooling

For problems where the gas is optically thin (the cluster/ISM problems, not the dense
collapse core), AthenaPK provides **tabular cooling**
(`src/hydro/srcterms/tabular_cooling.cpp`): a user-supplied cooling function Λ(T) is read
from a table and applied as an energy source term, with sub-cycling to handle stiff
cooling. For the dense collapse runs this is **not** the relevant thermodynamics — there
the barotropic EOS of §14 stands in for the (absent) radiation transport. The honest
statement of AthenaPK's thermodynamic capability for star formation is therefore: *no
radiative transfer; optically-thin tabular cooling, or a barotropic-EOS workaround for the
first-core regime.*

## 20. Super-time-stepping (RKL2)

The parabolic timestep Δt_diff ∝ Δx²/η of AD and Ohmic diffusion can be thousands of
times smaller than the hyperbolic (MHD) step on a refined grid. Sub-cycling every cell at
Δt_diff would be ruinous. **Runge–Kutta–Legendre super-time-stepping (RKL2)** solves this:
it takes a *sequence of s carefully-chosen sub-steps* whose composite is stable over a
step far larger than s·Δt_diff — specifically O(s²)·Δt_diff — while remaining
second-order accurate. It is the reason AthenaPK can integrate AD through collapse at all.

**The stage count.** Given the hyperbolic step τ (here ½Δt, because diffusion is
Strang-split), the number of RKL stages is chosen from the ratio of hyperbolic to
diffusive step (Meyer+ 2014):

```cpp
int s_rkl = static_cast<int>(0.5*(std::sqrt(9.0 + 16.0*tau/mindt_diff) - 1.0)) + 1;
if (s_rkl % 2 == 0) s_rkl += 1;            // ensure odd number of stages
```

**The recursion.** The Legendre coefficients are built per stage (Meyer+ 2012/2014) and
the state is advanced by the two-register recursion (`hydro_driver.cpp`):

```cpp
// first stage:  Y_1 = Y_0 + mu_tilde_1 * tau * M(Y_0)
Yjm1 = Y0 + mu_tilde_1 * tau * MY0;
// stage j:      Y_j = mu_j Y_{j-1} + nu_j Y_{j-2} + (1-mu_j-nu_j) Y_0
//                     + mu_tilde_j tau M(Y_{j-1}) + gamma_tilde_j tau M(Y_0)
const Real Yj = mu_j*Yjm1 + nu_j*Yjm2 + (1.0-mu_j-nu_j)*Y0
              + mu_tilde_j*tau*MYjm1 + gamma_tilde_j*tau*MY0;
```

where M(·) is the diffusive flux divergence (`CalcDiffFluxes` → `FluxDivHelper`). Each
stage exchanges ghost cells so the diffusive operator is consistent across blocks.

**The safety cap.** Even RKL2 has a practical limit: stretching the step too far (very
large s) loses accuracy. AthenaPK caps the ratio with `rkl2_max_dt_ratio`
(`src/hydro/hydro.cpp`):

```cpp
} else if (hydro_pkg->Param<DiffInt>("diffint") == DiffInt::rkl2) {
  const auto max_dt_ratio = hydro_pkg->Param<Real>("rkl2_max_dt_ratio");
  if (max_dt_ratio > 0.0 && dt_hyp / dt_diff > max_dt_ratio) {
    min_dt = std::min(min_dt, max_dt_ratio * dt_diff);     // don't let STS stretch beyond this
  }
}
```

So the global step is min(Δt_hyp, max_dt_ratio·Δt_diff). When AD becomes the binding
constraint near the core, this cap directly controls cost vs. fidelity: too small a cap
(e.g. 100) throttles the run far below the hydro limit; raising it (e.g. 1000) recovers
most of the speed, until the *bare* Δt_diff itself collapses as B → large and Δx → small.

---

# Part IV — Synthesis

## 21. Anatomy of a single timestep

We can now read the whole conductor's score — `HydroDriver::MakeTaskCollection` — as a
sequence. For each cycle, for each integration stage (VL2 has two):

1. **(Stage 1 only) Operator-split sources.** If RKL2 diffusion is active, the
   super-time-stepping sub-cycle for a ½Δt diffusive update is added first
   (`AddSTSTasks(..., 0.5*tm.dt)`); then any Strang-split sources get their first ½Δt.
2. **Register init.** Copy Uⁿ into the `u1` register for the corrector.
3. **Flux calculation.** `calc_flux_fun` reconstructs (§6) and solves Riemann problems
   (§7) for all faces, writing the hyperbolic MHD fluxes. Optional first-order flux
   correction guards against new extrema.
4. **Flux correction at coarse/fine boundaries.** Send/receive/set corrected fluxes so the
   scheme stays conservative across AMR levels.
5. **Update.** `UpdateWithFluxDivergence` advances the conserved variables by −∇·F with
   the stage coefficients (§9).
6. **Unsplit sources.** `AddUnsplitSources` adds cooling, the GLM/Dedner source (§10), and
   problem-defined terms; then (if a collapse run) `ApplyBarotropicCooling` (§14).
7. **(Final stage) Strang-split sources** get their closing ½Δt.
8. **Boundary exchange.** Ghost cells are filled (local + non-local), prolongated/
   restricted across levels, and physical BCs applied.
9. **(Final stage) Self-gravity.** `SelfGravity::SolvePoisson` runs the BiCGSTAB+MG solve
   (§12); `ApplyGravitySource` adds ρ**g** and the flux-weighted work.
10. **FillDerived.** Recompute primitives from conserved via the EOS.
11. **(Final stage) Second RKL2 ½Δt** diffusive update (Strang symmetry), then reset of the
    reduction scratch params (`mindx`, `dt_hyp`, `dt_diff`).
12. **(Final stage) New timestep.** `EstimateTimestep` reduces over all cells: the
    hyperbolic CFL Δt_hyp = CFL·Δx/(|v|+c_f), the diffusive Δt_diff, and the RKL2 cap
    combine into the next Δt; the GLM cleaning speed c_h is refreshed from it.
13. **(Final stage) Tracers, then AMR tagging.** Lagrangian tracer particles are advected;
    `Refinement::Tag` evaluates the Jeans criterion (§11) and flags blocks to refine or
    derefine before the mesh is remade.

The dependency edges between these tasks let Parthenon overlap MPI ghost exchange with
interior computation, which is what makes the scheme scale.

## 22. A collapse simulation end-to-end

Putting the pieces together, a magnetized core-collapse run proceeds like this:

1. **Initial condition** (`collapse_be` pgen): a Bonnor–Ebert sphere (mass, T, BE
   truncation set by the user), threaded by a uniform field of prescribed mass-to-flux
   ratio, perturbed by a turbulent velocity field of chosen Mach number and spectrum
   (§13). Units fixed by 4πG=1, c_s=1 (§23).
2. **Isothermal collapse.** Self-gravity (§12) overwhelms thermal + magnetic + turbulent
   support; the central density runs away. The barotropic EOS (§14) keeps the gas at
   ~10 K. The Jeans criterion (§11) drives ever-deeper AMR refinement to keep λ_J resolved.
3. **Flux dragging and AD.** The field is compressed and amplified with the gas (B ∝
   ρ^{~2/3}, flux freezing). Ambipolar diffusion (§17) lets neutrals slip across the
   field, regulating magnetic braking. As B grows, η_A ∝ B² grows and the AD timestep
   shrinks; RKL2 (§20) keeps the run feasible, with `rkl2_max_dt_ratio` trading cost
   against fidelity.
4. **First-core formation.** At ρ ≥ ρ_crit ≈ 10⁻¹³ g cm⁻³ the barotropic EOS stiffens to
   adiabatic; pressure halts the collapse and a pressure-supported first hydrostatic core
   forms (~AU scale, ~10⁻² M_⊙). The timestep, now set by the AD diffusion on the
   refined core, becomes the limiting cost.
5. **Analysis.** Mass-to-flux ratio, plasma β, flux-freezing exponent, angular-momentum
   transport, and disc properties are measured from the snapshots — the science output.

## 23. Units

Self-gravity in AthenaPK runs in code units with **4πG = 1** (set as `four_pi_G` in the
`self_gravity` block), and the collapse problem additionally normalizes c_s = 1 at 10 K
and the BE central density to 1. Magnetic fields are in **Heaviside–Lorentz** form: the
magnetic pressure is B²/2 (no 4π), and the Alfvén speed is v_A = √(B²/ρ). This has two
consequences worth flagging because they bite in analysis:

- Plasma β in code units is 2P/B² (not the Gaussian 8πP/B²); the two differ by 4π.
- The Alfvén speed and the Jeans criterion (§11) use √(B²/ρ); using the Gaussian
  B/√(4πρ) underestimates v_A by √(4π) ≈ 3.5.

Conversion factors to cgs (`code_length_cgs`, `code_density_cgs`, `code_bfield_cgs =
v₀√ρ₀`, etc.) are written to a `units.json` alongside the outputs for post-processing.

## 24. Performance portability: the GPU story

The reason AthenaPK exists is the last line of every kernel: it is a Kokkos `par_for` that
compiles to CPU threads or GPU kernels from one source. The same binary logic runs on a
workstation and on a GPU supercomputer. For star formation this is double-edged:

- The **hyperbolic MHD** work — reconstruction, Riemann solves, flux updates — is
  flop-dense, regular, and embarrassingly parallel: it maps beautifully onto GPUs.
- The **elliptic self-gravity** solve — multigrid V-cycles with frequent small ghost
  exchanges across coarse–fine and rank boundaries — is latency- and
  communication-bound, and GPUs do *not* dominate here.
- The **AMR** machinery (regridding, prolongation/restriction, load balance) is likewise
  bookkeeping-heavy.

The net effect on a deep-AMR self-gravitating collapse is that the GPU advantage over CPU
is modest compared to the raw-FLOP ratio — empirically, for this class of problem, a
single high-end data-centre GPU is worth on the order of ~10–15 CPU cores, not the ~100×
one might naively expect, precisely because the gravity solve and AMR dilute the
flop-dense MHD work the GPU is good at. This is not a deficiency of AthenaPK but an
intrinsic property of the algorithm mix in self-gravitating collapse, and it is the right
context in which to read any GPU-vs-CPU timing comparison.

---

## Coda

AthenaPK assembles, on a single performance-portable framework, every ingredient the
star-formation story requires: a robust high-resolution Godunov MHD core (PLM + HLLD +
VL2), GLM divergence control, Jeans-driven AMR, a scalable multigrid-preconditioned
self-gravity solver, stochastic turbulence forcing, a barotropic surrogate for the
first-core thermodynamics, and — the focus of the current development — a complete
non-ideal MHD suite (Ohmic, a faithfully ported and validated ambipolar diffusion, and a
genuinely working Hall term that goes beyond the parent Athena++), integrated either
explicitly or via RKL2 super-time-stepping. The one missing physical ingredient,
radiative transfer, is stood in for by the barotropic EOS up to first-core formation. The
result is a code that can follow a turbulent, magnetized molecular cloud core from its
quiescent equilibrium, through runaway gravitational collapse, to the birth of a
magnetically-regulated first hydrostatic core — on hardware ranging from a laptop to a
GPU supercomputer.

---

*Source references (all paths relative to the AthenaPK root):*
`src/main.hpp` (variables/enums) · `src/hydro/hydro.cpp` (package, timestep, c_h) ·
`src/hydro/hydro_driver.cpp` (task collection, VL2, RKL2) ·
`src/recon/plm_simple.hpp` (reconstruction) ·
`src/hydro/rsolvers/glmmhd_hlld.hpp` (HLLD) · `src/eos/adiabatic_glmmhd.hpp` (EOS) ·
`src/hydro/glmmhd/dedner_source.cpp` (GLM cleaning) ·
`src/refinement/jeans.cpp` (AMR criterion) ·
`src/self_gravity/{self_gravity,poisson_equation,self_gravity_driver}.*` (gravity) ·
`src/utils/few_modes_ft.cpp` (turbulence) ·
`src/hydro/diffusion/{resistivity,ambipolar,hall}.cpp` (non-ideal MHD) ·
`src/pgen/collapse_be.cpp` (BE sphere + barotropic EOS).
```
