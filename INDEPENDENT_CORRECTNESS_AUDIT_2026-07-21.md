# AthenaPK Independent Correctness Audit

Date: 2026-07-21  
Scope: AthenaPK application code under `src/`, build configuration, and collapse production inputs.  
Method: static, read-only, from-scratch audit. Existing comments, validation claims, run products, and documentation were not accepted as evidence.

## Instructions for the next auditor

Treat this report as a set of claims to reproduce, not as ground truth. Re-open every cited location, trace the relevant data flow, and independently confirm or reject each finding. Preserve the distinction between defects proven directly from code and hypotheses requiring runtime experiments. Do not modify code while auditing unless explicitly asked in a later task.

## Executive conclusion

The inspected code should not yet be treated as quantitatively validated for fossil-field or protostellar-collapse conclusions. The production configuration contains inconsistent ionization parameters across coupled packages, first-order/lagged self-gravity coupling, and inconsistent radiation temperature normalization. Separately, the collapse unit sidecar converts code magnetic field to gauss incorrectly by a factor of `sqrt(4*pi)`.

The ideal-MHD flux convention itself is internally Heaviside-Lorentz (`Pmag = B^2/2`). The Poisson operator and force signs are internally consistent. The configured GPU build genuinely enables CUDA and Hopper (`Kokkos_ARCH_HOPPER90=ON`). These clean observations do not validate the full coupled calculation.

## Findings, ranked by severity

### 1. HIGH — physics — chemistry and non-ideal MHD use inconsistent cosmic-ray ionization rates

Status: **CONFIRMED**

Evidence:

- `src/chemistry/chemistry.cpp:74`
- `src/hydro/diffusion/ionization.hpp:82`
- `runs/prod_t4_full/fhc.in:88`
- `runs/prod_t4_full/fhc.in:91`

The production chemistry network defaults to `zeta = 1e-16 s^-1`, while the equilibrium ionization model used by production Ohmic and Hall diffusion defaults to `1e-17 s^-1`. Production ambipolar diffusion consumes chemistry-evolved `x_e`, while Ohmic and Hall use the separate equilibrium model. The three non-ideal terms therefore do not describe the same charge population. In a gas-phase equilibrium limit, this alone produces an approximately `sqrt(10)` mismatch in electron abundance and resistivity scale before grain effects.

Failure scenario: ambipolar flux loss responds to gas ionized at ten times the cosmic-ray rate used for Hall rotation and Ohmic dissipation, silently changing relative non-ideal effects and fossil-flux retention.

Suggested fix: define a single shared ionization parameter, make all consumers use it, and reject initialization if chemistry and diffusion request inconsistent values.

### 2. HIGH — numerics/physics — self-gravity is a lagged first-order end-of-step source

Status: **CONFIRMED**

Evidence:

- `src/hydro/hydro_driver.cpp:553`
- `src/hydro/hydro_driver.cpp:597`
- `src/hydro/hydro_driver.cpp:606`
- `src/self_gravity/self_gravity.cpp:431`

The VL2 hydro flux update completes before Poisson is solved. Gravity is then applied once on the final stage using the full timestep. It is absent from predictor and corrector flux states. The first step likewise has no gravity in its predictor.

Failure scenario: central collapse carries an `O(dt)` splitting error in momentum and gravitational work; the mass flux used for gravity-energy work was computed without the contemporaneous force. Collapse time, kinetic energy, central entropy, and hydrostatic balance can be biased.

Suggested fix: use stage-consistent gravity or a time-centered kick-drift-kick scheme, with potential and mass flux at compatible temporal states.

### 3. HIGH — physics/units — collapse output converts code magnetic field to gauss incorrectly

Status: **CONFIRMED**

Evidence:

- `src/pgen/collapse_be.cpp:122`
- `src/units.hpp:92`
- `src/hydro/rsolvers/glmmhd_hlld.hpp:99`
- `src/hydro/diffusion/ionization.hpp:78`

The MHD kernels use Heaviside-Lorentz units:

`Pmag = B_code^2 / 2`

so the physical conversion is:

`B_cgs = B_code * sqrt(4*pi*rho0) * v0`.

The collapse problem instead writes `sqrt(rho0)*v0` to `units.json`, omitting `sqrt(4*pi)`. Any physical field inferred from that sidecar is too small by `sqrt(4*pi)`, i.e. it is reported as `0.282094...` of the correct value. This affects interpretation rather than the evolved code state.

Suggested fix: use `sqrt(4*pi*rho0)*v0`, preferably through the common `Units` implementation.

### 4. HIGH — physics — production radiation uses a temperature unit inconsistent with EOS, chemistry, and non-ideal MHD

Status: **CONFIRMED**

Evidence:

- `runs/prod_t4_full/fhc.in:109`
- `src/radiation/radiation.cpp:61-68`
- `src/chemistry/chemistry.cpp:50`
- `src/hydro/diffusion/ionization.hpp:73-78`

Radiation identifies `mu=2.29` as the normalization consistent with `v0=1.9e4 cm/s`, giving `T0 ~= 10.015 K`. Production overrides radiation to `mu=2.33`, giving about `10.19 K`, while chemistry and diffusion retain `10.015 K`. This shifts `a_R*T0^4/e0` by about seven percent and also shifts opacity temperatures and matter-radiation equilibrium.

Suggested fix: remove the override or replace duplicated package scales with one shared immutable unit system checked during initialization.

### 5. HIGH — code/physics — sink accretion removes magnetic energy without accreting magnetic field

Status: **CONFIRMED**. Sinks are disabled in the inspected tier-4 production input.

Evidence:

- `src/sinks/sinks.cpp:557-562`

Accretion multiplies total energy by `(1-phi)` but leaves `B` unchanged. Primitive recovery subsequently subtracts the unchanged `B^2/2`, so thermal energy is spuriously reduced by `phi*B^2/2`. A magnetically dominated accretion cell can acquire negative internal energy and be repaired only by a pressure floor.

Suggested fix: retain magnetic energy when retaining magnetic field, or implement a defined magnetic-flux accretion model and transfer the associated energy consistently.

### 6. HIGH — numerics — radiation accepts a failed implicit matter-coupling iteration silently

Status: missing failure handling is **CONFIRMED**; occurrence in production is **SUSPECTED**.

Evidence:

- `src/radiation/radiation_moments.cpp:333-365`

After `inner_iteration_max`, the code writes the last iterate regardless of residual. There is no convergence flag, global reduction, warning, timestep rejection, or fallback. Floors also change the system without an explicit energy-defect diagnostic.

Failure scenario: optically thick cells around dust sublimation or H2 dissociation can accept a non-solution and violate the intended gas-plus-reduced-radiation energy balance.

Confirming experiment: record final residuals, iteration counts, and floor activations per cell over representative production checkpoints; require every cell to meet tolerance.

Suggested fix: reject or subcycle failed cells, expose maximum residuals, and account explicitly for floor energy.

### 7. MEDIUM — numerics/code — collapse cooling and final split sources have conflicting task dependencies

Status: dependency error is **CONFIRMED**; an observed race is **SUSPECTED**.

Evidence:

- `src/hydro/hydro_driver.cpp:567-580`

The cooling task depends on `source_unsplit`, but on the final stage the Strang source also depends directly on `source_unsplit` rather than `after_cooling`. Both can modify `cons`. With radiation active, collapse cooling still changes ambient momentum.

Suggested fix: make the final Strang task depend on `after_cooling` and preserve a single state-mutating dependency chain.

Confirming experiment: vary task concurrency or insert a diagnostic delay and compare ambient momentum and total energy.

### 8. MEDIUM — physics — production gravity boundary is finite-box zero Dirichlet, not isolated

Status: quantitative impact is **SUSPECTED**.

Evidence:

- `runs/prod_t4_full/fhc.in:126-131`
- `src/self_gravity/self_gravity.cpp:39-45`
- `src/self_gravity/multipole.hpp:46-62`

Production imposes `phi=0` on all six box faces. This is not `phi -> 0` at infinity. The force for a perfectly centered monopole is mostly insensitive to the missing constant potential, but turbulence, rotation, asymmetry, and outflow generate multipoles affected by the finite cubic wall.

Confirming experiment: compare successively larger zero-boundary boxes and the available multipole boundary using central acceleration, virial terms, collapse time, and `-GM(<r)/r^2` force profiles.

### 9. MEDIUM — EOS/code — out-of-table states silently clamp to edge values

Status: **CONFIRMED**

Evidence:

- `src/eos/eos_table.hpp:37-56`
- `src/eos/eos_table.hpp:73-83`

The interpolation clamps every coordinate to the table edge without reporting an out-of-range state. The inspected binary spans approximately:

- `log10(rho_code) = [-1.738, 18.15]`
- `log10(esp_code) = [0.0704, 5.333]`
- `log10(T/K) = [0.9031, 5.477]`

Nominal second-core density is within the density domain, but shocked or floored cells can leave the energy or temperature domain. Pressure, sound speed, or heat capacity then become artificially constant along the clamped direction.

Suggested fix: count and report edge hits, fail production runs on them, and extend the table or supply thermodynamically consistent asymptotic continuations.

### 10. MEDIUM — physics/code — dust evolution uses an ideal-gas temperature under the hydrogen EOS

Status: **CONFIRMED**. Dust is disabled in the inspected production input.

Evidence:

- `src/dust/dust.cpp:110-119`

`ReactDust` computes `T=(gamma-1)*e/rho` even when `hydro/eos=hydrogen`. This is wrong in dissociation and ionization regimes and consequently moves grain growth/sublimation and radiation opacity transitions.

Suggested fix: use the shared EOS table temperature or reject this configuration.

### 11. MEDIUM — MPI/AMR/restart — tracer diagnostic values lack explicit restart metadata

Status: **SUSPECTED**

Evidence:

- `src/tracers/tracers.cpp:57-68`
- contrast with `src/sinks/sinks.cpp:41-45`

The tracer swarm has restart metadata, while tracer `rho`, pressure, velocity, and magnetic values are only marked `Real`. Sink values explicitly add `Restart`. A write/read round trip is needed to determine whether tracer values are dropped or stale immediately after restart.

### 12. LOW — MPI/restart — dynamically assigned sink ID state is not demonstrably restart-persistent

Status: **SUSPECTED**

Evidence:

- `src/sinks/sinks.cpp:78-79`
- `src/sinks/sinks.cpp:399-429`

`next_sink_id` is a mutable package parameter initialized from `nseed`; no application field stores it. If mutable package parameters are not serialized, newly created sinks can reuse restored IDs. Confirm by creating sinks, restarting, creating another, and checking uniqueness. A robust implementation derives the next ID from the global maximum restored ID.

## Clean observations supported by code read

1. The ideal-MHD core consistently uses Heaviside-Lorentz magnetic units. Magnetic pressure and energy are `B^2/2`, and no explicit `4*pi` belongs in MHD flux terms.
2. The Poisson implementation solves `laplacian(phi) = 4*pi*G*rho`. Its potential wells are negative and `-grad(phi)` points toward overdensities; no gravity sign reversal was found.
3. Hall is excluded from RKL2 super-time-stepping and is assigned a whistler timestep. Ohmic/AD and the optional Hall Ohmic floor use dimension-dependent parabolic factors.
4. Main hydro and non-ideal updates are expressed as face fluxes and conservative divergence updates.
5. GPU-shared sink accretion accumulators use Kokkos atomics, and global contributions use MPI reductions.
6. The inspected GPU CMake cache is Release, has CUDA enabled, and sets `Kokkos_ARCH_HOPPER90=ON`. The CPU cache is Release and enables Kokkos OpenMP. No fast-math flag was found in either cache.
7. Bitwise CPU/GPU or rank-count determinism should not be assumed because atomics, reductions, AMR ordering, and MPI reduction trees alter floating-point order.

## Coverage summary

| Subsystem | What was examined | Depth / gaps |
|---|---|---|
| Driver and `main.cpp` | Package registration, VL2/STS task graph, gravity, radiation, chemistry, dust, sinks, tracers, AMR tagging | Deep task-order read; scheduler race not dynamically instrumented |
| `hydro/` | State registration, primitive/conserved conversion, production HLLD, source dispatch, timestep bookkeeping | Moderate/deep; alternative Riemann solvers not read end-to-end |
| `hydro/glmmhd` | GLM Riemann coupling and Dedner source | Moderate; no runtime divergence measurement |
| `hydro/diffusion` | Ohmic, AD, Hall dispatch, EMFs, eta cache, fused timestep, RKL2 staging | Deep on production path; not every stencil independently re-derived |
| `self_gravity/` | RHS, operator sign, AMR flux correction, boundaries, multipoles, force and work | Deep; temporal coupling is the main confirmed defect |
| `radiation/` | M1 closure, HLL transport, causal clamp, CFL, opacity units, matter coupling, RSLA | Deep static read; no runtime convergence experiment |
| `chemistry/` | Scalar layout, reduced network, positivity, H/C caps, EOS temperature | Deep for production network; it is not the full GOW17 network |
| `eos/` | Binary domain, loader, interpolation, inverse pressure, sound speed, RT inverse | Moderate/deep; table was not independently regenerated |
| `sinks/` | Swarm metadata, creation, gravity, accretion, atomics, N-body, migration | Deep; disabled in inspected production input |
| `dust/` | Reaction kernel and opacity coupling | Moderate; disabled in inspected production input |
| `tracers/` | Registration, seeding, RK2 advection, migration, sampling | Moderate; incompatible with adaptive meshes by explicit check |
| Reconstruction | Production PLM selection and ghost requirement | Shallow; PPM/WENO variants not fully audited |
| Refinement | Jeans criterion and task placement | Moderate; no dynamic hysteresis test |
| Boundaries | Application and gravity boundary registration | Shallow/moderate; vector reflection parity not exhaustive |
| Collapse generator | BE normalization, ICs, units, turbulence, chemistry/radiation seeds, ambient control | Deep |
| Build configuration | CPU/GPU caches, architecture, backend, optimization flags | Deep enough to verify configured targets |
| Restart/AMR metadata | Major fields and swarms | Moderate; no bitwise round-trip experiment |

## Full-step data-flow trace

For the inspected tier-4 configuration:

1. `collapse_be` initializes cell-centered conserved GLM-MHD state, five scalar densities, LTE radiation energy, zero radiation flux, turbulent velocity, and a uniform code magnetic field.
2. Primitive variables are recovered with the hydrogen EOS table.
3. VL2 reconstructs hydro states and HLLD computes ideal-MHD/GLM face fluxes. In mixed diffusion mode, Hall is added unsplit while Ohmic/AD are advanced in Strang-split RKL2 halves.
4. Hydro flux divergence and unsplit sources update conserved gas state.
5. Collapse cooling is invoked, but with radiation present it does not overwrite thermal energy; it still zeroes momentum outside the initial sphere.
6. On the final VL2 stage, split sources run, then Poisson is solved from the final density and a full-step gravity source is applied.
7. M1 radiation transport subcycles at the reduced-light CFL, then implicit matter coupling changes radiation, gas thermal energy, momentum, and kinetic energy.
8. Chemistry reacts scalar abundances using gas temperature. Production AD subsequently consumes chemistry `x_e`; production Ohmic/Hall use a separate equilibrium charge calculation.
9. Derived primitives are refreshed, the second RKL2 half runs, timestep estimates are reduced, and AMR tagging occurs.
10. Independent radiation fields and hydro fields are eligible for checkpointing; sink values are explicitly restart-marked, while tracer value-level restart behavior remains uncertain.

## Recommended priority order

1. Unify ionization parameters and charge-state ownership across chemistry, Ohmic, Hall, and AD.
2. Correct magnetic-field physical conversion and re-evaluate every fossil-field value previously derived from `units.json`.
3. Unify temperature and unit calibration across EOS, radiation, chemistry, and diffusion.
4. Make gravity time-centered with the hydro integrator.
5. Add mandatory radiation nonlinear convergence diagnostics/rejection.
6. Fix the cooling task dependency.
7. Validate isolated gravity boundaries and perform box-size convergence.
8. Add EOS edge-hit diagnostics.
9. Repair sink MHD energy accounting before enabling sinks in magnetized production.
10. Perform restart round trips for all fields and swarms on CPU, one GPU, and multi-GPU MPI.

## Runtime experiments needed before scientific use

- Magnetized collapse convergence in spatial resolution, AMR level, CFL, GLM parameters, diffusion caps, and RKL2 freeze behavior.
- Collapse comparisons with stage-centered versus current gravity.
- Multipole versus enlarged-box gravity boundary convergence.
- Radiation equilibrium, free-streaming, diffusion-limit, radiative shock, and optically thick collapse tests with recorded implicit residuals.
- Ohmic, AD, and signed Hall wave/diffusion tests using the same ionization state as production chemistry.
- EOS thermodynamic derivative and entropy consistency checks across H2 dissociation and both ionization regimes.
- Exact restart round trips before/after AMR and particle migration.
- CPU/GPU statistical parity and conservation trends; bitwise identity is not a reasonable requirement for reduction-heavy paths.

## Detailed remediation and convergence plan

This section is the minimum program needed to turn the present implementation into a controlled numerical experiment. It does not make the physical model complete; it establishes that the equations actually selected are solved consistently and that reported fossil-field quantities are not dominated by coding defects, resolution, boundaries, or arbitrary numerical controls.

### Governing principles and publication gates

1. Freeze scientific production while any confirmed high-severity defect remains open.
2. Give every change a small analytic/unit test, a subsystem benchmark, and a coupled-collapse regression. A successful collapse alone is not validation.
3. Keep correctness changes separate from model changes. Never combine, for example, a gravity-integration fix and a new opacity in one comparison.
4. Record all runtime parameters, compile flags, dependency revisions, GPU model, MPI layout, and physical tables in machine-readable provenance.
5. Pre-register quantitative acceptance criteria before examining the preferred result.
6. Use dimensionless diagnostics wherever possible and attach units only through one authoritative unit system.
7. Treat caps, floors, reduced light speed, density floors, pressure floors, and refinement ceilings as modifications of the solved problem. Report where and how often each activates.
8. No production result passes if its central conclusion changes materially under one additional level of spatial, temporal, microphysical, or boundary convergence.

### Phase 0 — establish a reproducible baseline

#### 0.1 Capture the exact baseline

- Record source revision and all submodule revisions.
- Archive `CMakeCache.txt`, `compile_commands.json`, compiler and CUDA versions, MPI version, HDF5 version, Kokkos configuration, and the exact executable hash.
- Save the fully expanded runtime input after command-line and restart overrides.
- Record rank count, GPUs per node, block size, pack count, load-balancing settings, environment variables, and deterministic seeds.
- Generate a manifest containing hashes of EOS, opacity, chemistry, and any turbulence data.

#### 0.2 Define primary science observables before fixes

At fixed density thresholds and Lagrangian enclosed masses, measure:

- central density and temperature versus time;
- first-core and second-core formation times;
- core mass, radius, entropy, and accretion rate;
- maximum and mass-weighted magnetic field;
- magnetic flux `Phi(M) = integral B dot dA` through surfaces enclosing fixed mass;
- mass-to-flux ratio `M/Phi`;
- effective scaling `d ln B / d ln rho`;
- magnetic, kinetic, thermal, radiation, and gravitational energies;
- Ohmic, Hall, and ambipolar diffusivity distributions;
- Elsasser numbers, magnetic Reynolds numbers, Hall length, and ambipolar length;
- signed angular momentum and magnetic torque by enclosed mass;
- volume and mass norms of `div(B)`;
- total mass, momentum, angular momentum, gas-plus-radiation energy, and elemental budgets;
- floor/cap activation counts and the mass/energy/flux affected by them.

These definitions must be implemented once in a versioned analysis package. Diagnostics must use the corrected magnetic conversion, not the existing sidecar value.

#### 0.3 Build small restartable test configurations

Create cheap CPU and one-GPU inputs for:

- uniform advection;
- MHD linear waves;
- GLM divergence pulse;
- Ohmic diffusion;
- ambipolar C-shock;
- Hall whistler wave with both Hall signs;
- Jeans-stable and Jeans-unstable gravity;
- isolated gravitating sphere;
- radiation free streaming, shadowing, equilibrium, diffusion, and radiative shock;
- EOS pressure/temperature round trips;
- chemistry box reactions;
- a reduced-resolution collapse reaching first-core formation.

Each test must emit machine-readable norms and pass/fail thresholds, not image comparisons.

### Phase 1 — repair unit and shared-parameter inconsistencies

#### 1.1 Replace duplicated unit calibrations

- Introduce one immutable `PhysicalUnits` package created before all physics packages.
- Store base mass, length, time, temperature, density, energy-density, magnetic-field, opacity, and diffusivity units.
- Derive `B_unit = sqrt(4*pi*rho_unit)*v_unit` for the Heaviside-Lorentz code field.
- Require radiation, chemistry, EOS, diffusion, sinks, dust, problem generators, and output analysis to consume this package.
- Eliminate independent default values such as `rho_unit`, `T_unit`, `v_unit`, and `length_unit` from subsystem input blocks, or retain them only as checked aliases that must agree within a specified tolerance.
- Write the authoritative unit set into every restart and output file, not only a sidecar.

Acceptance tests:

- Round-trip randomly sampled code quantities to CGS and back at machine precision.
- Verify `B^2/2` maps to `B_cgs^2/(8*pi)` for at least 1,000 random states.
- Confirm that radiation LTE, EOS temperature, chemistry rates, and ionization rates all see the same physical temperature for identical code states.
- Reprocess historical fossil-field outputs with both conversions and explicitly publish the correction factor where relevant.

#### 1.2 Unify charge-state and cosmic-ray parameters

- Create one `IonizationEnvironment` containing cosmic-ray rate, attenuation prescription, radionuclide floor, composition, grain population, and thermal ionization parameters.
- Make chemistry and every non-ideal coefficient consume the same environment.
- Decide explicitly whether Ohmic, Hall, and AD use equilibrium charges or time-dependent chemistry. Do not mix them without a derived, documented hybrid closure.
- If chemistry evolves only a subset of charged species, compute the full conductivity tensor from a charge-neutral reconstruction of electrons, ions, and charged grains consistent with those evolved abundances.
- Add initialization checks for incompatible chemistry species, scalar indices, EOS, dust model, or ionization parameters.

Acceptance tests:

- Homogeneous boxes at fixed `(rho,T,B)` must return identical electron abundance and conductivity inputs from chemistry and the diffusion evaluator.
- Compare conductivities against an independent high-precision reference implementation over a logarithmic grid spanning the full collapse domain.
- Require relative errors below 1% away from sign changes and report absolute error around `eta_H = 0`.
- Verify exact charge neutrality and elemental conservation to roundoff after every chemistry step.

### Phase 2 — correct time integration and task dependencies

#### 2.1 Make gravity stage-consistent

Candidate implementation sequence:

1. Store `phi^n` and compute gravity at the beginning of the step.
2. Apply a half-step momentum/energy kick or incorporate gravity in the VL2 predictor.
3. Advance hydro to the appropriate predicted state.
4. Solve Poisson for the predicted/final density.
5. Complete the second half kick/corrector with temporally centered acceleration.
6. Derive gravitational work using mass fluxes and potentials at compatible time centering.
7. Preserve AMR refluxing and boundary communication between every state that uses neighbor potential.

Test hierarchy:

- Uniform acceleration: exact velocity and kinetic-energy change.
- Pressureless homologous sphere: compare collapse trajectory to a high-accuracy ODE reference.
- Stable Jeans mode: frequency and phase convergence.
- Unstable Jeans mode: exponential growth-rate convergence.
- Hydrostatic polytrope/BE sphere: measure spurious velocity for many sound-crossing times.
- Isolated two-body gas/sink test if sinks are enabled.

Acceptance criteria:

- Demonstrate second-order temporal convergence before shocks for smooth gravity tests.
- Hydrostatic test maximum Mach number must decrease at the formal order with timestep and resolution.
- Global energy error, including gravitational energy, must converge and remain below a preselected fraction of the physical energy change.

#### 2.2 Repair task dependencies

- Make collapse cooling precede final Strang and first-order split sources through a single dependency chain.
- Audit every task that writes `cons`, `prim`, radiation fields, particle swarms, diffusivity caches, and package reduction parameters.
- Draw a generated task dependency graph for every enabled package combination.
- Add debug metadata declaring read/write sets for tasks, then detect unordered writes to the same field.
- Test with alternate task scheduling, forced delays, and multiple pack counts.

Acceptance criteria:

- Thread/task sanitizer or the nearest available instrumentation reports no unordered host writes.
- Forced scheduling perturbations produce roundoff-level differences only, not changes in conservation or physical trajectories.

#### 2.3 Validate RKL2 and mixed Hall splitting

- Re-derive stage coefficients and the `tau/dt_diff` stage-count expression in a design note with executable symbolic checks.
- Test both frozen and stage-refreshed diffusivities against a high-accuracy explicit reference.
- Quantify the model error introduced by frozen, state-dependent diffusivities.
- Keep Hall entirely outside parabolic STS unless a mathematically valid Hall integrator is introduced.
- Verify ghost exchange and primitive refresh before every stencil evaluation.

Acceptance criteria:

- Constant-coefficient diffusion reaches the expected second-order temporal behavior of RKL2.
- State-dependent tests demonstrate convergence as the hydro timestep is reduced.
- `rkl2_freeze_eta=true` may be used only if the difference from refreshed coefficients is below the declared science tolerance for flux retention.

### Phase 3 — repair local robustness defects

#### 3.1 Radiation nonlinear solve

- Return per-cell convergence state, iteration count, final residuals, and floor activation.
- Reduce maximum residual and number of failures globally every step.
- On failure, retry with local or global source subcycling; reject the hydro step if convergence remains unsuccessful.
- Use a bracketed or safeguarded nonlinear method across opacity discontinuities and EOS dissociation zones.
- Include opacity derivatives consistently or adopt a robust Jacobian-free/bracketed solve.
- Track gas-plus-RSLA-radiation energy before and after coupling, including explicit floor injection.

Acceptance criteria:

- Zero unconverged cells in production.
- Equilibrium boxes preserve LTE to solver tolerance for at least 1,000 coupling times.
- Stiff heating/cooling boxes converge to the same state from both hotter and colder initial conditions.
- Coupling results converge as source substep size and nonlinear tolerance are tightened.

#### 3.2 EOS bounds and thermodynamic consistency

- Add device-side counters for low/high hits on every table axis.
- Abort production on the first unapproved edge hit.
- Extend the table beyond all observed states with a safety margin.
- Verify monotonic pressure with energy, positive `c_s^2`, positive heat capacity, Maxwell consistency where applicable, and smooth derivatives.
- Replace independent forward/inverse interpolation surfaces with a thermodynamically consistent free-energy formulation if discrepancies are significant.

Acceptance criteria:

- No edge hits in the final production ensemble.
- `rho,e -> P -> e` and `rho,T -> e -> T` round trips satisfy a preselected tolerance over the entire table.
- Entropy along adiabatic trajectories remains constant to interpolation accuracy.
- Shock-tube solutions converge without table-induced oscillations.

#### 3.3 Sink MHD energy and restart state

- Correct accretion energy partition so retained magnetic field retains its magnetic energy.
- Define whether mass leaves magnetic flux behind, accretes it, or transfers it through a subgrid reconnection/decoupling prescription.
- Recompute `next_sink_id` from the global maximum restored ID.
- Store all sink integrator and feedback state in restart data.
- Make creation/accretion decisions deterministic under MPI ownership changes where practical.

Acceptance criteria:

- Single-cell and spherical accretion tests conserve mass, momentum, angular momentum, and the chosen energy budget.
- Magnetized accretion never creates negative thermal energy from bookkeeping.
- Restart immediately before creation/accretion reproduces IDs and integrated quantities within the expected reduction-order tolerance.

#### 3.4 Dust and tracer fixes

- Obtain dust temperature from the selected EOS rather than the ideal-gas identity.
- Mark tracer values as restartable or explicitly derived and force recomputation before any consumer.
- Add AMR-aware particle migration/prolongation support before using tracers in adaptive collapse.

### Phase 4 — subsystem verification matrix

#### 4.1 Ideal MHD and GLM

Run at least 32, 64, 128, and 256 cells per wavelength where affordable:

- all seven MHD linear-wave families at multiple propagation angles;
- Brio-Wu and multidimensional shock tubes;
- circularly polarized Alfven wave;
- field-loop advection;
- Orszag-Tang vortex;
- rotor and blast tests;
- a controlled divergence pulse.

Measure L1/L2/Linf errors, conservation, wave speed, damping, phase error, positivity interventions, and `div(B)` norms. Repeat on CPU and GPU.

Acceptance criteria:

- Expected formal convergence on smooth waves.
- No resolution-independent magnetic-flux drift in field-loop advection.
- GLM-cleaning parameter variations do not change collapse flux retention within the declared tolerance.

#### 4.2 Non-ideal MHD

- Ohmic decay of sinusoidal fields with analytic exponential rate.
- Ambipolar diffusion of a Gaussian/sinusoidal field.
- Steady C-shock profiles over weak and strong coupling regimes.
- Hall whistler dispersion and polarization for both signs of `eta_H`.
- Oblique multidimensional tests to expose tensor/index sign errors.
- Combined Ohmic+AD+Hall comparisons against an independent 1D reference solver.

Sweep resolution, CFL, coefficient magnitude, STS ratio, eta-cache mode, and AMR placement. Confirm dissipative heating/energy flux where applicable.

#### 4.3 Gravity and AMR

- Point-mass and uniform-sphere potentials with zero, enlarged-box, and multipole boundaries.
- Off-center and quadrupolar mass distributions.
- Repeat with the density peak crossing block, rank, and refinement boundaries.
- Measure solver residual separately from physical discretization error.
- Confirm refinement/restriction does not create self-force or discontinuous acceleration.

#### 4.4 Radiation

- Uniform LTE equilibrium and matter-radiation relaxation.
- Free-streaming pulse and beam.
- Crossing beams, explicitly documenting the known M1 limitation.
- Shadow casting.
- Optically thick diffusion with analytic diffusion coefficient.
- Static and moving radiative shocks.
- Radiation pressure acceleration.
- Tests across every Bell-Lin regime boundary.
- RSLA sweep against full-`c` reference in reduced problems.

#### 4.5 Chemistry, ionization, and EOS

- Reaction-box comparisons with a high-order stiff ODE solver.
- Positivity and elemental/charge conservation under extreme timesteps.
- Equilibrium approach from multiple initial compositions.
- Conductivity tables cross-checked against an independent microphysics calculation.
- EOS comparisons with published Saumon-Chabrier-van Horn or equivalent hydrogen/helium tables in their shared domain.

### Phase 5 — AMR, MPI, GPU, and restart demonstrations

#### 5.1 AMR convergence

- Repeat smooth waves and diffusion with a stationary refinement interface and with the feature crossing it.
- Compare uniform-grid results to AMR results having the same finest resolution.
- Track mass, momentum, energy, scalar, radiation, and magnetic flux across coarse-fine boundaries.
- Test refine/derefine cycling to expose accumulated interpolation error.
- Confirm every evolved field has appropriate conservative restriction/prolongation; document why nonconservative variables such as potential use different operators.

#### 5.2 Restart round trips

For hydro, GLM scalar, radiation, chemistry scalars, dust, gravity potential/initial guess, sinks, tracers, RNG state, and mutable counters:

1. Write a restart.
2. Read it without advancing and compare every stored value.
3. Advance continuously for `N` steps and compare with a run restarted at `N/2`.
4. Repeat across a remesh and particle migration.
5. Repeat with one and multiple MPI ranks and GPUs.

Bitwise equality is required for the no-advance write/read state. Continued trajectories may use norm-based tolerances if collective ordering changes.

#### 5.3 GPU correctness

- Run compute-sanitizer memory, race, and initialization checks on reduced tests.
- Compile a debug GPU build with bounds and finite-value checks.
- Add optional per-kernel NaN/Inf detection at package boundaries.
- Compare CPU and GPU convergence curves, not only final values.
- Verify all host/device transfers explicitly for table loads, particle gathers, and package views.

### Phase 6 — coupled collapse convergence campaign

#### 6.1 Define a reference ladder

Use geometrically increasing effective resolutions, for example at least four levels in finest cell size. Keep all physical parameters fixed. At each level ensure the Jeans length, Hall length, resistive length, ambipolar length, first-core pressure scale height, radiation mean-free path where transport requires it, and relevant disk scale height are measured against cell size.

Do not call a run converged merely because the Jeans criterion is satisfied. Report the number of cells resolving every controlling physical length.

#### 6.2 Independent convergence axes

Vary one axis at a time:

- root-grid and maximum AMR resolution;
- block size and refinement interface placement;
- Jeans refinement threshold;
- CFL and integrator timestep;
- Poisson residual tolerance and solver choice;
- radiation CFL, reduced-light factor, nonlinear tolerance, and subcycling;
- GLM wave speed and damping;
- RKL2 maximum ratio and frozen/refreshed eta;
- Ohmic/Hall/AD caps and floors;
- domain size and gravity boundary condition;
- density/pressure/radiation floors;
- EOS table resolution and bounds;
- opacity model and resolution;
- chemistry tolerance/substeps and network choice;
- turbulence mode count, seed, and spectral realization.

#### 6.3 Convergence measures

Compare runs at common physical landmarks, not only equal wall time:

- when central density crosses fixed thresholds;
- at first-core formation;
- immediately before and after second collapse;
- at fixed enclosed protostellar mass;
- at fixed time after second-core formation.

For each landmark compare radial and angular distributions of density, temperature, entropy, ionization, diffusivities, velocity, angular momentum, and magnetic field. Compute relative differences in `Phi(M)`, `M/Phi`, core field, and magnetic-energy spectra.

Proposed minimum numerical acceptance criterion:

- primary fossil-flux observables change by less than 5% between the two highest resolutions;
- formation times and core masses/radii change by less than 5%;
- global conservation errors are at least an order of magnitude below the physical differences being interpreted;
- no cap/floor affects more than a predeclared negligible fraction of the mass or magnetic flux used in the conclusion;
- boundary, RSLA, and GLM variations change primary observables by less than the quoted systematic uncertainty.

If 5% cannot be achieved, report the measured convergence order and extrapolate with an uncertainty large enough to cover the unresolved trend. Do not silently widen acceptance after seeing results.

#### 6.4 Turbulent ensemble

One seed is not a physical prediction. At the final two feasible resolutions run enough independent turbulent realizations to estimate intrinsic scatter. A reasonable starting design is 8–16 seeds, expanded until confidence intervals on the median flux-retention statistic stabilize. Use paired seeds across physics variants to reduce variance.

Report distributions, not a preferred realization. Separate numerical convergence uncertainty, microphysical uncertainty, and turbulent sample variance.

### Phase 7 — release gate for the corrected current model

The corrected present model may be used scientifically only after:

- all confirmed high findings are fixed and reviewed;
- all suspected high/medium findings have decisive experiments;
- unit, subsystem, AMR, restart, and GPU test matrices pass automatically;
- nonlinear radiation failures and EOS edge hits are zero;
- isolated gravity and RSLA convergence are demonstrated;
- primary fossil-field observables meet the predeclared convergence threshold;
- a turbulent ensemble quantifies realization scatter;
- all cap/floor activations are disclosed and their effect bounded;
- an independent reviewer reproduces selected tests from archived inputs and containers.

## State-of-the-art physics roadmap for a definitive fossil-field study

No simulation can be literally perfect, and “most reliable ever” cannot be guaranteed by adding modules. The defensible target is the most complete, transparent, and uncertainty-quantified calculation feasible: multiple independently verified numerical formulations, microphysics tied to laboratory/observational constraints, and conclusions stable under a declared hierarchy of approximations. This section is intentionally more ambitious than the remediation plan above.

### Program architecture: use a hierarchy, not one monolithic run

Build four linked model tiers:

1. **Microphysics tier:** zero-dimensional and one-dimensional high-accuracy EOS, chemistry, grain charging, conductivity, and opacity calculations.
2. **Core-collapse tier:** global three-dimensional radiation-non-ideal-MHD from molecular-cloud core to second core.
3. **Zoom tier:** nested ultra-high-resolution calculations around the first/second core resolving magnetic decoupling and reconnection scales as far as feasible.
4. **Long-term tier:** protostellar evolution, disk, outflow, accretion, and surface-field evolution beyond second-core birth.

Pass tables and uncertainty distributions upward from the microphysics tier. Cross-check overlapping epochs between tiers. Avoid forcing one calculation to span every scale with uncontrolled subgrid assumptions.

### Workstream A — thermodynamics and composition

#### A.1 Adopt a thermodynamically consistent free-energy EOS

- Use an H/He mixture with realistic mass fractions and metals.
- Include translational, rotational, and vibrational H2 states, ortho/para H2, H2 dissociation, H and He ionization stages, electron degeneracy where relevant, Coulomb corrections, and radiation contributions.
- Derive pressure, entropy, internal energy, heat capacities, adiabatic derivatives, and sound speed from one Helmholtz free energy so Maxwell relations hold.
- Cover molecular-cloud through stellar-interior regimes with smooth table joins.
- Include composition as a state variable where non-equilibrium chemistry matters.
- Quantify interpolation error and table uncertainty against established EOS products.

#### A.2 Non-equilibrium thermal chemistry

- Evolve H, H2, H+, He species, electrons, key metal ions, molecular ions, CO/C/C+, O-bearing species, and grain charge populations where their timescales compete with collapse.
- Include cosmic-ray ionization, attenuation with column density, radionuclide ionization, thermal ionization of alkalis and hydrogen, recombination, grain adsorption/desorption, molecular formation/destruction, and relevant three-body reactions at high density.
- Couple reaction enthalpies to gas energy and radiation consistently.
- Use a robust implicit stiff solver with analytic/sparse Jacobians, positivity, exact elemental constraints, and error control.
- Generate reduced networks through sensitivity analysis, and verify them against a larger reference network throughout representative collapse trajectories.

Deliverable: a published open microphysics table/code with uncertainty bands for EOS and charged-species abundances across `(rho,T,B,column,grain distribution)`.

### Workstream B — grains and opacity

#### B.1 Evolve a grain-size and charge distribution

- Use multiple dust fluids or moments/bins for size distribution rather than one characteristic radius.
- Include coagulation, fragmentation, sublimation, condensation, ice mantles, charging, and possibly drift/settling.
- Evolve charge states per bin or solve a charge distribution consistent with the chemistry.
- Couple grains to neutral gas, charged fluids, magnetic fields, and radiation.
- Track how grain surface area changes recombination and how charged grains alter Hall-sign reversals.

#### B.2 Frequency-dependent opacity

- Construct composition- and size-dependent absorption and scattering opacities.
- Include gas opacities after dust sublimation, molecular bands, H-minus, bound-free/free-free, electron scattering, and line contributions where important.
- Use distinct Planck, Rosseland, flux, and energy means derived from the same monochromatic data.
- Propagate uncertainties from grain composition and size distribution.
- Validate against benchmark opacity databases and laboratory constraints.

### Workstream C — radiation transport

#### C.1 Move beyond grey M1

- Implement multigroup radiation hydrodynamics spanning dust infrared through stellar/ionizing bands.
- Include mixed-frame terms consistently through the required order in `v/c`: radiation advection, work, Doppler shifts, aberration, and momentum exchange.
- Treat scattering, including anisotropy where relevant.
- Preserve total gas-radiation energy and momentum to solver tolerance.
- Use implicit or IMEX transport/source integration suitable for optically thick cells.

#### C.2 Hybrid transport

M1 cannot represent crossing beams and can produce artificial radiation interactions. For maximum credibility:

- combine moment transport for diffuse radiation with ray tracing or Monte Carlo for direct protostellar irradiation;
- compare selected snapshots against a high-accuracy Monte Carlo post-processing/reference solution;
- use an asymptotic-preserving method that recovers the correct diffusion limit without resolving the mean free path;
- verify reduced-speed-of-light validity locally using diffusion, dynamical, and light-crossing timescale ratios.

### Workstream D — magnetic induction and charged-fluid physics

#### D.1 Constrained transport

- Implement staggered constrained transport so discrete `div(B)=0` is preserved to roundoff.
- Use AMR-compatible divergence-preserving prolongation, restriction, and reflux-curl.
- Retain GLM only as an optional comparison, not the primary production method.
- Compare CT and GLM flux-retention measurements to quantify historical numerical divergence transport.

#### D.2 Conductivity tensor from consistent microphysics

- Compute Ohmic, Hall, and Pedersen conductivities from all relevant charged species and grain bins.
- Include sign changes naturally rather than clipping them.
- Couple conductivity to the same time-dependent chemistry and grain model.
- Include cosmic-ray attenuation and thermal-ionization transitions.
- Propagate microphysical rate uncertainties into diffusivity distributions.

#### D.3 Test single-fluid validity and add multifluid capability where required

- Evaluate collision frequencies, drift speeds, gyrofrequencies, charged-species inertia, and strong-coupling assumptions throughout the run.
- Where assumptions fail, evolve ions/electrons/charged grains as separate fluids or use a generalized multifluid Ohm law.
- Include ion-neutral and grain-neutral drag heating.
- Resolve or accurately model ion-neutral drift and Hall scales.
- Consider electron pressure/battery terms if they influence field generation/topology.

#### D.4 Reconnection and unresolved transport

- Measure numerical resistivity with controlled tests at every AMR level.
- Resolve current sheets where possible using refinement based on current density, field curvature, and non-ideal length scales.
- Use physically motivated explicit dissipation rather than relying on grid-scale reconnection.
- If turbulent reconnection remains unresolved, bracket it with multiple justified subgrid models and treat spread as systematic uncertainty.

### Workstream E — gravity, rotation, and dynamics

- Use stage-centered self-gravity with an isolated multipole or Green-function boundary validated against analytic solutions.
- Include gas, dust, and particle mass consistently in Poisson.
- Use angular-momentum-conserving formulations and report torque budgets by source: gravity, Maxwell stress, Reynolds stress, radiation, outflow, and numerical transport.
- Resolve first-core, disk, pseudodisk, second-core, and launching regions with physics-based AMR criteria, not density/Jeans length alone.
- Refine on Jeans length, pressure scale height, radiation gradients, current sheets, Hall/AD/resistive lengths, disk scale height, and vorticity.

General relativity is unnecessary for first-core and newborn low-mass second-core conditions, but the Newtonian validity range should be stated. If evolution proceeds toward compact/high-mass regimes, reassess this assumption.

### Workstream F — initial and environmental conditions

#### F.1 Replace a single idealized sphere with a controlled ensemble

- Retain the BE sphere as a clean reference problem.
- Add cores extracted from larger molecular-cloud simulations with self-consistent turbulence, magnetic topology, accretion environment, and filamentary structure.
- Sample core mass, rotation, magnetization, field-core misalignment, turbulence spectrum/mode mixture, metallicity, cosmic-ray environment, and external pressure.
- Use observational priors rather than arbitrary parameter ranges.

#### F.2 Magnetic topology

- Explore uniform, turbulent, hourglass, and reversally structured initial fields.
- Quantify net versus unsigned flux and topology-dependent reconnection.
- Report results versus initial mass-to-flux ratio and field alignment.

### Workstream G — first-core to protostar and beyond

Fossil-field conclusions cannot stop exactly at second-core creation if the claimed observable is a stellar field.

- Continue through second-core accretion with sufficient surface resolution.
- Couple to a stellar evolution model that supplies radius, internal structure, convection, luminosity, and accretion energy deposition.
- Model accretion shocks and their radiative efficiency.
- Include magnetically launched outflows/jets when the launching region is resolved or use a validated subgrid model.
- Include protostellar irradiation and, when relevant, ionizing radiation.
- Track magnetic flux crossing the stellar surface, remaining in the disk, expelled by outflows, and dissipated through non-ideal/reconnection processes.
- Distinguish inherited large-scale flux from small-scale dynamo-generated field.
- Follow convective stability and dynamo action long enough to determine whether the birth field survives, reorganizes, or is overwhelmed.

### Workstream H — sinks and subgrid models

For a highest-credibility birth-field study, avoid sinks until the physical stellar surface is resolved through the epoch used for the central conclusion. If sinks are later unavoidable:

- define a control-volume conservation law for mass, momentum, angular momentum, total energy, magnetic flux, and radiation;
- calibrate accretion and magnetic-flux handling against resolved zoom simulations;
- carry stellar magnetic multipoles and spin;
- couple accretion luminosity and feedback;
- perform sink-radius convergence over several radii;
- never interpret sink-contained field without a validated flux-transfer prescription.

### Workstream I — numerical methods for extreme scale separation

- Use high-order, positivity-preserving Godunov schemes with robust fallback near shocks.
- Use CT-compatible Riemann/EMF construction.
- Adopt IMEX or fully implicit integration for stiff radiation, chemistry, and non-ideal terms where explicit/STS splitting produces excessive error.
- Use local time stepping or subcycling only with conservative synchronization across AMR.
- Develop error estimators tied to science observables, particularly enclosed magnetic flux.
- Maintain double precision for conserved state, gravity reductions, EOS inversions, and global budgets; evaluate mixed precision only after error bounds.
- Add reproducible reductions for selected diagnostics even if the full evolution is not bitwise deterministic.

### Workstream J — verification by independent formulations

The strongest publication should not depend on one code.

- Reproduce a common idealized collapse with at least one independent radiation non-ideal-MHD code.
- Cross-check microphysics tables with an independently written solver.
- Compare CT versus GLM, multigroup versus grey, equilibrium versus time-dependent chemistry, and single-fluid versus multifluid approximations on overlapping problems.
- Conduct a blinded comparison where analysis definitions are fixed before code labels are revealed.
- Publish discrepancies, not only agreements.

### Workstream K — uncertainty quantification and experimental design

#### K.1 Separate uncertainty classes

Quantify separately:

- discretization and timestep error;
- AMR and boundary error;
- radiation closure/RSLA error;
- microphysical rate and opacity uncertainty;
- initial-condition/turbulent variance;
- subgrid/sink uncertainty;
- analysis and unit-conversion uncertainty.

#### K.2 Ensemble design

- Use Latin hypercube, sparse-grid, or Bayesian experimental design for expensive parameters.
- Use paired random seeds when comparing physical models.
- Construct emulators only after numerical convergence is established.
- Report posterior/predictive distributions for flux retention, not a single curve.
- Perform global sensitivity analysis to identify which uncertain inputs dominate the fossil field.

#### K.3 Primary causal diagnostics

Track a Lagrangian magnetic-flux budget:

- ideal advection;
- Ohmic contribution;
- ambipolar contribution;
- Hall redistribution;
- resolved reconnection/topological change;
- AMR/numerical transport estimate;
- flux expelled by outflows;
- flux crossing the protostellar surface.

The Hall term redistributes magnetic flux/geometry but is nondissipative in the ideal continuum limit; diagnostics must not mislabel Hall work as irreversible flux destruction.

### Workstream L — observational connection

- Generate synthetic polarized dust emission with evolved grains and radiative alignment assumptions.
- Generate Zeeman observables where applicable.
- Predict core/disk field morphology, outflow polarization, rotation measures, and protostellar surface-field distributions.
- Forward-model instrumental resolution and selection effects.
- Compare to ensembles of observed cores/protostars rather than one object.
- Clearly distinguish directly observable quantities from simulation-internal constructs such as flux through a chosen isodensity surface.

### Workstream M — open and reproducible publication package

Release:

- source and exact dependency revisions;
- containers or environment recipes;
- all initial conditions and parameter manifests;
- EOS/opacity/conductivity tables with generation code;
- raw reduced diagnostics and representative full checkpoints;
- analysis software and figure scripts;
- convergence data, failed tests, and negative results;
- machine-readable uncertainty budget;
- a documented procedure for reproducing at least one reference run on another system.

Arrange an independent replication before submission. Archive immutable artifacts with persistent identifiers.

### Proposed staged publication sequence

1. **Methods paper:** verified CT radiation non-ideal-MHD, microphysics, unit system, restart/AMR tests, and public benchmark suite.
2. **Controlled BE-sphere paper:** resolution/boundary/RSLA/microphysics convergence and a turbulent ensemble through second-core birth.
3. **Environmental ensemble paper:** cores drawn from cloud simulations with observational priors and uncertainty quantification.
4. **Protostellar survival paper:** evolution of inherited flux through accretion, convection, outflows, and dynamo competition.

### Final readiness criteria for the flagship fossil-field claim

A flagship claim should be made only if:

- two independent codes agree within the combined numerical uncertainty on a shared benchmark;
- primary flux-retention observables are spatially and temporally converged or have a credible extrapolation;
- CT and an alternate divergence treatment demonstrate that divergence control does not set the result;
- gravity boundary, domain size, and Poisson tolerance are converged;
- full-speed/reduced-speed or justified RSLA sweeps bound radiation error;
- grey/multigroup and M1/hybrid-transport differences are quantified;
- conductivity uncertainties from chemistry and grains are propagated;
- explicit caps/floors do not control the reported flux, or their effect is included in the uncertainty;
- a statistically adequate turbulent/environmental ensemble exists;
- the field is followed to a physically defined protostellar surface, not inferred from a sink with unvalidated flux handling;
- numerical, microphysical, and sample-variance uncertainties are reported separately;
- all headline conclusions survive reasonable alternate analysis surfaces and definitions of core/stellar flux.

The desired end product is therefore not a claim of perfection. It is a reproducible causal and uncertainty budget showing where the initial magnetic flux went, which physical processes moved or dissipated it, and how confidently the surviving protostellar field can be predicted.

---

## Final whole-tree recheck and GPU performance addendum

This addendum records a second independent static pass over all 121 application source files (26,374 lines under `src/`), the build caches and compile databases, production inputs/scripts, test inventory, and retained run telemetry. It does not supersede the findings above; it adds three confirmed correctness defects and a performance assessment.

### Additional confirmed findings

#### A1. HIGH — physics / MPI-AMR — radiation moments are not refluxed across coarse-fine interfaces

- **Status:** CONFIRMED.
- **Location:** `src/radiation/radiation_moments.cpp:450-454`; compare the hydro flux-correction sequence at `src/hydro/hydro_driver.cpp:528-551` and the Poisson sequence at `src/self_gravity/poisson_equation.hpp:49-54`.
- **One-sentence defect:** Every radiation subcycle computes fluxes and immediately applies the update without Parthenon's flux-correction receive/send/set tasks, so coarse and fine representations use different net radiation fluxes at an AMR interface.
- **Concrete failure:** When a radiation pulse, first-core photosphere, or radiative shock crosses a refinement boundary, the sum of fine-face `Er`/`Fr` fluxes need not equal the coarse-face flux. The update therefore creates or destroys radiation energy and momentum at the interface; subsequent matter coupling transfers that error into gas energy and momentum. The error depends on refinement placement and can look like physical heating/cooling.
- **Evidence:** The radiation fields are registered with `Metadata::WithFluxes` at `src/radiation/radiation.cpp:165-184`, but `AddRadiationTasks` contains only `CalculateRadFluxes -> ApplyRadUpdate -> AddBoundaryExchangeTasks`. No radiation call to `StartReceiveFluxCorrections`, `LoadAndSendFluxCorrections`, `ReceiveFluxCorrections`, or `SetFluxCorrections` exists.
- **Fix:** For every radiation subcycle, insert a radiation-only reflux sequence after flux construction and before `ApplyRadUpdate`; ensure the correction pack contains exactly `rad.Er` and `rad.Fr1-3`.
- **Demonstration:** Advect/free-stream a radiation packet and run an optically thick diffusion pulse across a fixed coarse-fine interface. Compare domain-integrated radiation-plus-gas energy, interface flux balance, and L1 error for uniform and AMR meshes through at least three resolutions. Repeat with the interface displaced by half a coarse domain.

#### A2. HIGH — numerics / MPI — matter coupling leaves stale radiation ghosts

- **Status:** CONFIRMED.
- **Location:** `src/radiation/radiation_moments.cpp:451-459`.
- **One-sentence defect:** The final transport boundary exchange occurs before `MatterCoupling`, although coupling modifies all four radiation moments, and no exchange follows it.
- **Concrete failure:** On the next hydro step, `CalculateRadFluxes` reconstructs an interface state from newly coupled interior moments on one side and pre-coupling ghost moments on the other. MPI and MeshBlock boundaries receive a one-step, decomposition-dependent radiation flux impulse. Strong absorption/emission in the first-core photosphere maximizes the error.
- **Evidence:** Line 453 exchanges `md_rad`; line 459 subsequently changes the radiation moments through `MatterCoupling`. The line-457 statement that no exchange is needed is false for the radiation variables consumed by the next transport flux calculation.
- **Fix:** Exchange the radiation-only container after coupling, or make the first operation of every subsequent radiation transport step a completed radiation ghost exchange. Preserve a dependency linking coupling to that exchange.
- **Demonstration:** Run an optically thick matter-radiation equilibration front with the front centered first inside a block and then on a rank boundary. Compare 1, 2, and 4 MPI decompositions cell-by-cell and monitor the first-step interface flux after coupling.

#### A3. HIGH — numerics / MPI — initial RKL2 non-ideal flux can use stale post-chemistry scalar ghosts

- **Status:** CONFIRMED for configurations in which diffusivity reads a chemistry scalar, including production `ionization_chem` ambipolar diffusion.
- **Location:** `src/hydro/hydro_driver.cpp:636-663`, `src/hydro/hydro_driver.cpp:239-265`, `src/hydro/diffusion/diffusion.cpp:287-301,336-351`, and `src/hydro/diffusion/ambipolar.cpp:165-176`.
- **One-sentence defect:** Chemistry changes interior conserved scalar abundances immediately before STS, but the initial STS diffusivity/flux task is not dependent on completion of a scalar boundary exchange and reads one ghost layer.
- **Concrete failure:** At a block, rank, or AMR boundary, `PrecomputeNonidealEta` obtains `x_e` from the reacted interior on one side and an old ghost value on the other; ambipolar face diffusivity averages those two cells. The first RKL2 half-step therefore applies a spurious, decomposition-dependent non-ideal EMF exactly where ionization changes steeply.
- **Evidence:** Chemistry is added at lines 642-644 and only a local `FillDerived` follows at 652-657. `AddSTSTasks` begins at line 663. Within STS, receive tasks start at 244-246, but `CalcDiffFluxes` depends only on `ResetFluxes` at 257-258 and precedes completion of boundary exchange. The eta kernel explicitly spans `interior +/- 1`, and the face kernel reads/averages `i-1` and `i` chemistry values.
- **Fix:** Complete a conserved-scalar boundary exchange plus `FillDerived` after chemistry and before the initial STS flux, or make initial `CalcDiffFluxes` depend on an STS-entry exchange restricted to variables required by diffusivity.
- **Demonstration:** Construct a stationary uniform MHD state with a sharp but analytically prescribed `x_e` transition at a rank boundary and constant physical diffusivity after the prescribed mapping. Compare the first STS update with the transition inside a block and across 1/2/4-rank decompositions. Require convergence to roundoff for equivalent decompositions before testing physical diffusion convergence.

### Final recheck disposition

The second pass did not provide evidence that invalidates the earlier confirmed findings. It did identify historical comments in `runs/prod_t4_full/submit.sh:50-54` describing an older broad `{Independent}` STS-pack overwrite of radiation and gravity. The current source copies only `cons` and `prim` into `u1` at `src/hydro/hydro_driver.cpp:207-220`, so that particular historical defect is not reported as a current-source finding. This distinction is important: archived comments are provenance clues, not validation evidence.

No new confirmed host-pointer/device-lambda violation was found in the second pattern-and-call-site pass. Host mirrors and MPI calls found in EOS-table loading, gravity reductions, turbulence setup, and sink synchronization occur outside device lambdas. That does not establish race freedom for every optional cluster/sink path; those paths remain less deeply runtime-tested than the collapse path, as stated in the coverage table.

## GPU speed assessment

### What was and was not measured

Live Nsight Systems/Compute profiling could not be performed in this audit environment: the frontend has no accessible NVIDIA driver (`nvidia-smi` cannot communicate with a driver), and `nsys`, `ncu`, and `nvprof` are not installed in the active environment. A direct CPU smoke launch also failed during frontend MPI/PSM2 initialization before application execution. No retained `.nsys-rep`, `.qdrep`, `.ncu-rep`, or profiler SQLite artifact exists under `runs/`. Consequently:

- scheduler telemetry and AthenaPK throughput below are **measured archived observations**;
- source-derived bottlenecks are **mechanistic optimization candidates**;
- speedups quoted only in shell-script comments are **unverified historical assertions** until the underlying profiler report or a new controlled A/B run is produced.

The active GPU build cache is internally consistent with the stated target: `build_gpu/CMakeCache.txt` is `Release`, uses `-O3 -DNDEBUG`, has `Kokkos_ENABLE_CUDA=ON`, and `Kokkos_ARCH_HOPPER90=ON`. The CPU cache is also `Release -O3 -DNDEBUG` with CUDA disabled. No project-wide fast-math flag appears in either compile database.

### Observed production performance

| Evidence | Workload / allocation | Throughput | GPU use | Memory and imbalance |
|---|---|---:|---:|---|
| `runs/prod_t4_full/prod_t4_full_2360418.out:19,96-124` | full tier-4 collapse, 1 node, 5 H100, 5 MPI ranks | `4.12e6` zone-cycles/s | scheduler average 0.7/5 GPU equivalents (14%) | 55.6-56.3 GiB/GPU; per-device GPU load 28%, 33%, 33%, 39%, 62% |
| `runs/prod_t4_full/prod_t4_full_2367303.out:93-121` | same campaign family, 5 H100 | terminal throughput not printed in this file | scheduler average 0.7/5 (14%) | 56.7-58.0 GiB/GPU; load 22%, 22%, 23%, 39%, 77% |
| `runs/prod_t5_smallbox/t5_smoke_2372692.out` | small-box smoke, 5 H100 | about `3.0e6` zone-cycles/s on ordinary steps; about `1.89e6` on an output step | roughly 28-29% per GPU; scheduler aggregate 18% | about 12 GiB/GPU |

Across retained production summaries, reported terminal throughput spans approximately `0.943e6` to `4.12e6` zone-cycles/s; recent full-tier runs commonly report `3.27e6-4.12e6`. These values are not a controlled scaling curve because mesh state, AMR population, timestep, output, and physics workload differ.

The primary performance conclusion is nevertheless robust: the campaign is not saturating five H100s. The scheduler reports only about 14% aggregate GPU-equivalent use in the representative full runs, while one rank/GPU reaches 62-77% and others remain near 22-39%. This is severe work/critical-path imbalance plus launch/communication latency, not a simple lack of peak H100 arithmetic. GPU memory, at 56-58 GiB per card, constrains rank-count experiments but is not at the 80-GiB hardware ceiling. Only five of forty allocated CPU cores are materially used; that is allocation waste rather than the main GPU-speed limiter.

The production script at `runs/prod_t4_full/submit.sh:24-41` records a prior claim that `PrecomputeNonidealEta` consumed 65.7% of GPU kernel time and that frozen eta reduced step time from 43.9 s to 20.1 s. Because the actual profile is absent, this is classified as **SUSPECTED/archived**, although the source makes the mechanism plausible: the kernel performs EOS inversion and a conductivity calculation over interior plus ghost cells at every refresh.

## GPU optimization plan, ordered by expected scientific value

All performance changes must be benchmarked only after the correctness defects above are fixed. Each optimization must preserve a frozen reference result within a predeclared tolerance and must report both zone-cycles/s and simulated-time/s; the latter catches changes that make a step cheaper by shrinking the physical timestep.

### P0 — obtain an actionable profile and load-balance trace

1. Build a `RelWithDebInfo` Hopper binary with CUDA line information and Kokkos profiling hooks, without fast math.
2. Profile 20-50 ordinary steps from a copied representative checkpoint, excluding startup and output, using one Nsight Systems run. Capture CUDA kernels, CUDA API, MPI, NVTX/Kokkos regions, OS runtime, and per-rank traces. Repeat a separate short window containing AMR and a separate output step.
3. Run Nsight Compute only on the top 5-10 kernels identified by Systems, collecting launch dimensions, achieved occupancy, registers, L1/L2/DRAM throughput, branch efficiency, eligible warps, and roofline data. Do not profile the full application with all metric sets.
4. Record per-rank MeshBlock counts by level and per-package wall time. A five-GPU average hides the observed 3x rank-load spread.
5. Establish controlled 1/2/4/5/8-GPU strong-scaling and weak-scaling points with identical physics, output disabled, fixed mesh, identical cycle interval, and at least five repeats. Report median and dispersion.

Acceptance: the trace must account for at least 90% of wall time as kernel, CUDA/MPI wait, boundary/AMR, gravity, radiation, chemistry, diffusion, or I/O time. Archive the profiler files with the commit, binary hash, input, restart hash, module list, and GPU clocks.

### P1 — eliminate physics-weighted rank imbalance

The 22-77% per-device utilization spread makes load balance the first likely end-to-end win.

- Instrument cost per block separately for hydro, gravity iterations, radiation subcycles, chemistry, non-ideal eta, STS stages, and particle work.
- Replace a zone-count-only block weight with a measured cost model including AMR level and enabled/stiff physics. Refit the weights at remesh intervals, not every step.
- Verify that partitions contain comparable active block counts and that no rank owns the critical finest-level/chemistry-heavy blocks disproportionately.
- Test 4 and 8 ranks/GPUs. Five-way placement may be necessary for memory today, but reducing scratch/duplicate state can make four GPUs feasible and eight may reduce critical-path imbalance.

Acceptance: maximum per-rank step time divided by median below 1.10 for fixed-mesh windows, and parallel efficiency above 70% from 1 to 4 GPUs for a workload that fits one GPU.

### P2 — remove redundant STS memory and work

- `src/hydro/hydro_driver.cpp:232-235` allocates `MY0` and `Yjm2` from the full base container even though only the STS-updated state is required. Construct narrow registers containing only hydro conserved/primitive data actually used by RKL2.
- `ResetFluxes` at `src/hydro/hydro_driver.cpp:49-90` selects all `Metadata::Independent` flux variables and launches once per direction. Use a cached descriptor restricted to the hydro variables affected by enabled diffusion terms; zero each required flux while computing its first contribution where safe.
- Boundary exchange in every STS stage currently sends the broad base container (`src/hydro/hydro_driver.cpp:279-290,348-358`). Exchange only the variables whose stencil is consumed in the next diffusion stage.
- If the profile is launch-bound, capture the fixed RKL stage graph with CUDA Graphs/Kokkos Graphs and fuse coefficient/update/positivity work where data dependencies permit.

Acceptance: lower peak GPU memory by at least 20%, no change to the chosen RKL2 eigenmode convergence order, and bitwise or tightly bounded agreement for a fixed stage count.

### P3 — accelerate EOS inversion and non-ideal coefficients

The tabulated EOS and conductivity path is a likely arithmetic/branch-divergence hotspot.

- Profile the number and cost of `TemperatureKFromPres`/internal-energy inversions and conductivity tensor solves per cell, face, STS stage, and hydro stage.
- Provide monotone inverse EOS tables directly in the variables needed by hot kernels, e.g. `(rho,P) -> (T,e,c_s^2)`, rather than executing a fixed iterative inversion at every face/cell. Preserve thermodynamic derivatives and table continuity.
- Cache cell-centered `T`, sound speed, and all three diffusivities once per accepted physical state. Reuse them only over a quantitatively bounded state-change interval; freezing eta is a numerical approximation, not a free optimization.
- Consider an adaptively refined 3D/4D conductivity table or interpolation surrogate in `(rho,T,|B|,x_e/grain state)`. Generate it from the same microphysics and bound interpolation error in log eta, including sign for Hall.
- Skip recomputation only for cells whose relative changes in the controlling variables prove an eta-error bound; force refresh near shocks, chemistry transitions, refinement changes, and core thresholds.

Acceptance: maximum and percentile errors for `T`, `c_s`, `eta_O`, `eta_A`, and signed `eta_H` below predeclared physics tolerances over the entire table domain; identical convergence class in shock, diffusion, and collapse tests; a sweep of refresh thresholds showing the fossil-flux observable is insensitive within its numerical error bar.

### P4 — reduce launch and synchronization overhead

- Fuse the three timestep reductions (hyperbolic, parabolic, Hall) where their inputs overlap, and combine the MPI minima into one collective.
- Cache mesh `dx_min` used at `src/radiation/radiation_moments.cpp:408-423` and hydro minimum-cell-size reductions between remesh events. Geometry does not change every cycle.
- Overlap nonblocking global reductions and halo exchange with interior kernels; measure whether GPU-aware MPI is active rather than assuming it.
- Increase useful work per launch by testing 64^3 blocks against the current smaller blocks, subject to AMR localization and memory limits. Tune Kokkos team/vector policies from NCU occupancy and memory-transaction evidence.
- Add Kokkos profiling regions around every major package and NVTX ranges around task regions so scheduling gaps are attributable.

Acceptance: CUDA API/launch plus device-idle gaps below 15% of fixed-mesh step time, with no regression in AMR accuracy or memory capacity.

### P5 — gravity and radiation algorithms

- For self-gravity, profile smoother, residual, restriction/prolongation, reflux, boundary multipoles, and global norms separately. Warm-start potential and reuse hierarchy data; select smoother parameters from measured spectral behavior. Terminate against a discretization-error-scaled residual, not an unnecessarily tight fixed residual.
- Remove atomic contention from multipole accumulation (`src/self_gravity/self_gravity.cpp:301-336`) using hierarchical/team reductions if isolated multipoles are enabled.
- Radiation currently incurs several kernels and a boundary exchange per subcycle. First fix reflux/ghost correctness; then investigate an asymptotic-preserving or implicit transport formulation that relaxes the reduced-light CFL in optically thick cells. Kernel fusion alone cannot compensate for an excessive subcycle count.
- Cache radiation `dx_min`; overlap radiation interior fluxes with halo receipt; communicate only four radiation moments.

Acceptance: gravity force/error benchmarks remain unchanged at the accepted Poisson tolerance; radiation free-streaming and diffusion limits retain their convergence; reduced-speed and subcycle changes remain inside the radiation error budget defined earlier in this report.

### P6 — chemistry, particles, and I/O

- Group chemistry cells by stiffness or use an active-cell mask to reduce warp divergence, while preserving positivity and elemental conservation. Tabulate only rate combinations demonstrated to dominate.
- Sink routines perform host mirrors and global all-gathers/all-reduces (`src/sinks/sinks.cpp:128-146,486-583,676-685`). Replace all-gathered global sink state with spatial ownership/neighbor communication when sink counts grow, and keep integration/accretion state device-resident.
- A representative full slot wrote 423.7 GiB. Separate output-step timing from compute timing, reduce duplicated/unneeded fields, use collective/parallel HDF5 settings established by measurement, and consider asynchronous staging. Never reduce restart cadence below the recoverability requirement.
- Allocate only the CPU cores actually used, or deliberately use spare cores for MPI progress/asynchronous I/O after proving a gain.

Acceptance: output overhead below the campaign's declared fraction of wall time; restart files remain complete and round-trip tests pass; chemistry/sink conservation is unchanged.

### Required optimization experiment matrix

For every candidate, run A/B/A in one allocation to control clock and filesystem variation. Use the same checkpoint, fixed cycle range, GPU clocks, rank mapping, and output policy. Report:

- median step wall time and 95% interval;
- zone-cycles/s and simulated-time/s;
- per-package and per-rank time;
- kernel count, average duration, occupancy, DRAM bandwidth, and device-idle fraction;
- maximum GPU memory;
- conservation residuals, `div B`, Poisson residual, radiation-plus-gas energy balance, species conservation, and minimum positivity margins;
- the primary fossil-field diagnostic at the end of a short sensitive-window restart.

Promote an optimization only when it is faster in simulated-time/s, passes the full regression suite, and introduces no uncontrolled approximation. In particular, eta caps, eta freezing, reduced light speed, looser Poisson tolerances, and lower output precision are physics/numerics changes and must never be presented as purely computational optimizations.

### Performance conclusion

The current H100 run is fast in raw throughput—up to about 4.12 million zone-cycles/s in retained full-production telemetry—but inefficient in hardware use. Five H100s deliver roughly 0.7 GPU-equivalents of average scheduler load in representative runs, with pronounced rank imbalance. The most credible near-term gains are therefore physics-aware block rebalancing, narrower STS storage/communication, removal or controlled caching of repeated EOS/conductivity work, and reduction of small-launch/global-synchronization overhead. A defensible numerical speedup cannot be quoted until a new archived Nsight trace and controlled scaling/A-B campaign are completed.

---

## Executed H100 profile — 2026-07-21

After the static audit, an actual Nsight Systems profile was run on an H100 80 GB HBM3 compute node (driver 580.159.04). The profiler was Nsight Systems 2025.1.3; the application remained the CUDA-12.5 Hopper build `athenaPK_eos_v7`, MD5 `505491e86a9e90f8a470f2bc298939d6`. The valid trace and derived data are:

- `profiling/final_gpu_audit/nsys_direct1_2025.nsys-rep` (53 MB);
- `profiling/final_gpu_audit/nsys_direct1_2025.sqlite` (144 MB);
- `profiling/final_gpu_audit/direct1_stats_cuda_gpu_kern_sum.csv`;
- `profiling/final_gpu_audit/direct1_stats_cuda_api_sum.csv`;
- `profiling/final_gpu_audit/direct1_stats_osrt_sum.csv`;
- `profiling/final_gpu_audit/direct1_analyze.txt`.

The profiled state was the earliest v7 small-box checkpoint: 512 MeshBlocks, full tabulated EOS, gravity, M1 radiation, chemistry, Ohmic/ambipolar/Hall physics, and RKL2 enabled. It ran two completed post-restart timesteps. It is a valid execution profile of the current code and physics, but it is **not** representative of the late-collapse STS stiffness: this state used three RKL2 stages per half, whereas the later checkpoint uses thirteen. The trace also includes restart/setup work, so API call totals must not be interpreted as pure per-step rates. A five-H100 late-state trace is queued separately as Slurm job `2375318`.

### Measured kernel-time breakdown

Total CUDA kernel time in the trace was approximately 16.0 s. The leading kernels were:

| Kernel family | GPU kernel time | Calls | Interpretation |
|---|---:|---:|---|
| `SelfGravity::FillPoissonRHS` | 19.7% | 17 | Gravity RHS construction is the single largest named kernel in this state; 186 ms per launch. |
| fused non-ideal timestep reduction | 9.4% | 3 | About 502 ms per reduction; conductivity/EOS work in timestep control is load-bearing. |
| radiation update | 7.9% | 948 | About 1.34 ms per radiation subcycle. |
| radiation flux, three directions | 21.4% total | 2,844 | About 1.19-1.23 ms per direction/subcycle. |
| Parthenon `SetBounds` | 8.4% | 1,011 | Boundary processing is a major GPU consumer. |
| Parthenon send/restrict boundary packing | 4.0% | 1,011 | Adds to the boundary/communication cost. |
| `PrecomputeNonidealEta` | 4.0% | 8 | About 80.8 ms per refresh even in the early, relatively easy state. |
| chemistry reaction reduction | 1.7% | 2 | 77-190 ms per call with large dispersion. |
| RKL2 recurrence update | 1.5% | 8 | About 30.5 ms per stage update. |
| principal HLLD hydro flux kernels shown in top table | about 4% | 8 | Roughly 75-80 ms per direction/reconstruction variant. |

The dominant result is radiation subcycling: 948 update calls over two completed timesteps means 474 radiation subcycles per timestep in this early state. Its three flux kernels plus update consume 29.3% of all measured GPU kernel time before including the associated boundary work. This directly confirms that the reduced-light CFL/subcycling design, not merely implementation details within one radiation kernel, is a primary performance limiter.

### Measured launch, allocation, and synchronization pressure

The trace recorded:

- 468,900 `cudaLaunchKernel` calls;
- 170,801 `cudaDeviceSynchronize` calls;
- 227,854 `cudaStreamSynchronize` calls;
- 24,898 `cudaEventSynchronize` calls;
- 101,858 `cudaMalloc` and 101,858 `cudaFree` calls;
- 101,852 asynchronous copies and 53,150 synchronous copies.

These totals include checkpoint reconstruction and setup, but their scale is unequivocally excessive for a trace containing only two accepted timesteps. CUDA API time was led by device synchronization (30.8%), stream synchronization (19.4%), allocation (12.9%), free (12.5%), and event synchronization (11.1%). Kernel-launch API overhead itself was only 8.2%; the larger cost is synchronization and dynamic resource management around the launches.

Nsight's automated analysis independently flagged pageable-memory asynchronous copies, synchronous host/device copies, and blocking CUDA synchronization. Examples include repeated synchronous pageable-to-device transfers of 91.386 MB, 52.221 MB, 43.706 MB, and many 5.225 MB blocks. These are concrete transfer paths to replace with persistent device storage, pinned staging buffers, and asynchronous dependency-driven copies where application ownership permits.

No NVTX ranges were present even though NVTX tracing was enabled. Package/task-level attribution therefore had to be reconstructed from demangled kernel names. Adding Kokkos Tools/NVTX regions around driver stages, radiation subcycles, gravity solves, STS halves, chemistry, AMR, and output is a prerequisite for a maintainable performance audit.

### Additional confirmed GPU-performance findings

#### P1. HIGH — GPU/numerics — explicit radiation subcycling dominates launch count and at least 29.3% of kernel time

- **Status:** CONFIRMED by profile.
- **Location:** `src/radiation/radiation_moments.cpp:408-459`.
- **Mechanism:** `nsub = ceil(dt/dt_rad)` and each subcycle launches three directional flux kernels, an update kernel, and boundary work. The observed state executed 474 subcycles per hydro step.
- **Failure scenario:** As AMR reduces `dx_min` or a less aggressive reduced light speed is required for accuracy, the global finest-cell radiation CFL multiplies work on every block, including blocks that do not need that local subcycle rate. The cost rises independently of local optical stiffness.
- **Fix/experiment:** First repair radiation reflux and ghost ordering. Then compare global subcycling with level/local subcycling and an asymptotic-preserving or implicit transport method. Sweep reduced light speed against physical radiation convergence, reporting subcycle count, energy error, and simulated-time/s.

#### P2. HIGH — GPU/code — extreme synchronization and allocator churn prevents asynchronous execution

- **Status:** CONFIRMED by profile; separation of setup versus steady-state contributions requires a delayed late-state trace.
- **Location:** application-wide task/container use, with prominent broad register creation at `src/hydro/hydro_driver.cpp:224-235` and repeated stage/boundary scheduling at `src/hydro/hydro_driver.cpp:239-358`.
- **Mechanism:** The trace contains over 100,000 device allocations/frees and nearly 400,000 device/stream synchronizations. These operations serialize the host/device pipeline and make many otherwise asynchronous kernels execute in short synchronous bursts.
- **Fix/experiment:** Preallocate and reuse MeshData/register/communication buffers; narrow STS containers; cache pack descriptors; use a device memory pool; replace host synchronization with stream/event dependencies. Reprofile setup and ten fixed-mesh steps separately and require allocation count after setup to approach zero.

#### P3. HIGH — GPU/physics — non-ideal timestep estimation is nearly as expensive as the eta update it controls

- **Status:** CONFIRMED by profile.
- **Location:** `src/hydro/diffusion/diffusion.cpp:215-284` and its fused ionization evaluator.
- **Mechanism:** Three fused timestep reductions consumed 9.4% of kernel time at about 502 ms each; eight `PrecomputeNonidealEta` kernels consumed another 4.0% at about 80.8 ms each. The timestep path therefore recomputes expensive microphysics independently of the eta cache used by fluxes.
- **Fix/experiment:** Compute eta extrema and timestep candidates while filling the eta cache, reduce the cached values once, and refresh only when the state validity criterion requires it. Verify identical `dt_par`, `dt_strict`, chosen RKL stage count, and diffusion convergence before accepting the fusion.

#### P4. MEDIUM — GPU/build — no package-level GPU profiling ranges are emitted

- **Status:** CONFIRMED.
- **Location:** driver and package task construction generally; the valid trace's `nvtx_sum` report is empty.
- **Mechanism:** Without NVTX/Kokkos profiling ranges, idle intervals, task dependencies, MPI wait, and kernel groups cannot be robustly assigned to physics packages, which impedes regression detection and rank comparison.
- **Fix:** Add low-overhead named regions around every package/stage and preserve them in production-capable builds. Add a CI check that a one-step trace contains expected ranges.

### Revised optimization order from measured evidence

1. Repair radiation correctness, then reduce the 474 global radiation subcycles per hydro step; this is the largest directly measured algorithmic opportunity.
2. Fuse/cache non-ideal timestep and eta evaluation; together they account for 13.4% of kernel time in the early state and will grow with late-collapse STS work.
3. Eliminate post-setup device allocation/free and broad STS registers/communications.
4. Remove unnecessary device/stream synchronizations and pageable/synchronous transfers.
5. Reduce boundary kernel multiplicity and exchange only package-specific variables; radiation plus boundary kernels dominate much of the remaining trace.
6. Optimize gravity RHS construction and the multigrid solve, using late-state/package ranges to distinguish physical work from restart initialization.
7. Only after those changes, tune individual HLLD, radiation, chemistry, and RKL kernels with Nsight Compute. Kernel micro-tuning before launch/synchronization/algorithmic fixes would attack the smaller term.

This profile strengthens, rather than weakens, the earlier conclusion: H100 utilization is limited primarily by algorithmic subcycling, task/kernel multiplicity, synchronization, allocation, and package boundary work. It is not limited by lack of nominal Hopper compute throughput.
