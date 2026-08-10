# AthenaPK — magnetized BE-sphere / first-hydrostatic-core collapse (flagship tree)

A fork of [AthenaPK](https://github.com/parthenon-hpc-lab/athenapk) (Parthenon + Kokkos),
reduced on **2026-08-10** to the flagship configuration: the magnetized Bonnor-Ebert sphere
gravitational-collapse problem, the physics it needs, and its production runs.

Branch `flagship-phase2-ct`. The upstream README is preserved at
[`docs/UPSTREAM_README.md`](docs/UPSTREAM_README.md).

## What this tree is

One problem generator — `collapse_be` — plus a substantial custom physics stack on top of
upstream AthenaPK:

| package | `src/` | notes |
|---|---|---|
| self-gravity | `self_gravity/` | multigrid Poisson, ported from Artemis. Multipole exterior BCs (`*_bc=multipole`) validated to 0.54%; zero-Dirichlet is the fallback; **Neumann is broken** |
| non-ideal MHD | `hydro/diffusion/` | Ohmic, ambipolar (η_A = Q_A·B²), Hall. AD works under `unsplit` and RKL2 STS; **Hall is dispersive ⇒ `unsplit` only**, and needs an Ohmic floor |
| constrained transport | `hydro/ct/` | face-centered B. GLM/Dedner is the alternative `divergence_control` |
| M1 radiation transport | `radiation/` | multigroup (`n_group`), tabulated + Bell & Lin opacity, Planck/Rosseland split, RSLA (`creduc`) |
| tabulated EOS | `eos/` | multi-Saha + H₂ dissociation, second-core capable |
| ionization / conductivity | `hydro/diffusion/ionization.hpp` | NICIL-class MRN-grain Wardle tensor, x_e(ρ,T) → η_A/η_H/η_O |
| chemistry | `chemistry/` | GPU H₂ and reduced-GOW17 networks via passive scalars |
| sinks | `sinks/` | Parthenon swarm, restart- and AMR-safe |
| dust | `dust/` | growth, sublimation, opacity coupling |

Everything past ideal MHD is **gated off by default**, so an ideal or AD-only production
run is unaffected by it.

## Layout

```
src/          flagship source. src/pgen/ contains collapse_be.cpp and nothing else.
inputs/       collapse_be.in, collapse_be_fhc_production.in
runs/         production only -- see runs/README.md
docs/         flagship documentation; HANDOFF_2026-08-10.md is the cold-start entry point
build_cpu/    Kokkos OpenMP (front-end builds are fine here)
build_gpu/    Kokkos CUDA, Hopper90 -- must be built through SLURM
DEV_LOG.md    the running development log
```

**The validation material is not here.** Problem generators, input decks, the regression
suite, the WP01–WP22 reports and 99 validation run directories live in
`/beegfs/u/bbg6470/validation/` — start at its `README.md`. The complete pre-split tree is
in git at tag `validation-complete-2026-08-10`.

## Building

```bash
source ~/athenapk_env.sh          # GCC 13.3 / OpenMPI 5.0.7 / HDF5 / CMake

# CPU
make -C build_cpu athenaPK -j

# GPU: the front end has no GPU and a 1-CPU compile cgroup, so go through SLURM
sbatch --export=ALL,TMPDIR=/tmp runs/submit_build_gpu.sh
# -> build_gpu/bin/athenaPK
```

## Running

```bash
sbatch runs/root_ladder/submit_root.sh        # never launch MPI on the front end
```

Submit scripts are timeout-safe and self-resuming: on start each picks the newest restart
(`parthenon.out2.NNNNN.rhdf`) and continues, or starts fresh from t=0 if there is none.

## Conventions that bite

- **Magnetic units are Heaviside-Lorentz**: `v_A = B/√ρ`, `P_mag = B²/2`, `β = 2P/B²`.
  Athena++ is Gaussian. The *code-unit number* for B is the same in both codes — each
  absorbs its own 4π — so `B0z = 0.15` (HL) ↔ Athena++ `mu = 31.9` ↔ 7.47 µG. In analysis,
  the √(4π) in `v_A` is the only line that differs between the two codes' notebooks.
  `units.json`'s `code_bfield_cgs` is the **Heaviside-Lorentz** unit; converting to Gauss
  needs another √(4π), a 3.545× error if taken at face value.
- **Code units are shared with Athena++** for the FHC IC (mass 6, T 10, f 5):
  `rho0 = 5.467e-19 g/cm³`, `v0 = 1.9e4 cm/s`, `l0 = 2.81e16 cm`, `t0 = 1.48e12 s`,
  `four_pi_G = 1`.
- **Output layout**: science is `parthenon.out1.NNNNN.phdf` (`prim`, `grav.phi`); restarts
  are `parthenon.out2.NNNNN.rhdf`. `prim` is `[block, component, k, j, i]` with
  0 = density, 1–3 = velocity, 4 = pressure, 5–7 = B, 8 = ψ, 9+ = scalars. Align raw reads
  against the `prim_*` entries of `Info/ComponentNames` — that attribute concatenates *all*
  datasets, `grav.phi` first.
- **Parameter defaults**: registering the same `(block, key)` twice with *different*
  defaults hard-aborts at startup, even when the deck sets the key explicitly. Paired
  defaults (e.g. the EOS table named in both `hydro.cpp` and `radiation.cpp`) must be
  changed together.
- An unknown key passed on the command line is **silently inert**. Check the binary
  (`strings`) actually knows a diagnostic before spending GPU time on it.
