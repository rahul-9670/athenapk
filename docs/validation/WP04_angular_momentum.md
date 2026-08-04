# WP-4 — Angular-momentum diagnostic

> **PRODUCTION-CONDITION DRIFT MEASURED 2026-08-04 (also closes B7's production leg).** From the
> deep AMR run (job 2453201, 204 rows, t = 0 → 1.01503, blocks 64 → 792, diode BCs, B ≠ 0,
> self-gravity + AMR on, `hydro/angmom_diag=true`):
>
> | axis | L(0) → L(end) | drift |
> |---|---|---:|
> | L_x | −1.6105e+02 → −1.5754e+02 | **+2.18 %** |
> | L_y | +2.2290e+01 → +2.0107e+01 | **−9.79 %** |
> | L_z | −1.1023e+02 → −1.0754e+02 | **+2.44 %** |
>
> **B7 context:** the previously quoted 11.7 / 10.0 / 7.9 % at 32/64/128 were from a *smoke* deck.
> Under production conditions the drift is **2–10 %** over a full t_ff — same order, and better on
> two of three axes. Angular momentum is still not conserved by Cartesian finite volume; this is
> now quantified where it matters.
>
> **The budget does not close, but the sampling cannot settle it.** Integrating
> −am-FT + am-Tgrav + am-Tmag over the history rows accounts for −1.13e-01 of L_x's actual
> +3.51e+00 change, leaving >100 % unexplained (and the wrong sign). That is *consistent* with
> B7 — but the history spacing is **1.78e-03**, far coarser than the timescale on which the
> surface flux varies, so a trapezoid over these rows is not an adequate budget integral. The
> honest claim is the **drift magnitudes above**; a closure statement would need per-cycle output.
> (Note B6 separately established that on this deck |dL/dt| exceeds every accounted term by ~10⁴,
> which is the same fact seen from the other side.)

**Status: instrument IMPLEMENTED, gated, and validated three independent ways.
Acceptance criterion in `VALIDATION_PLAN.md` CORRECTED (it was unachievable as written).
Two findings that change how L must be reported. 2026-07-31.**

Binaries: control `95ad54a7816ee3975435eb726215886d`; WP-4 v1 `d89fdcecc0247f8f08325290d4da318d`;
WP-4-corrected + WP-6 `7664869694cace9af38a35eb1c43b409` (all `build_cpu`, GCC 13.3, Kokkos
OpenMP, `OMP_NUM_THREADS=1`).

## The problem

The production history has 31 columns and **no angular momentum**. For a magnetic-braking
result L(t) and the torque budget *are* the measurement, so every braking statement in the
current output is an inference from `rho_max` and `ME/E`, not a measurement of L.

## What was implemented

`src/diagnostics/angmom_diag.{hpp,cpp}`, registered from `hydro.cpp`, gate key
**`<hydro> angmom_diag`, default `false`**.

| column | meaning |
|---|---|
| `am-Lx/Ly/Lz` | ∫ρ(r×v)dV |
| `am-FTx/FTy/FTz` | ∮ r×(T·n̂) dA — **the budget flux term**, full stress |
| `am-Tgravx/y/z` | ∫ρ r×(−∇Φ)dV — gravitational torque |
| `am-FLx/FLy/FLz` | ∮ρ(r×v)(v·n̂)dA — advective part alone (interpretive) |
| `am-Tmagx/y/z` | ∫ r×(J×B)dV — volume Lorentz torque (interpretive) |

with `T_ij = ρv_iv_j + (P + B²/2)δ_ij − B_iB_j` (Heaviside-Lorentz, no 4π). J = ∇×B by the
same 2Δx central difference `mag_diag` and the `jeans_nonideal` criterion use, so refinement
trigger and both budgets see one current. Moments are about the **box centre** from
`mesh_size`, echoed in the rank-0 banner. All flux columns are **outward-positive**.

Budget: **d(am-L)/dt = −am-FT + am-Tgrav.**

## Finding 1 — the budget specified in the plan misses by >2 orders of magnitude

`VALIDATION_PLAN.md` WP-4 asks for "the magnetic torque ∫(r×(J×B)) and the advective flux of
L through the outflow boundaries, so the budget can close". Implemented literally, on
`eos_smoke` (outflow BCs, B₀z=0.15):

| quantity | value |
|---|---|
| measured dLz/dt | **+12.5** code |
| max\|am-Tmagz\| | 8.1e−02 |
| max\|am-FLz\| | 7.5e−04 |

The two proposed terms are **~150× too small** to account for dL/dt. Two things were wrong:

1. **The advective flux is not the flux.** The angular-momentum flux through a surface is
   r×(T·n̂) with the *full* stress. Under outflow BCs the **surface pressure torque**
   −∮P(r×n̂)dA is the dominant term and the advective-only column omits it entirely.
   Hence `am-FT`.
2. **Adding `am-Tmag` to a surface flux double-counts.** r×(J×B) is the volume form of
   exactly the Maxwell part of that surface term. The budget is complete in the surface
   form; the volume form is the different, equivalent grouping
   `dL/dt = ∫r×(J×B)dV + ∫r×(−∇P)dV + Tgrav − am-FL`. Mixing terms across the two groupings
   is the easiest way to produce a budget that appears not to close.

`am-Tmag` is kept because the magnetic torque is the quantity a braking paper reports — but
it is a **diagnostic, not a budget term**, and the header says so.

## Gate — OFF-state bit-identical: **PASS**

Method per `VALIDATION_PLAN.md` §2 landmine 9: explicitly built control binary, matched
thread count, never a number from an old log. Case `runs/eos_smoke/fhc.in`, 12 cycles, 32³
root, self-gravity + EOS table + chemistry + RT + non-ideal MHD.

| comparison | `.hst` | snapshots |
|---|---|---|
| ctrl vs **off** (`angmom_diag` absent ⇒ default false) | **byte-identical** | **15/15 bit-identical** at `00000` *and* `final` |
| ctrl vs **on**, physics columns 1–13 | **byte-identical** | **15/15 bit-identical** |

Bit-identical both off *and* on, as expected for read-only reductions.

> **A false gate failure, recorded because it will recur.** The first attempt reported the
> OFF gate failing with `max|rel| = 2.6e4` on `/prim`. It was not a physics difference: a
> 10-minute tool timeout killed the launching shell but **orphaned the `athenaPK` child**,
> which was still writing into the run directory when the leg was relaunched — two processes
> appending one `.hst`. The tell was that the rows were out of order (cycle 0, 3, 1) while
> every value present matched the control exactly. **Before rerunning into a directory,
> check for a live writer** (`ps -eo pid,cmd -u $USER | grep athenaPK`).

## Instrument validation — three independent routes, all PASS

1. **Against numpy on the snapshot** (`docs/validation/scripts/wp4_angmom_check.py`, no
   shared code path with the Kokkos kernel). On `wp4_gate/on`:

   | | t = 0 | t = 1.121 (final) |
   |---|---|---|
   | `am-Lx` | hst 1.83286e+02 / numpy 183.2864 | hst 1.73153e+02 / numpy 173.1530 |
   | `am-Ly` | hst 3.93706e+01 / numpy 39.37060 | hst 3.62247e+01 / numpy 36.22471 |
   | `am-Lz` | hst −2.60503e+02 / numpy −260.5032 | hst −2.46510e+02 / numpy −246.5097 |

   Every value agrees to the **full 6-significant-figure precision the `.hst` prints**. (The
   apparent "rel ≈ 1e−6" is text quantisation, not disagreement — the WP-5 trap.)

2. **Axisymmetry.** On the rotating, turbulence-free acceptance IC, `am-Lx` and `am-Ly` must
   vanish identically. Measured `|Lx|,|Ly| ≤ 7e−14` against `|Lz| = 281`, i.e. **2.6e−16
   relative — machine precision**. In the numpy recomputation they are exactly `0.0`.

3. **Zero-flux.** With reflecting walls, `max|am-FTz| = 1.3e−12`. The flux kernel's boundary
   detection fires on the right cells and returns zero when nothing crosses.

## Finding 2 — angular momentum is NOT conserved by the scheme, and does not converge

### The acceptance criterion in the plan cannot be met, for a structural reason

The plan asks for "L conserved to machine precision" on a rotating, B=0, reflecting run. A
Cartesian finite-volume scheme cannot deliver that. Linear momentum telescopes exactly (the
flux leaving cell *i* is the flux entering *i+1*). Angular momentum does not: summation by
parts leaves `Σ_faces [A_x F^x_{ρv_y} Δx − A_y F^y_{ρv_x} Δy]`, which cancels only if the
numerical stress is symmetric, `T_xy = T_yx`. Those come from **different Riemann problems at
different locations** (an x-face and a y-face) and agree only to O(Δx²). So L is conserved at
best to truncation error. Expecting roundoff would make a correct instrument look broken.

Replaced with: **(a)** exact instrument validation (above), **(b)** measured drift required
to fall under refinement.

### The measurement — and (b) fails

Rotating BE sphere, reflecting walls, no turbulence, **uniform grid, AMR off**
(`refinement=none`, `numlevel=1`, `nbtotal` constant), `t = 0 → 1`.

| leg | N³ | Lz(0) | Lz(1) | drift | max\|am-FTz\| |
|---|---|---|---|---|---|
| `r32_nograv` | 32 | −2.890950e+02 | −2.553000e+02 | **+11.69 %** | 5.4e−12 |
| `r64_nograv` | 64 | −2.810380e+02 | −2.528210e+02 | **+10.04 %** | 1.3e−12 |
| `r128_nograv` | 128 | −2.694380e+02 | −2.481280e+02 | **+7.91 %** | 3.4e−11 |

Observed order **p = 0.22 (32→64) and 0.34 (64→128)** against 2.0 for a second-order scheme.
**The drift is very nearly resolution-independent.**

Everything that could source it was measured and excluded, not argued away:

| candidate | measured | verdict |
|---|---|---|
| boundary flux | max\|am-FTz\| ≤ 3.4e−11 | excluded |
| gravitational torque | run repeated with `self_gravity=false` | excluded |
| mass injection / density floor | mass drift **exactly 0** (2.7e−14 from snapshots) | excluded |
| AMR prolongation/restriction | `nbtotal` constant at 64, AMR off | excluded |
| the diagnostic itself | numpy on snapshots reproduces the same 10.04 % | excluded |
| linear-momentum bug | 1-mom, 3-mom conserved to ~1e−13 | excluded |

With self-gravity **on** at 64³ the drift is +4.37 %, with max|am-Tgravz| = 3.9e−2 and
max|am-FTz| = 3.4e−5 — sources ~300× too small to account for dLz/dt ≈ +12.

**The `am-Tgrav` stencil was verified against the code's own gravity kick.**
`self_gravity.cpp:519-539` forms `dpl+dpr = −(φ_{i+1} − φ_{i−1})` and applies
`ρ·(0.5·βΔt/Δx)·(dpl+dpr)`, i.e. exactly the 2Δx central difference `am-Tgrav` uses. The
diagnostic therefore measures precisely the torque the code applies, so "gravity is not the
source" is a measurement, not an assumption.

**Confidence.** *Verified*: the drift, its independence of gravity/flux/floor/AMR, and the
numpy reproduction. *Inferred*: that the residual is the Cartesian-FV stress-asymmetry term
argued above — consistent with linear momentum being exact while angular momentum is not,
but not yet isolated term-by-term. *Caveat, stated plainly*: the no-gravity legs blow the
pressure-supported sphere apart and drive shocked flow into the walls, where the scheme is
first-order at best; poor convergence there is not proof of poor convergence in the
production collapse. **The number that matters for the paper is the drift under production
conditions, which needs a GPU leg and is not yet measured.**

## Finding 3 — global L is carried by a tiny fraction of the mass (same pathology as WP-8)

Decomposition of the acceptance IC at 64³:

| cut | mass fraction | **Lz fraction** | volume fraction |
|---|---|---|---|
| ρ > 0.5 code | 0.0162 | **0.6385** | 5.2e−03 |
| ρ > 1.0 code | 0.0103 | **0.2409** | 2.1e−03 |
| R < 6.45 (the BE sphere) | 0.0617 | **1.0000** | — |

Rotation is applied **only inside the sphere**, which holds 6.2 % of the mass, and within it
**1.6 % of the total mass carries 64 % of Lz**. So a single global `am-Lz` is dominated by a
thin, highly-weighted shell — structurally the same ill-conditioning WP-8 found for
`mag-dissO`/`mag-dissA`, arriving here through the R² moment weighting rather than through
η(ρ).

Corroborating: across the ladder ∫ρR²dV converges cleanly (2.2664e7 → 2.2680e7 → 2.2685e7,
0.07 % then 0.02 %) and mass converges (5.0913e4 throughout), while **Lz(0) itself moves 7 %**
(−289.10 → −281.04 → −269.44) and is not converging. The moment of inertia and the mass are
resolved; the rotating velocity field at the sphere edge is not.

**Consequence for the paper: do not quote a global `am-Lz`.** The needed follow-up is a
density-split L, `am-Lz-hi/lo` about a `angmom_diag_rho_split`, exactly as `mag_diag` gained
`mag-dissO-hi/lo` — so the core budget and the envelope budget converge, or visibly fail to,
independently. Specified below as the next increment.

## Finding 3 follow-up — density-split L columns **IMPLEMENTED + GATED**

Binary `09e68f75f776d5c12dfb9374bb5a5059` (`build_cpu`, GCC 13.3, `OMP_NUM_THREADS=1`).

Added behind `hydro/angmom_diag_rho_split` (default **0 = off**): `am-Lx-hi/lo`,
`am-Ly-hi/lo`, `am-Lz-hi/lo`, and `am-Mhi` (= ∫ρdV over the hi side — without it the split
columns cannot be normalised into a "x % of the mass carries y % of L" statement).

**Gate — OFF-state bit-identical: PASS.** `runs/wp4_split_gate/off` vs `runs/wp4_gate/ctrl`
(binary `95ad54a7…`): `.hst` **byte-identical**, and **15/15 datasets bit-identical** at both
`00000` and `final`. All seven columns sit behind `rho_split > 0`, deliberately — appending
history columns shifts every downstream index and would silently break parsers reading `.hst`
files already on disk.

**Reconstruction check.** `runs/wp4_split_gate/on` (32³, `glmmhd`, `rho_split = 5.0`):
`|hi + lo − L|` is at most **0.28× the `.hst` print-resolution floor** on all three axes
(0.199 / 0.282 / 0.173 for x/y/z). That is *consistent with* exact summation, but the `.hst`
writes 6 significant figures and cannot resolve any finer — **exactness is not claimed**, and
the code comment says so. The two sides are separate Kokkos reductions, so bit-exactness is
not guaranteed under a non-deterministic reduction order.

**What it shows on this deck** (`rho_split = 5.0`, i.e. a much tighter cut than Finding 3's
ρ > 0.5): the concentration grows monotonically as the collapse proceeds —

| t | `am-Mhi`/mass | `am-Lz-hi`/`am-Lz` |
|---|---|---|
| 0.000 – 0.255 | 0.000 % | 0.000 % |
| 0.380 | 0.045 % | 0.074 % |
| 0.597 | 0.193 % | 1.122 % |
| 0.855 | 0.352 % | 3.025 % |
| 1.121 | **0.602 %** | **8.380 %** |

By t = 1.12 **0.6 % of the mass carries 8.4 % of Lz**, and the ratio is still climbing — the
instrument now measures the ill-conditioning directly instead of it having to be inferred from
`phdf` post-processing. Note this is a *different* cut from Finding 3's table, so the numbers
are not comparable to the 1.6 %/64 % figures there; the production cut should be chosen against
`rhocrit`, not carried over from this smoke deck.

## Remaining work

1. ~~Density-split columns~~ — **DONE**, see above.
2. **Production-condition drift** on a GPU leg with outflow BCs, B ≠ 0 and AMR on. The
   CPU-scale non-convergence above must not be extrapolated to it.
3. **Pressure torque** `∫r×(−∇P)dV` if the volume-form budget is ever wanted; the surface
   form is complete without it.
4. If the production drift is confirmed at the several-percent level, **it is the same order
   as the magnetic-braking signal**, and every braking number needs it quoted as a systematic.

## Reproducing

```bash
source ~/athenapk_env.sh
export OMPI_MCA_pml=ob1        # else MPI_Init dies in psm2_ep_open (Bus error)
export OMPI_MCA_io=romio341    # else every MPI_File_open fails, incl. reading the deck
export FI_PROVIDER=tcp
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false PMIX_MCA_gds=hash
env -C <rundir> build_cpu/bin/athenaPK -i <rundir>/fhc.in    # `cd` does not reach the process
```

Run dirs: `runs/wp4_gate/{ctrl,off,on}`, `runs/wp46_gate/{off,on}`,
`runs/wp4_accept/{r32,r64,r32_nograv,r64_nograv,r128_nograv}`,
`runs/wp4_split_gate/{off,on}` (density-split gate).
Analysis: `docs/validation/scripts/wp4_angmom_check.py` (`--drift <hst>` for the drift table).

Note the acceptance deck must use `fluid = euler`: AthenaPK rejects reflecting boundaries
under `glmmhd` ("Reflecting boundary conditions for MHD need special treatment"), and
`eos = hydrogen` requires `<physics> radiation = true`. Both are hard failures, not silent.
