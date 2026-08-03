# WP-14 — order of accuracy on a smooth problem, in the production configuration

**Status: PASS on all four families — fast 2.13, Alfven 2.08, slow 2.06, entropy 2.06.**
2026-08-02/03. Reproduce:

```bash
sbatch --export=ALL,WAVES="0 1 2 3" runs/wp14_order/submit_wp14.sh
/beegfs/u/bbg6470/venvs/analysis_env/bin/python docs/validation/scripts/wp14_order.py
```

## Why this WP existed, and why the stock test does not answer it

Before this, **AthenaPK had never been shown to achieve its formal order of accuracy in the
configuration production actually runs.** That is not a bookkeeping gap: every convergence claim
in the campaign (WP-7's root ladder, WP-8's njeans ladder) assumes the scheme converges at its
design rate, and none of them tests it — they measure whether *physics results* stop moving,
which a first-order scheme would also eventually do, just more slowly and to a different answer.

`inputs/linear_wave3d.in` exists but does **not** answer WP-14: it runs `fluid = euler` with
`riemann = hlle`. Production runs `fluid = glmmhd` with `riemann = hlld`. A convergence result
for the hydro HLLE path says nothing about the MHD HLLD path that carries the science.

`runs/wp14_order/lw_mhd.in` therefore copies every numerical knob from
`runs/root_ladder/fhc_rootladder.in`:

| knob | production | this test |
|---|---|---|
| `fluid` | `glmmhd` | `glmmhd` |
| `riemann` | `hlld` | `hlld` |
| `reconstruction` | `plm` | `plm` |
| `integrator` | `vl2` | `vl2` |
| `cfl` | 0.3 | 0.3 |
| `nghost` | 2 | 2 |

PLM + vl2 is formally **2nd order** in space and time, so the acceptance target is p → 2.

## Method

A smooth MHD eigenmode of amplitude 1e-6 is advected exactly one period, so the exact solution at
`tlim` **is** the initial condition — which makes `linear_wave_mhd::UserWorkAfterLoop`'s L1 norm a
true error rather than a difference against a finer run (no Richardson step, no extrapolation
assumption). Resolutions 16 → 32 → 64 → 128 (nx2 = nx3 = N, nx1 = 2N), observed order
`p = log2( L1(h) / L1(h/2) )`.

## Results — all completed families

Observed order `p = log2( L1(h) / L1(h/2) )`. Every family reaches p ~ 2 and the approach is
monotone from below, which is the signature of a genuinely 2nd-order scheme rather than a
coincidence at one resolution pair.

**fast-magnetosonic (`wave_flag = 0`)**

| N | RMS L1 | p(RMS) | p(rho) | p(M1) | p(E) | p(B2c) |
|---|---|---|---|---|---|---|
| 16 | 1.897e-07 | — | — | — | — | — |
| 32 | 5.588e-08 | 1.76 | 2.10 | 1.88 | 1.75 | 1.93 |
| 64 | 1.343e-08 | 2.06 | 2.12 | 1.96 | 2.07 | 2.07 |
| 128 | 3.070e-09 | **2.13** | 2.08 | 2.11 | 2.13 | 2.18 |

**Alfven (`wave_flag = 1`)**

| N | RMS L1 | p(RMS) | p(rho) | p(M1) | p(E) | p(B2c) |
|---|---|---|---|---|---|---|
| 16 | 2.067e-07 | — | — | — | — | — |
| 32 | 6.005e-08 | 1.78 | 2.52 | 1.78 | 2.42 | 1.85 |
| 64 | 1.526e-08 | 1.98 | 2.58 | 1.94 | 2.62 | 1.99 |
| 128 | 3.613e-09 | **2.08** | 2.50 | 2.05 | 2.66 | 2.08 |

**slow-magnetosonic (`wave_flag = 2`)**

| N | RMS L1 | p(RMS) | p(rho) | p(M1) | p(E) | p(B2c) |
|---|---|---|---|---|---|---|
| 16 | 1.954e-07 | — | — | — | — | — |
| 32 | 6.669e-08 | 1.55 | 1.55 | 1.54 | 1.46 | 1.83 |
| 64 | 1.751e-08 | 1.93 | 1.93 | 1.92 | 1.92 | 1.97 |
| 128 | 4.212e-09 | **2.06** | 2.06 | 2.04 | 2.06 | 2.05 |

**entropy (`wave_flag = 3`, `vflow = 1.0`)**

| N | RMS L1 | p(RMS) | p(rho) | p(M1) | p(E) | p(B2c) |
|---|---|---|---|---|---|---|
| 16 | 1.097e-07 | — | — | — | — | — |
| 32 | 4.091e-08 | 1.42 | 1.42 | 1.42 | 1.42 | — |
| 64 | 1.053e-08 | 1.96 | 1.96 | 1.96 | 1.96 | — |
| 128 | 2.531e-09 | **2.06** | 2.06 | 2.06 | 2.06 | — |

`p(B2c)` is omitted for the entropy family rather than reported: the entropy mode carries no
transverse-field perturbation, so `L1(B2c)` sits at **4.4–4.8e-16** — machine epsilon — on every
rung. The "order" computed from roundoff (−0.07) is meaningless. That the transverse field stays
at 1e-16 through the run is itself a correct-behaviour check.

**PASS on all four.** Asymptotic p(RMS): fast **2.13**, Alfven **2.08**, slow **2.06**, entropy
**2.06**. Every family approaches 2 monotonically from below. The slow and entropy families
converge more slowly at coarse resolution (1.55 and 1.42 at 16→32) and recover by 64→128 —
expected, since both are the least well resolved per wavelength at fixed grid.

### The entropy family had to be run twice

The first attempt was a **null test**, and the failure mode is worth recording because nothing in
the output flagged it. `src/pgen/linear_wave_mhd.cpp:164-168` reinterprets `tlim` as a number of
wave *periods*:

```cpp
Real ntlim = lambda / std::abs(ev[wave_flag]) * tlim;
pin->SetReal("parthenon/time", "tlim", ntlim);
```

The entropy mode's eigenvalue **is** the background flow speed, and `lw_mhd.in` sets
`vflow = 0.0`. So `ev[3] = 0` and `tlim = lambda/0 = inf`. The run.logs say it outright —
`tlim=inf nlim=100000` — and the rungs terminated on `nlim` at t = 1406 (N=16) and t = 703
(N=32), with N=64 still running at t = 71.7 after 20396 cycles when the 4 h job wall-clocked.
Those runs reported L1 ~ 1e-12 with **no** convergence (observed order −0.16), which is a real
number — how far a stationary profile drifts over 100k cycles — but not a truncation error.

The rerun (`submit_wp14_entropy.sh`) sets `vflow = 1.0`, giving `ev[3] = 1`, `lambda = 1` (the box
is built so the oblique wavelength is exactly 1.0), hence `tlim = 1.0` = one period. Note the
entropy ladder therefore sits on a **different background state** than the other three
(`vflow = 0`); that is unavoidable, since the entropy mode is identically static without a
background flow, and it is the standard way this family is tested.

Families 4/5/6 (slow+, Alfven+, fast+) are deliberately not run: they are the mirror-image
branches of 2/1/0, exercising the same reconstruction and flux path with the eigenvalue sign
flipped, and add no independent information about the order.

## Scope — what this does and does not establish

**Does:** the base MHD update (PLM reconstruction + HLLD + vl2 + GLM divergence cleaning), as
configured for production, is 2nd-order accurate on a smooth solution.

**Does NOT:** say anything about the order of the *physics packages*. Self-gravity, M1 radiation,
non-ideal diffusion, the tabulated EOS and the chemistry network are all deliberately OFF here.
Several are operator-split and therefore **1st order in time by construction**, so the order of
the full production system is expected to be lower than 2 — establishing the base scheme's order
is the prerequisite for attributing any shortfall to a specific package, not a claim about the
whole.

**Does NOT:** address AMR. This is a uniform grid; prolongation/restriction order at refinement
boundaries is a separate question.

**Feeds WP-22.** The Alfven ladder doubles as a measurement of NUMERICAL RESISTIVITY, at no extra
compute: a wave with zero physical resistivity decays only through the scheme's own dissipation,
so `eta_num ~= 2 (L1/A0)/(k^2 T)`. Over the four rungs that gives **`eta_num ~ dx^1.97`** — the
numerical resistivity vanishes at the scheme's full second order rather than sitting at a floor
refinement cannot beat — with the dimensionless coefficient `eta_num/(v_A dx)` falling
0.0591 → 0.0329 → 0.0165 → 0.0078. Because L1 also contains dispersion error (which does not damp
the wave) this is a rigorous UPPER bound. See `scripts/wp22_eta_numerical.py`.

## Note on running this

The front-end shell is in a **1-CPU cgroup** (`nproc` = 1). The 256x128x128 rung is impractical
there and the load average was already >2 while the pilot ran, so the ladder belongs on the `std`
partition. Account is `banerjee_std` for `std` and `banerjee_gpu` for `gpu`/`gputest` — using the
wrong pairing fails at submit time with "Invalid account or account/partition combination".
