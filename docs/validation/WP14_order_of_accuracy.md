# WP-14 — order of accuracy on a smooth problem, in the production configuration

**Status: fast-magnetosonic family PASSES (p = 2.06). Remaining families queued (job 2448312).**
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

## Result — fast-magnetosonic (`wave_flag = 0`)

| N | RMS L1 | rho | M1 | E | B2c | p(RMS) | p(rho) | p(M1) | p(E) | p(B2c) |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 1.897e-07 | 4.402e-08 | 6.618e-08 | 1.513e-07 | 1.893e-08 | — | — | — | — | — |
| 32 | 5.588e-08 | 1.026e-08 | 1.796e-08 | 4.513e-08 | 4.974e-09 | 1.76 | 2.10 | 1.88 | 1.75 | 1.93 |
| 64 | 1.343e-08 | 2.356e-09 | 4.632e-09 | 1.078e-08 | 1.184e-09 | **2.06** | 2.12 | 1.96 | 2.07 | 2.07 |

**PASS.** The observed order approaches 2 from below and reaches **2.06** between the two finest
completed rungs, with every individual conserved variable in 1.96–2.12. The trend is monotone
toward 2, which is the signature of a genuinely 2nd-order scheme rather than a coincidence at one
resolution pair.

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

**Not yet measured:** Alfvén (1), slow (2) and entropy (3) families — queued as job 2448312. The
Alfvén family is the one to watch: it is the mode most sensitive to numerical resistivity, and it
doubles as the input to WP-22.

## Note on running this

The front-end shell is in a **1-CPU cgroup** (`nproc` = 1). The 256x128x128 rung is impractical
there and the load average was already >2 while the pilot ran, so the ladder belongs on the `std`
partition. Account is `banerjee_std` for `std` and `banerjee_gpu` for `gpu`/`gputest` — using the
wrong pairing fails at submit time with "Invalid account or account/partition combination".
