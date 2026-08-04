# WP-9 — Multigroup radiation: is the production group count converged?

**Status: PASS.** Production's `n_group = 3` is converged. Doubling to 6 groups changes nothing
above **3.4e-06** relative in any history column. The gray approximation would have cost
**−0.0108 %** on `MEtor/MEpol`.

Job 2452599 (2 GPUs, `build_gpu_v3` = `6b1fe753`), 3 legs × 90 cycles, 128³ uniform,
deck `runs/root_ladder/fhc_rootladder.in`.

## A correction to the campaign record first

`CAMPAIGN_STATUS` and the validation plan both described production as *gray, with multigroup
untested*. **That is wrong.** `runs/root_ladder/fhc_rootladder.in:192` sets `n_group = 3`.
Production has been multigroup all along. WP-9 is therefore not "does multigroup change the
answer" but the sharper question **is 3 groups enough** — plus, as a separate and separately
useful number, what the gray answer would have been.

## Result, read at t = 0.90

Never at t = 1.0: WP-7 established that the uniform-grid (`refinement = none`) endpoint is a
singular stall where dt falls 1–2 decades in the last 0.01 t₀.

| `n_group` | `MEtor/MEpol` | vs production | ME | KE |
|---|---|---|---|---|
| 1 (gray) | 7.471130e-03 | **−0.0108 %** | 5.00634e+01 | 2.42424e+03 |
| **3 (production)** | 7.471938e-03 | — | 5.00650e+01 | 2.42546e+03 |
| 6 | 7.471938e-03 | +0.0000 % | 5.00650e+01 | 2.42546e+03 |

Acceptance was 3-vs-6 within σ = 16 %. Measured: below the 6 printed digits.

## The "+0.0000 %" was verified, not reported

An exact zero is the signature of a **null test** — a parameter that never reached the code —
and this campaign has produced several (WP-1's STS cap peaking at 4.04 so no cap in 60–1000 ever
bound; the entropy-wave family running with `tlim` reinterpreted as wave periods; the first B1
blast falsification never reaching a boundary; WP-17 attempt 1, below). So `n_group = 6` was
checked for engagement rather than assumed:

```
cmp ng3 ng6 -> files DIFFER
column-by-column: 2 of 31 columns differ
   relDivB        3.3910e-06
   mag-dissO      1.7644e-06
```

The knob engaged. The effect is simply smaller than the headline column's printed precision.
That is a convergence result, not an inert parameter.

## The 6-group opacity table had to be built, and had to nest

`n_group` is **not** a free runtime knob under `opacity_model = tabulated`. Attempt 1 (job
2452594) ran groups 1 and 3 and then aborted:

```
what():  GroupOpacityTable: file n_group mismatch (3 != 6)
```

The shipped `opacity_table.bin` encodes the production 3-group structure. A 6-group table was
generated with the tree's own `src/radiation/gen_opacity_table.py`, **keeping both production
edges and log-splitting the finite decades between them**:

```
production (3):  0, 1e12,                                   1e15, 1e30
new        (6):  0, 1e12, 5.62341e12, 3.16228e13, 1.77828e14, 1e15, 1e30
```

The nesting is what makes this a convergence test: the 3-group partition is a strict subset of
the 6-group one, so 3→6 is a genuine refinement of the *same* spectral partition. An arbitrary
6-way binning would compare two unrelated partitions and any difference would be uninterpretable.

Written to `src/radiation/opacity_table_g6.bin`. `n_group = 1` runs against the 3-group file:
`radiation.cpp:156` documents a single group spanning [0,∞) as the gray path, and the loader does
not enforce the match on that branch (confirmed — the leg completed).

## What this does and does not license

- **Does**: production `n_group = 3` needs no change; the spectral discretization is not a
  contributor to the flagship error budget.
- **Does**: comparisons against gray literature runs can quote −0.0108 % as the spectral cost on
  `MEtor/MEpol` for *this* configuration and phase.
- **Does not**: say anything past t = 0.90 or at first-core densities, where the opacity is
  decades higher and the group structure could matter more. The legs are 90 cycles of the
  pre-collapse phase on a uniform grid.
