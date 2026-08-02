# WP-11 — tabulated-EOS interpolation error

**Status: RESOLVED — PASS for the production table, with one characterized caveat. 2026-08-02.**
Pure analysis, no simulation, no GPU. Reproduce:

```bash
/beegfs/u/bbg6470/venvs/analysis_env/bin/python docs/validation/scripts/wp11_eos_interp.py
/beegfs/u/bbg6470/venvs/analysis_env/bin/python docs/validation/scripts/wp11_eos_selfconv.py
```

## Correction to the WP-11 specification

`VALIDATION_PLAN.md` describes the production EOS table as "400×1000 (ρ,e)". Two things are
wrong with that as a starting point, and both were found by reading the artifacts:

1. **There are two tables.** `src/eos/eos_table.bin` is 180×220×200; `src/eos/eos_table_hires.bin`
   is 400×1000×920. Only the second is 400×1000.
2. **Production loads the hi-res one** — `runs/prod_flagship_test/fhc_flagship.in:59` and
   `runs/root_ladder/fhc_rootladder.in:77` both set
   `eos_table_file = .../src/eos/eos_table_hires.bin`. The 180×220 table is merely the
   compiled-in default (`src/hydro/hydro.cpp:790`), used only if a deck omits the key.

So the WP splits in two: the coarse table can be measured directly against the hi-res one, but
the *production* table has no finer reference on disk and needs a self-convergence estimate.

**Axis units:** both tables are indexed in **code units**, not cgs — `gen_eos_table.py:139`,
`rho_code = rho_phys/rho0`. The ρ span is log₁₀ρ_code ∈ [−1.738, 18.262], i.e.
**1.0e−20 … 1.0 g/cm³** with `rho0 = 5.467e-19`. (I initially binned the results as though the
axis were cgs, which put the worst errors in the wrong physical regime entirely. The numbers
below are after the conversion.)

## Part 1 — the DEFAULT table (not used in production)

`eos_table.bin` bilinearly interpolated, `eos_table_hires.bin` as truth, coarse-grid nodes
excluded (exact there by construction), N ≈ 4.0e5 off-node interior points:

| table | max | p99 | median | rms |
|---|---|---|---|---|
| `P` | 1.50e−1 | 2.66e−2 | 7.34e−3 | 9.41e−3 |
| `cs2` | 3.77e−1 | 8.34e−2 | 6.72e−3 | 1.96e−2 |
| `logT` | 1.70e−2 | 1.82e−3 | 1.77e−5 | 5.07e−4 |
| `espT` | 1.25e−1 | 2.10e−2 | 2.67e−4 | 4.87e−3 |

**If a deck ever omits `eos_table_file` it silently gets an EOS with 8.3 % p99 / 37.6 % worst-case
error in the sound speed.** That is the practical finding of Part 1: the default is not fit for
production, and nothing in the code warns about it.

## Part 2 — the PRODUCTION table, by self-convergence

No finer table exists, so decimate the hi-res table by stride s = 2,3,4,5, measure each
decimation's error against the hi-res values at off-node points, fit the observed order p in
err ∝ hᵖ, and Richardson-extrapolate back to s = 1.

| table | p(rms) | p(p99) | p(max) | est. production rms | est. p99 | est. max |
|---|---|---|---|---|---|---|
| `P` | 1.91 | 2.01 | 1.61 | 0.152 % | 0.197 % | 1.51 % |
| `cs2` | 1.29 | 1.38 | **0.55** | 0.166 % | 0.509 % | 20.9 % |
| `logT` | 1.73 | 2.10 | 1.46 | 0.0035 % | 0.0072 % | 0.167 % |
| `espT` | 1.92 | 1.95 | 2.00 | 0.028 % | 0.117 % | 0.605 % |

Bilinear interpolation of smooth data is O(h²). `P`, `logT`, `espT` recover p ≈ 1.9–2.1 on the
rms and p99 statistics, so those extrapolations are trustworthy. **`cs2` does not: p(max) = 0.55.**
A max-statistic order that far below 2 means the table is *not smooth on its own grid scale* at
the worst point — which for a protostellar EOS is expected, at the H₂-dissociation and
H-ionization kinks. Consequence: the 20.9 % cs2 max entry is an **extrapolation that is not
justified by its own convergence order** and should be read as "large and not converging",
not as a calibrated number.

## Where the error actually is — the part that decides the WP

A max is meaningless without its location. Errors at stride 2, binned by physical density:

**`cs2`**

| ρ band (g/cm³) | regime | N | max | p99 | rms |
|---|---|---|---|---|---|
| 1e−20 … 1e−13 | envelope → first core | 104361 | 1.19e−1 | 2.73e−2 | 5.44e−3 |
| 1e−13 … 1e−8 | first core → second core | 74900 | 5.35e−2 | 8.51e−3 | 2.08e−3 |
| 1e−8 … 1 | second core / protostar | 118841 | 3.06e−1 | 1.97e−3 | 3.51e−3 |

**`P`** is uniform and benign across all three bands: max ≤ 4.6e−2, p99 ≈ 8e−3, rms ≈ 5.7e−3.

**The global cs2 worst point is at ρ = 0.63 g/cm³, T = 300 000 K** — the extreme top corner of the
table, deep protostellar interior. The FHC target is first-core formation at ρ ≈ 1e−13 and the
flagship reaches ~1e−8 at most, so **the worst error in the table is nowhere near the physics
being simulated.** In the band the collapse actually traverses the p99 is 2.7 % (envelope) and
0.85 % (core), at stride 2 — the production table is finer still.

## Acceptance

**PASS.** WP-11 asks whether the EOS is a leading error term. Against WP-18's realization scatter
σ ≈ 16 % on `MEtor/MEpol`:

- production-table rms error is ≤ 0.17 % on every table;
- p99 ≤ 0.51 %;
- within the collapse-relevant density range the stride-2 p99 is ≤ 2.7 % *before* the
  factor-≈4 refinement to the production grid;
- the one large number (cs2 max ≈ 21 %) is localized to ρ ≈ 0.6 g/cm³ / T = 3e5 K, outside the
  simulated regime.

The EOS interpolation is one to two orders of magnitude below the noise floor that WP-18 set.

## Not claimed

- That the *table generator* is correct. WP-11 measures interpolation error only — how faithfully
  the interpolant reproduces the tabulated values. Whether the tabulated multi-Saha physics is
  right is a different question and is not tested here.
- That p ≈ 2 holds band-by-band. The orders above are fitted globally; the bands were compared at
  fixed stride 2, not refit per band.
- Anything about the opacity table (`src/radiation/opacity_table.bin`), which is a separate
  interpolant and is not covered by this WP.

## Follow-up worth doing

Make the missing-key case loud: `hydro.cpp:790` silently falls back to the 180×220 table. Given
Part 1, that fallback should either warn or be removed, so a deck typo cannot quietly degrade the
sound speed by several percent.
