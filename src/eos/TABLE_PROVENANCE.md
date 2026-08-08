# EOS / opacity table provenance — what production depends on, and how to rebuild it

Written 2026-08-08 because of a reproducibility hole found while closing audit item N14: **the
flagship's EOS table is not in version control, and nothing recorded how to regenerate it.**

`runs/root_ladder/fhc_rootladder.in:91` and all 24 `runs/ensemble/design01/point*/fhc_ens.in`
name `src/eos/eos_table_hires_v2.bin`. That file is 12.5 MB and `src/eos/.gitignore` — correctly,
and deliberately — excludes the hi-res tables as "regenerable via gen_eos_table.py and too large
to version". The policy is fine. What was missing is the *record* that makes "regenerable" true in
practice: dimensions, hashes, and the exact command. Without it, "regenerable" is an assertion, not
a procedure.

## The tables

| file | tracked | ver | nr | ne | nT | md5 | used by |
|---|---|---|---|---|---|---|---|
| `eos_table.bin` | yes | 1 (legacy) | 180 | 220 | 200 | — | nothing; superseded |
| `eos_table_v2.bin` | yes | 2 (cgs) | 180 | 220 | 200 | — | the **default** since 2026-08-08 |
| `eos_table_hires.bin` | **no** (ignored) | 1 (legacy) | 400 | ? | ? | `6b8e3999eca19806d8f4d43054e0447c` | superseded |
| `eos_table_hires_v2.bin` | **no** (ignored) | 2 (cgs) | 400 | 1000 | 920 | `80746eb5924487f2ab350ea050c8bbf4` | **the flagship + all 24 ensemble members** |
| `../radiation/opacity_table.bin` | yes | 1 (legacy) | — | — | — | — | 28 completed decks, deliberately |
| `../radiation/opacity_table_v2.bin` | yes | 2 (cgs) | — | — | — | — | flagship; the **default** since 2026-08-08 |

Version is the first `int64` of the file: `0x454F535441424C31` ("EOSTABL1") for v2; a legacy v1
file starts with `nr` instead. `flags` bit 0 set = arrays stored in **cgs**, converted at load by
the run's own units.

## Rebuild

```bash
/beegfs/u/bbg6470/venvs/analysis_env/bin/python src/eos/gen_eos_table.py \
    src/eos/eos_table_hires_v2.bin 400 1000 920
```

`build_table(path, nr, ne, nT)` is a deterministic physics computation (multi-Saha + H2
dissociation on a fixed log grid), so the same generator at the same commit reproduces the same
bytes. **Check the md5 above after regenerating.** If it does not match, the generator changed and
the flagship's EOS changed with it — that is a result-changing event, not a build detail.

## Why v1 vs v2 matters (audit N13/N14)

A v1 table stores axes and arrays in the generator's own **rounded** code-unit constants; a v2
table stores cgs and lets the loader convert with the run's exact normalization. Measured:

* opacity: v1 is uniformly **+0.13538 %** high (all of kP/kR/ks, every one of 7200 nodes).
* EOS: the runtime warning reports `rho0 = 5.467e-19` vs the run's `rho_unit = 5.46683e-19`,
  **rel diff 3.15e-05**, so a v1 table is consulted at slightly wrong `(rho, esp)`.

Both loaders read v1 files correctly and **warn** (verified firing, with the numbers above, and
verified silent on v2), so completed decks that name a v1 file explicitly stay reproducible and
are left alone on purpose. Only the *defaults* were changed, which closes the silent path where a
deck that omits the key inherits a biased table with nothing in the deck to show it.
