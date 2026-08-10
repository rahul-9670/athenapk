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

---

## 2026-08-10 — the opacity table was 2 decades short of the EOS, silently (item A3)

**The defect.** `OpacityTable::bilin` (`radiation/radiation_opacity.hpp:62-77`) clamps *both* the
grid index and the interpolation weight. A lookup past an edge therefore returns the edge opacity:
no warning, no error, no NaN, and nothing in the output to show it happened. Measured domains:

| table | rho [g/cm^3] | T [K] | grid |
|---|---|---|---|
| `eos_table_hires_v2.bin` (flagship EOS) | 1e-20 … **1e+00** | 8 … **3.000e5** | 400 x 1000 x 920 |
| `opacity_table_v2.bin` (flagship opacity) | 1e-20 … **1e-02** | 3 … **1.000e5** | nr=60, nT=120 |

So the EOS stayed valid two decades in density and half a decade in temperature past the point
where the opacity quietly froze. `rhocrit = 1e-13`, so the density edge sits at ~1e11 x rhocrit —
**the second core forms right there**, which makes this a hard ceiling on the second-collapse
campaign rather than a corner case. The temperature edge is the broader of the two: it binds at
*every* density, not only second-core densities.

**Measured cost of the clamp** (`bilin` ported to python exactly as written, compared against
`gen_opacity_table.group_opacities` evaluated directly — script kept with the verification below):

| rho [g/cm^3] | T [K] | exact kappa_R | v2 (clamped) | error |
|---|---|---|---|---|
| 1e+00 | 2.0e3 | 8.000e1 | 3.741e0 | **-95.3 %** |
| 1e+00 | 1.0e4 | 1.000e4 | 2.315e3 | -76.8 % |
| 1e-02 | 2.0e5 | 8.385e4 | 4.743e5 | **+466 %** |
| 1e-06 | 2.0e5 | 8.385e0 | 4.841e1 | **+477 %** |
| 1e+00 | 3.0e5 | 3.043e6 | 4.743e5 | -84.4 % |

### The new table

```bash
/beegfs/u/bbg6470/venvs/analysis_env/bin/python src/radiation/gen_opacity_table.py \
    --out src/radiation/opacity_table_v3.bin \
    --nr 101 --nT 201 --rho_min 1e-20 --rho_max 1e0 --T_min 3.0 --T_max 3.0e5 \
    --edges 0,1e12,1e15,1e30
```

| file | ver | ng | nr | nT | rho | T | md5 |
|---|---|---|---|---|---|---|---|
| `opacity_table_v2.bin` | 2 (cgs) | 3 | 60 | 120 | 1e-20…1e-2 | 3…1e5 | `fa5cb4aee65e962dd3a1237a721d75a7` |
| `opacity_table_v3.bin` | 2 (cgs) | 3 | 101 | 201 | 1e-20…**1e0** | 3…**3e5** | `c8c7c7c91d064b9b509dbfc5926cadf6` |

`--edges 0,1e12,1e15,1e30` is not a free choice: it is what `BuildRadGroups` produces for the
flagship's `n_group = 3` with the default `nu_min_hz = 1e12` / `nu_max_hz = 1e15`
(`radiation/radiation_groups.hpp:338-353`). A table generated with different edges is silently
wrong for the deck that reads it. 497 KB, so it is **tracked**, like `opacity_table_v2.bin` — the
gitignore policy covers only the 12.5 MB EOS hires tables.

### Verification (not assumed — measured)

1. **Coverage:** v3 covers the full `eos_table_hires_v2.bin` domain in both axes; v2 covers
   neither. Confirmed at runtime by the new startup report, against the flagship EOS file.
2. **Accuracy did not regress to buy coverage.** 4000 random points inside *v2's own* domain,
   both tables' interpolants compared against the exact generator function:

   | quantity | v2 median | v3 median | v2 p95 | v3 p95 |
   |---|---|---|---|---|
   | kappa_P | 3.483e-3 | **1.455e-3** | 2.060e-1 | **1.090e-1** |
   | kappa_R | 3.387e-3 | **1.443e-3** | 2.060e-1 | **1.043e-1** |

   v3 is **2.4x more accurate** where they overlap (dlogrho 0.305 -> 0.200, dlogT 0.038 -> 0.025),
   so the extension is strictly a gain. The residual p95 ~10% and max ~50% sit on the dust
   sublimation cliff (Bell & Lin regime 3, kappa ~ T^-24), where bilinear interpolation on any
   practical log grid cannot track the function; that limitation is **pre-existing and unchanged
   in kind** — v2's max error on the same sample was 109%, v3's is 51%.
3. **Sanity:** all 5 arrays finite, non-negative; `ks` takes only {0, 0.348}; the Bell & Lin
   anchor `kappa_R,gray == kappa_BL` holds at 10/100/1000/1e4 K by construction.

### The clamp is no longer silent

`radiation/radiation.cpp` now prints the loaded opacity table's domain at startup on rank 0, and
warns explicitly when it does not cover the EOS table's domain, naming the shortfall in decades.
Verified firing on both branches with the flagship EOS table:

```
### Radiation opacity table: .../opacity_table_v2.bin
###   domain rho[g/cm^3] 1e-20 .. 0.01  (nr = 60, dlogrho = 0.305085)
###   domain T[K]        3 .. 100000  (nT = 120, dlogT = 0.0380074)
###   lookups outside this box are CLAMPED to the edge value.
### WARNING Radiation (A3): the opacity table does NOT cover the thermodynamic range the EOS
###   rho: opacity ends at 0.01 g/cm^3, EOS reaches 1 g/cm^3 (short by 2 decades)
###   T:   opacity ends at 100000 K, EOS reaches 300000 K (short by 0.477121 decades)
```
```
### Radiation opacity table: .../opacity_table_v3.bin
###   domain rho[g/cm^3] 1e-20 .. 1  (nr = 101, dlogrho = 0.2)
###   domain T[K]        3 .. 300000  (nT = 201, dlogT = 0.025)
###   covers the EOS table domain (rho <= 1 g/cm^3, T <= 300000 K): OK.
```

It **warns, never aborts**: a hard failure here would be restart-hostile for every completed deck
that legitimately names a narrower table (the lesson of audit N3).

### Residual, not fixed here

`radiation/coupling_tmax_K` defaults to `opacity_table_Tmax_K` = **1e6 K**
(`radiation.cpp:263,321`), which is an independent deck key and is *not* read from the table file.
So the multigroup rtsafe bracket still reaches above even the v3 ceiling of 3e5 K. Swapping the
table does not change it (verified: no coupling side-effect), but a run that genuinely gets to
1e6 K would be extrapolating the bracket past both tables.

### Adoption: NONE yet — deliberately

**No deck points at `opacity_table_v3.bin`.** Swapping the flagship's opacity table is
result-changing and is a user decision, so v3 is staged and recorded, not adopted. `v2` remains
the default in `radiation.cpp:185-187` and remains what every production deck names, so all
completed results stay reproducible. Adopting it means editing
`runs/root_ladder/fhc_rootladder.in:203` **and** all 24 `runs/ensemble/design01/point*/fhc_ens.in`
— check the adoption count afterwards, the way the 2026-08-08 default flips had to be.
