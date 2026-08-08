# build_gpu_v5 gate — first GPU compile of the 2026-08-05 audit batch

**Status: the INERT claim PASSES. The other two legs were badly designed by me and prove nothing;
A1 and N2/N3 remain UNVERIFIED ON GPU.** Recorded in full because two of the three legs are
instructive failures.

| item | value |
|---|---|
| binary | `build_gpu_v5/bin/athenaPK` md5 `6d5b9895` |
| source | **`9f406ce` PLUS one uncommitted edit — see the provenance warning below. NOT a clean commit.** |
| build | job 2460514, `BUILD_DONE exit=0`, 171 objects, ~6 min on `std` |
| gate | job 2460585 (6 legs) |
| control | job 2460695 (v4 against itself) — **the leg that decided the verdict** |
| baseline | `build_gpu_v4/bin/athenaPK` md5 `869c1d34` (production) |

The build was verified genuine rather than assumed: all 171 objects freshly compiled, and
`d1_meminfo` / `ct_2d_outofplane_checked` are present in v5 and absent in v4 (`strings`).

## PROVENANCE WARNING — v5 corresponds to NO commit

**A second session was editing `src/` while this build and gate ran.** The tree was clean when
the build was submitted; by the time the gate finished, `src/radiation/radiation.cpp` (mtime
16:11:23) and `src/hydro/hydro.cpp` (16:24:19) were both modified and uncommitted. Established by
direct check, not inference:

- `git show 9f406ce:src/radiation/radiation.cpp | grep -c "audit N13"` → **0**. The N13 opacity
  guard is **not in the commit**; it is a working-tree edit.
- `strings build_gpu_v5/bin/athenaPK | grep -c "audit N13"` → **1** (v4 → 0). The build picked it
  up, because the edit landed at 16:11:23 and the build compiled `src/` at ~16:15.
- It **fired during the gate**: `nodust_new/run.log` has exactly one more line than
  `nodust_old` (244 vs 243), and that line is the N13 warning.

So **`build_gpu_v5` = `9f406ce` + the uncommitted N13 guard.** The N14 EOS guard (`hydro.cpp`,
16:24:19) landed *after* the build and is **not** in v5.

The N13 guard is warn-only, so the gate result stands — but the binary is not reconstructable
from git, which is exactly the condition Gate F exists to prevent. **Do not promote v5 to
production until the tree is committed and, ideally, rebuilt from a clean commit.**

Its warning also independently confirms N13 at runtime:
`opacity_unit = 0.0153623` (generator) vs `0.0153415` (PhysUnits), **ratio 1.00135, κ 0.135079 %
off** — matching the 1.001351 / +0.1351 % re-derived by hand from `be_normalization.hpp`.

## The result that matters: the audit batch is inert on the production path

`nodust_old` (v4) vs `nodust_new` (v5), 128³, `nlim=40`, dust off, production physics otherwise:
**40 of 44 history columns bit-identical over 41 rows, including every state column** (mass,
momenta, KE, tot-E, ME, relDivB, and all the cap/mag/cons diagnostics).

Four columns differed — `cons-Poutx/y/z` and `cons-Mout-out` — and the raw numbers looked too big
to wave away: `cons-Poutx` row 3 differed by 6.7e-10 on a value of 6.55e-06, i.e. **~1e-4
relative, ~67 units in the last of the 6 printed significant figures.** The v4 gate had dismissed
the same class of thing as "one unit in the last printed digit"; at 67 units that explanation does
not hold, and accepting it would have been assuming the conclusion.

### The control settles it

Job 2460695 ran **v4 twice**, byte-identical invocation:

| column | control (v4 vs v4) | gate (v4 vs v5) |
|---|---|---|
| `cons-Poutx` | **23** rows, absmax 1.0e-09 | 16 rows, absmax 1.0e-09 |
| `cons-Mout-out` | 1 row, absmax 1.0e-05 | 1 row, absmax 1.0e-05 |
| `cons-Poutz` | 6 rows, absmax 1.0e-09 | 7 rows, absmax 1.0e-09 |
| `cons-Pouty` | 9 rows, absmax 1.0e-09 | 7 rows, absmax 1.0e-09 |

**The same binary does not reproduce itself on these four columns — and differs in MORE rows than
v5 differs from v4.** `cons-Pout*` are surface-integral reductions; GPU reduction/atomic ordering
is not reproducible between runs. So the difference is run-to-run nondeterminism in the
diagnostics, not a code change.

**Correct statement of the result:** *v5 is bit-identical to v4 on every column that is itself
reproducible.* A bare `array_equal` gate cannot express that, which is why it returned FAIL.

**Correction to the v4 gate record:** it called this a *printing artefact*. That is wrong — it is
genuine run-to-run nondeterminism that merely surfaces at the last printed digit. Same verdict,
wrong mechanism, and the wrong mechanism predicts a bound (1 unit) that this run violates.

## Leg 2 — dust (N2/N3): INCONCLUSIVE, my design error

Expected to DIFFER, since N2 (dust sublimation T from the table) and N3 (dust cgs scales) are
*intended* to change numbers. It "differed" — in the same four nondeterministic columns plus
`maxRelDivB`, all at ~1e-6. **Nothing dust-related changed, and nothing could have:**

- the history file has **no dust and no temperature column** (44 columns, all hydro/MHD/diagnostic);
- N2 only bites at the 1500 K sublimation threshold, and this leg sits at T ≈ 10 K, t ≈ 0.17.

So the configuration cannot express the fix. **The script reports PASS here; that PASS is false
and must not be quoted.** A real GPU test of N2/N3 is the existing CPU execution decks
`runs/audit_fix_regress/{H_dust_ideal,I_dust_hydrogen}.in`, or a deep restart hot enough to
sublimate.

## Leg 3 — A1 refinement: NULL TEST, my design error

`refinement=adaptive`, `numlevel=3`, `njeans=8`, `nlim=25`. `nbtotal` is a **constant 64 in both
legs**, so the verdict script's "nbtotal IDENTICAL ⇒ A1 is inert here" is not evidence of
anything: **the mesh never refined.** `fhc_rootladder.in`'s own header puts the njeans=8 trigger
at t ≈ 0.9994 and these legs stopped at **t = 0.365**.

This is exactly the trap already recorded for `A_multipole.in`/`B_swindle.in` — *"numlevel pins
nbtotal at 64 for every cycle, so 'identical' there is a null result, not evidence"* — and it was
walked into again. **A1 is untested on GPU.** The cheap correct design already exists:
`F_a1discrim.in`, with `njeans = 4.7` chosen to sit *between* the ideal criterion (min n_J 4.537)
and the tabulated one (4.925), so the two formulas provably disagree — 120 blocks vs 64.

## Performance

No regression. Mean `wsec_step` over the last 8 cycles: nodust 1.914 (v4) → **1.936** (v5),
dust 1.920 → **1.929**, i.e. ~1 %. *(An earlier reading of "~14 % slower" was taken from an
in-flight leg during startup transients and is withdrawn — measure completed legs only.)*

## Bottom line

- **v5 is safe to promote on the evidence that exists**: the ~28 changes claimed inert are inert
  on the exact path production takes.
- **Two of the batch's substantive fixes (A1, N2/N3) still have no GPU evidence.** They have CPU
  gates, and A1 is argued inert today on physics grounds (below 85 K the table gives Γ₁ ≈ 5/3,
  cs_ideal/cs_table = 0.917 ⇒ over-refinement, safe). That argument is sound but it is not a GPU
  measurement.
