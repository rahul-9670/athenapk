# WP-8 follow-up — why `mag-dissO` / `mag-dissA` do not converge, and what was changed

**Status: DIAGNOSED (root cause identified, two competing hypotheses falsified) + FIX
IMPLEMENTED & GATED. 2026-07-31.**

## The observation

On the njeans ladder at matched epoch ρ = 1e-12 g/cm³:

| quantity | nj4 → nj8 | nj8 → nj16 | behaviour |
|---|---|---|---|
| `mag-ME` | −4.2% | **+0.2%** | converged |
| `mag-Jsq` | −50.0% | −39.6% | monotone, NOT converged |
| `mag-dissO` | −76.5% | **+40.5%** | not even monotone |
| `mag-dissA` | −76.3% | −72.9% | monotone, NOT converged |
| `dissA/dissO` | +0.8% | −80.7% | meaningless |

Magnetic energy converges to 0.2% while the dissipation integrals move by tens of percent
with no sign of settling. So the field is converging; the dissipation *diagnostic* is not.

## Two hypotheses, both FALSIFIED

**H1 — the diagnostic's curl is not the scheme's curl.** `mag_diag.cpp` takes J = ∇×B with a
**2Δx cell-centered** central difference; the flux kernel (`resistivity.cpp:139`) differences
across **one** cell at a face. A 2Δx central difference is blind to the Nyquist mode, so it
could systematically miss grid-scale current.

*Test:* both stencils were evaluated on identical saved snapshots (block interiors only,
since phdf carries no ghosts) at each rung.

| leg | Jsq (2Δx centered) | Jsq (1-cell) | ratio |
|---|---|---|---|
| nj4 | 1.91299e4 | 1.94129e4 | 1.015 |
| nj8 | 9.11550e3 | 9.18093e3 | 1.007 |
| nj16 | 5.53743e3 | 5.54524e3 | 1.001 |

Non-convergence is identical either way (−52.3%/−39.3% vs −52.7%/−39.6%). **The stencil is
fine**, and the ratio → 1 with refinement shows the field is getting *smoother* relative to
the grid, not rougher. H1 rejected.

**H2 — stale `eta` from the `rkl2_freeze_eta` cache.** `Jsq` carries no `eta`; `dissO` does.
If the cache were inconsistent with `prim`, `dissO` would jitter where `Jsq` does not.

*Test:* `dissO` moves in the **opposite** direction to `Jsq` in only **3–7%** of history
rows, and the row-to-row jitter of both falls together with resolution. H2 rejected.

## Root cause — CONCENTRATION

The integrals are carried by a vanishing, resolution-dependent fraction of the domain.
Volume fraction supplying 90% of the integral, measured at matched epoch:

| leg | 90% of ∫\|J\|²dV | 90% of ∫ρ\|J\|²dV (proxy for ∫η_O\|J\|²dV) |
|---|---|---|
| nj4 | 7.57e-7 | **3.55e-8** (6,085 cells of 3.18e7) |
| nj8 | 1.56e-6 | **7.68e-8** (103,030 of 3.56e7) |
| nj16 | 4.67e-6 | **1.26e-7** (575,345 of 3.71e7) |

(ρ is used as the weight because η_O climbs steeply with density, so ∫η_O J²dV is even more
concentrated than ∫ρJ²dV.)

`mag-dissO` and `mag-dissA` are therefore **not volume integrals in any useful sense — they
are point samples of the innermost core**, taken in precisely the region whose resolution
changes between rungs. Each refinement resolves a previously-unresolved region and changes
the integrand there by orders of magnitude. A single scalar summing an integrand that spans
~7 decades across the domain cannot converge until that ~1e-7 of the volume is itself
converged. Corroborating evidence: the effective volume-weighted ⟨η_O⟩ = dissO/Jsq spans
**6.7 / 4.3 / 3.0 decades** over the run, narrowing with resolution exactly as the carrying
region becomes better resolved.

**This is not a code defect.** The kernel computes ∫η|J|²dV correctly, non-negative,
volume-weighted, mesh-wide. The defect is in treating an extremely concentrated integral as
a convergent global scalar.

### The ratio is worse than either part

`dissA/dissO` is a quotient of two integrals dominated by **different regions** — ambipolar
in the diffuse envelope, Ohmic in the dense core — so it is not a physical quantity at all.
It is additionally undefined on the first history row: at t = 0 the field is uniform B0z, so
J ≡ 0 and dissO = 0 exactly (confirmed: exactly 1 row per leg has dissO = 0). The
qualitative statement "ambipolar dominates Ohmic by 4–5 decades" is robust. **The number is
not, and must never be quoted as a measurement.**

## What was changed

All new columns are behind `hydro/mag_diag_rho_split` (code density, **default 0 = OFF**).
Default keeps the original column set byte-for-byte — deliberately, because appending
columns would shift `maxRelDivB` from index 31 to 38 and silently break every existing
analysis script and any comparison with the ladder history files already on disk.

| column | meaning |
|---|---|
| `mag-dissO-hi` / `mag-dissO-lo` | Ohmic dissipation above / below the density split |
| `mag-dissA-hi` / `mag-dissA-lo` | ambipolar dissipation above / below the split |
| `mag-Vhi` | volume with ρ > split — how much volume the "hi" budget covers |
| `mag-dissOsq` / `mag-dissAsq` | ∫q² dV, for the concentration measure below |

Each split piece is an integral over a region defined by **physics**, not by the grid, so
core and envelope budgets converge (or visibly fail to) independently.

**Concentration measure.** From ∫q dV and ∫q² dV, analysis forms the inverse participation
ratio

> f_eff = (∫q dV)² / (V_box · ∫q² dV)

= the fraction of the box that would carry the whole integral if q were uniform over it.
f_eff ≈ 1 means volume-filling; f_eff ≈ 1e-7 means the number is a point sample and will not
converge. This makes the pathology impossible to miss in future output.

*A literal "volume where q > 0" was implemented first and is useless* — η and J are nonzero
almost everywhere, so it returns the box volume (measured: 1.40608e5 = exactly 52³ on the
L=52 smoke deck). Replaced with the ∫q² form.

## Guidance

1. Use `mag-dissO-hi/lo`, **not** the global `mag-dissO`, for any convergence work.
2. Report `f_eff` alongside any dissipation number. If f_eff ≲ 1e-3, say the number is a
   core point sample, not a global budget.
3. Never quote `dissA/dissO`. `analyze.py` deliberately no longer offers it.
4. The remaining non-convergence of the core budget is **WP-22** (numerical vs physical
   resistivity in the first core), not a diagnostic problem: it is the real physics question
   of whether the fossil-field result is grid-determined.

## Analysis-script hardening (same class of bug)

`analyze.py` now parses hst column indices **from the header** instead of hardcoding them —
without this, enabling the split would shift `maxRelDivB` and an index-based reader would
keep running while silently reporting the wrong quantity under the right name.

## ON-state verification (gate deck `runs/eos_smoke/fhc.in`, L=52, 12 cycles, rho_split=1.0)

| check | result |
|---|---|
| OFF state (`rho_split=0`) vs the pre-change build | **byte-identical**; columns still `[30]=mag-dissAcap [31]=maxRelDivB` |
| `mag-dissO` = hi + lo | 9.402320e-10 vs 9.402318e-10 — agree to **7 significant figures**, which is the history file's *text* precision, not a computational discrepancy (with `rho_split=100`, where hi is exactly 0, the two matched bit-for-bit) |
| `mag-dissA` = hi + lo | 4.901630e-02 vs 4.901633e-02, same 7-figure agreement |

**The split immediately confirms why the ratio was meaningless:**

- Ohmic: **89.60%** of `mag-dissO` comes from ρ > 1 code — the core.
- Ambipolar: only **16.13%** of `mag-dissA` does — it is envelope-dominated.
- `mag-Vhi` = 218.84 code = **0.156% of the box** carries that 89.6% of the Ohmic budget.

So numerator and denominator of `dissA/dissO` genuinely live in different places, exactly as
predicted. Concentration measure on the same row:

| | f_eff | verdict |
|---|---|---|
| `mag-dissO` | 1.88e-4 | concentrated |
| `mag-dissA` | 7.87e-3 | concentrated (40× less so than Ohmic) |

Both ≪ 1, and Ohmic is the more concentrated of the two — consistent with the core/envelope
split above. Note this is the *coarse smoke deck*; at production resolution the carrying
fraction was measured at ~1e-7, i.e. three to four decades worse.

---

# 2026-08-06 — the density split was never tested on the ladder. It is now, and for `Jsq` it does NOT work.

**Status: the remedy above is PARTIALLY FALSIFIED. Splitting by density does not restore
convergence for `Jsq`, because `Jsq`'s concentration is not organised by density.**

## Why this was reopened

Everything above — root cause, the two falsified hypotheses, the split columns, `f_eff` — was
implemented, gated and documented, but the split was verified only for *self-consistency* on a
12-cycle L = 52 smoke deck (`hi + lo == global` to the history file's text precision). **It was
never applied to the njeans ladder**, i.e. to the measurement that showed the non-convergence in
the first place. A remedy that is implemented and documented but never shown to work on the
failing case is not a closed finding.

## Method

`docs/validation/scripts/wp8_split_convergence.py`, run at the same matched epoch as the original
table (ρ_max = 1e-12 g cm⁻³), split at ρ_crit = 1e-13 g cm⁻³ — a physical boundary, not a tuned
one. Re-running the ladder with `mag_diag_rho_split` enabled would cost GPU time this project does
not have spare; all 135 snapshots (986 GB) are on disk, so the split is recomputed offline instead.

**Why `Jsq` and not `dissO`/`dissA`.** The ladder snapshots carry only `prim` — there is no
`nonideal_eta` field (that was added later, for WP-22 part 3, on a different run), and
reconstructing η_A offline is only an *upper bound* because the equilibrium NICIL/Wardle ceiling
needs the full grain + Saha charge solve. `Jsq = ∫|J|²dV` needs only **B** and the grid, so it is
exact offline — and it is itself one of the non-converging quantities. J is taken by the 1-cell
face difference on block interiors only (phdf carries no ghosts, so 17.6 % of cells — the outer
layer of each 32³ block — are dropped rather than differenced against wrong data).

**Validation of the script itself:** it reproduces the published global numbers exactly —
**−52.3 % / −39.3 %**, against −52.3 %/−39.3 % recorded above for the 1-cell stencil.

## Result

| | nj4 | nj8 | nj16 |
|---|---:|---:|---:|
| `Jsq-glob` f_eff | 1.137e-7 | 1.384e-7 | 3.102e-7 |
| `Jsq-hi` (ρ>ρ_crit) share of global | 2.58 % | 1.03 % | 2.14 % |
| `Jsq-hi` f_eff | **0.575** | **0.465** | **0.642** |
| `Jsq-lo` (ρ<ρ_crit) share of global | 97.42 % | 98.97 % | 97.86 % |
| `Jsq-lo` f_eff | 1.111e-7 | 1.369e-7 | 3.303e-7 |
| V(ρ>split)/V | 4.53e-9 | 3.34e-9 | 2.20e-9 |

| quantity | nj4→nj8 | nj8→nj16 | verdict |
|---|---:|---:|---|
| `Jsq-glob` | −52.3 % | −39.3 % | monotone, not converged |
| `Jsq-hi` | −80.9 % | +25.6 % | **NOT MONOTONE** — no convergence claim possible |
| `Jsq-lo` | −51.6 % | −39.9 % | monotone, not converged |

## What this means

1. **The split does not rescue `Jsq`.** The low-density bin carries 97–99 % of the integral *and*
   keeps f_eff ~ 1e-7, so it simply *is* the original pathology under a new name. Its convergence
   history (−51.6 %, −39.9 %) is the global one to within a percent.
2. **The core bin is nevertheless a well-posed integral.** `Jsq-hi` has **f_eff ≈ 0.47–0.64** — it
   is volume-filling *within its own region*, exactly what the split was designed to produce. It
   still does not converge (it is not even monotone), but that is now a statement about the
   physics being under-resolved, not about the diagnostic being a point sample. That is a
   qualitatively different and much more tractable failure.
3. **The premise "the concentration is in the core" is wrong for `Jsq`.** The smoke-deck
   measurement above found 89.6 % of `dissO` above ρ = 1 code, and the split was designed around
   that. But `dissO` is weighted by η_O, which climbs steeply with density; `Jsq` is unweighted.
   At production resolution the unweighted current is **envelope**-dominated and its pathological
   concentration lives *below* ρ_crit — presumably in grid-scale current sheets, which a density
   threshold cannot separate by construction.

**Consequence:** a density split is the right tool for `dissO` (core-dominated by η weighting) and
the wrong tool for `Jsq`. `Jsq` needs a split on a current-sheet indicator, not on density — or it
should be retired as a convergence metric in favour of the split dissipation budgets.

## What is still open

- The same test **for `dissO`/`dissA`**, which is what the split was actually built for. It needs
  either a ladder re-run with `mag_diag_rho_split` enabled and `nonideal_eta` in the output list,
  or the ladder snapshots regenerated with that field. **Until then the split's core claim is
  supported for `Jsq`'s core bin (f_eff → 0.5) and refuted for its envelope bin, and untested for
  the dissipation integrals themselves.**
- Guidance item 1 above ("use `mag-dissO-hi/lo`, not the global") remains sound but is now known to
  be *insufficient on its own*: report `f_eff` per bin, and treat any bin with f_eff ≲ 1e-3 as a
  point sample regardless of which side of the split it is on.

---

# 2026-08-08 — the current-sheet split, implemented

**Status: the remedy the 2026-08-06 reopening ASKED FOR is now in the code. It is a diagnostic
addition, not a physics change; the OFF state is bit-identical.**

The reopening ended with a specific instruction: *"`Jsq` needs a split on a current-sheet
indicator, not on density — or it should be retired as a convergence metric."* Density cannot work
for `Jsq` by construction, because `Jsq`'s pathological concentration sits **below** the density
threshold (the low-density bin carries 97–99 % of the integral at f_eff ~ 1e-7, i.e. `Jsq-lo` is
the original pathology renamed).

## What was added

A second, independent split on the **dimensionless grid-scale current**

> s = |J| · dx_min / |B|

= the fraction of the local field that reverses across one cell. `s → 1` means B flips over a
single zone: the current is at the grid scale, and is a resolution artefact as much as a physical
structure. `s ≪ 1` means a current spread over many cells — a resolved sheet. This is the right
variable precisely because it is measured **in units of the grid**, which is the thing that changes
between ladder rungs; density is blind to it.

| column | meaning |
|---|---|
| `mag-Jsq-sheet` | ∫\|J\|²dV over cells with s > `hydro/mag_diag_sheet_thresh` |
| `mag-Jsq-smooth` | ∫\|J\|²dV over cells with s ≤ threshold |
| `mag-Vsheet` | volume carrying `mag-Jsq-sheet` |
| `mag-Jsqsq` | ∫\|J\|⁴dV |

`dx_min = min(dx,dy,dz)`, so on an anisotropic cell the indicator reports the most grid-limited
direction rather than an average that could hide it. Cells with `|B| = 0` and `J ≠ 0` count as
sheet: s is then formally infinite, and a current with no field is a pure grid artefact.

## Why `mag-Jsqsq` matters on its own

Until now **`Jsq` had no `sq` companion**, unlike `dissO` and `dissA`. So `f_eff(Jsq)` could not be
formed from the history file at all — the 2026-08-06 test had to recompute it offline from 986 GB
of ladder snapshots, which is why that test could be run exactly once and only for `Jsq`. With this
column,

> f_eff(Jsq) = `mag-Jsq`² / (V_box · `mag-Jsqsq`)

is available on **every history row, for free, in every run**. The pathology that took a dedicated
offline campaign to measure is now a column.

## Gating

`hydro/mag_diag_sheet_thresh`, default **0 = OFF**, registers no columns — same discipline as the
density split, so `maxRelDivB` does not move and every existing analysis script and on-disk history
file stays valid. Registered **outside** the `have_eta` block on purpose: `Jsq` needs only B and the
grid, so an ideal run can measure it too. Suggested value 0.1–0.3 (s ≥ 0.3 ⇒ B turns over in ~3
cells).

## What is still open

Unchanged: the split has **not** been applied to `dissO`/`dissA` on the ladder, which needs either a
ladder re-run with `mag_diag_rho_split` enabled and `nonideal_eta` in the output list, or the
snapshots regenerated with that field. The current-sheet columns are likewise **implemented and
gated but not yet measured on the ladder** — this entry claims the instrument, not a result.
