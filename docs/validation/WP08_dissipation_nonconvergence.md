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
