# WP-20 — `turb_alpha = 3.667`: which spectral convention?

**Status: RESOLVED (source-verified) 2026-07-31.**
Binary/source pair: `docs/provenance/binary_5ebddce0/` (md5 `5ebddce0…`, AthenaPK `29a7174` + patches).
Source of truth: `src/pgen/collapse_be.cpp:313–431`. No compute was run — this is a source read,
as WP-20 specifies.

## Answer

**`turb_alpha` is NOT a 3D power-spectral-density slope. As implemented it is the
1D shell-integrated slope:  E(k) ∝ k^(−α).**

## Derivation from source

Two lines fix the convention.

1. **Mode wavenumbers are drawn uniformly in |k|** (`collapse_be.cpp:335–343`):
   ```cpp
   std::uniform_real_distribution<Real> ukmag(k_min, k_max);
   ...
   const Real kmag = ukmag(mrng);
   ```
   with the direction drawn separately and isotropically (`ucos` uniform in cos θ,
   `uang` uniform in φ, lines 344–349). Therefore the mode density is **flat in |k|**:

   > dN/d|k| = nmodes / (k_max − k_min) = const.

   This is the load-bearing step. A sampler that filled a 3D k-space *volume* uniformly
   would give dN/dk ∝ k², and the answer below would flip.

2. **Per-mode velocity amplitude** (`collapse_be.cpp:354`):
   ```cpp
   t_amp[m] = std::pow(kmag, -0.5 * alpha_s);          // |v_k| ∝ k^(−α/2)
   ```
   so per-mode energy ∝ |v_k|² ∝ k^(−α).

Combining, the energy per unit wavenumber is

> E(k) dk ∝ (dN/dk) · k^(−α) dk ∝ **k^(−α) dk**

because dN/dk is flat. The overall normalization
(`vrms_analytic = sqrt(0.5 · Σ t_amp²)`, lines 426–431) only rescales to the requested
`turb_mach`; it does not touch the slope. Nothing else in the block applies a k-space
volume weight.

## Consequence — the production value is probably a convention mix-up

| intent | E(k) | required `turb_alpha` here |
|---|---|---|
| Kolmogorov | k^(−5/3) | **1.667** |
| Burgers / supersonic | k^(−2) | **2.000** ← the code's *default* |
| what production runs | k^(−3.667) | 3.667 |

`α = 11/3 = 3.667` is the textbook **3D PSD** slope for Kolmogorov — E(k) = 4πk²P(k) gives
k²·k^(−11/3) = k^(−5/3). Under this sampler that conversion does not apply, so setting
3.667 does **not** produce Kolmogorov. It produces E(k) ∝ k^(−11/3), which is steeper than
Kolmogorov by a factor k^(−2), and steeper even than Burgers.

Corroborating (inferred, not proven): the pgen **default** is `turb_alpha = 2.0`
(line 319) = Burgers in the 1D convention — i.e. the author wrote α as the 1D slope. The
production deck value of 3.667 looks like it was chosen from the 3D-PSD table.

## How much does it actually matter?

Production band is narrow: `turb_kmin = 1`, `turb_kmax = 8` (under one decade), `turb_nmodes = 128`,
`turb_mach = 0.5`, `turb_zeta = 0.5`, `turb_seed = 42`. Amplitude ratio between the k=8 and k=1 modes:

| α | |v_k=8| / |v_k=1| = 8^(−α/2) |
|---|---|
| 3.667 (production) | 1 / 49 |
| 2.0 (Burgers) | 1 / 8 |
| 1.667 (Kolmogorov) | 1 / 5.4 |

So at α = 3.667 essentially **all** the turbulent power sits at k = 1–2. The production
"turbulent" IC is in practice a large-scale distortion of the BE sphere, not a multi-scale
cascade. Since `amp = 0.0` (the m = 2 perturbation is off) this large-scale mode **is** the
entire source of asymmetry in the collapse — and asymmetry is what the fossil-field /
magnetic-braking story is about.

## Actions

1. **Paper text must not call this "Kolmogorov."** Quote it as E(k) ∝ k^(−3.667), or restate
   in whichever convention the paper defines — but define the convention explicitly.
   (`fossil-field-paper1` standing directive.)
2. **`BOOK/` ch. on the turbulent IC**: this closes that verification to-do item; update it.
3. **Feeds WP-18**: the seed ensemble is run at this α. If the α convention is later corrected,
   the ensemble σ must be re-measured — σ is being used as the acceptance threshold for
   WP-1/2/3/7/8/9, so it is not portable across a change in α.
4. **Not a bug in the code** — the sampler is self-consistent and the default is sensible.
   It is a deck-value/convention question. Whether to re-run production at α = 2.0 is a
   science call for the user, not a numerical fix.

## Confidence

*Verified* — the sampler analysis (flat dN/dk ⇒ E(k) ∝ k^(−α)) follows directly from
lines 335–354 and there is no competing k-weight in the block.
*Inferred* — that 3.667 was entered under the 3D-PSD convention by mistake (based on the
default being 2.0). The user may have intended k^(−11/3).
