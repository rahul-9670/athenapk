# WP-22 — numerical vs physical resistivity in the first core

**Status: PASS on the first core. Measured against the APPLIED diffusivity, not an estimate.**

- **The applied ambipolar diffusivity exceeds the numerical-resistivity bound by ~20–130×** at
  first-core densities. The fossil-field argument is not compromised there.
- **At second-core densities (ρ ≳ 10⁴ ρ_crit) the margin is only ~2–4×.** The code is not
  validated for the second core — the same boundary WP-2 drew independently for the RSLA.
- **The offline estimate was optimistic by two to three decades and had to be replaced.**
  Reconstructing η_A from the written x_e gives 10⁴–10⁵×, but the kernel's `eta_eq` ceiling
  **suppresses it by 1486× and 2678×** in exactly the first-core shells. §4 has the measurement.
  The conclusion survives; the headroom does not.

This is the gap flagged in `CAMPAIGN_STATUS_2026-08-02.md` §5 as "the one most likely to affect
the paper's conclusion", since a fossil-field result rests on physical resistivity dominating
numerical resistivity and that had never been checked.

---

## 1. Numerical resistivity — measured, not estimated

An Alfvén wave with **zero** physical resistivity decays only through the scheme's own
dissipation, so the WP-14 Alfvén ladder doubles as a measurement of η_num at no extra compute:

```
B_perp(t) = B_perp(0) exp(-eta_num k^2 t / 2)
=>  L1(B_perp) over one period T  ~=  A0 * eta_num k^2 T / 2
=>  eta_num  ~=  2 (L1/A0) / (k^2 T)
```

L1 also contains **dispersion** error, which does not damp the wave, so this over-attributes error
to dissipation: every number below is a rigorous **upper bound** on η_num — the conservative
direction for this WP. (`docs/validation/scripts/wp22_eta_numerical.py`; geometry follows
`src/pgen/linear_wave_mhd.cpp`, background `d0 = 1`, `bx0 = 1` ⇒ `v_A = 1` in Heaviside-Lorentz —
**no 4π**.)

| N | dx | η_num (bound) | η_num/(v_A dx) |
|---|---|---|---|
| 16 | 0.09375 | 5.54e-03 | 0.0591 |
| 32 | 0.04688 | 1.54e-03 | 0.0329 |
| 64 | 0.02344 | 3.87e-04 | 0.0165 |
| 128 | 0.01172 | 9.14e-05 | 0.0078 |

**η_num ~ dx^1.97** — the numerical resistivity vanishes at the scheme's full second order rather
than sitting at a floor that refinement cannot beat. That is the single most important structural
fact here: it means the margin below *improves* with resolution.

The dimensionless coefficient is not constant because a 2nd-order scheme's numerical diffusivity
carries a factor of the wavenumber. Writing `η_num = C k v_A dx²` and using k = 2π/λ:

| N | k·dx | η_num/(v_A dx) | C |
|---|---|---|---|
| 16 | 0.589 | 0.0591 | 0.100 |
| 128 | 0.0736 | 0.0078 | 0.106 |

`C = 0.103 ± 0.003` across a factor of 8 in dx — a genuine constant, so the formula transfers:

> **η_num ≈ 0.103 · k · v_A · dx²  =  0.647 · v_A · dx / N_cell**, for a structure resolved by
> `N_cell` cells per wavelength.

## 2. Physical η_A — first reconstructed offline (superseded by §4, kept for the comparison)

η is not a written field, but production runs `ambipolar_coeff = ionization_chem` with
`xe_scalar_index = 4`, and x_e **is** written (passive scalar 4). So η_A can be recomputed exactly
as `Ionization::AmbipolarEtaFromXe` does (`docs/validation/scripts/wp22_eta_physical.py`):

```
eta_A[cgs] = B^2 mu_n / (4 pi gamma_AD m_ion x_e rho^2)      mu_n = 2.33 (NEUTRAL, not the
eta_A[code] = eta_A[cgs] * eta_unit                          thermal 2.29 -- audit finding #4)
```

The explicit 4π is the cgs Langevin-drag factor, **not** a unit conversion — B_unit already
carries the Heaviside-Lorentz √(4π).

Snapshot: `runs/prod_v9/parthenon.out1.00096.phdf`, t = 1.09456, 2402 blocks, AMR levels 0–12,
finest dx = 1.526e-05 code = 4.29e+11 cm, ρ_max = 6.66e+10 code = **3.6e5 × ρ_crit**.

## 3. The offline comparison — SUPERSEDED by §4

> **Read §4 instead for the result.** This section reconstructs only the *chemistry branch* of
> η_A. The kernel additionally applies an `eta_eq` ceiling that §4 measures at **1486×** and
> **2678×** in the two first-core shells, so every ratio below is too high by ~3 decades. It is
> retained because the *shape* it shows is correct and because it is what the ceiling is measured
> against.

η_num evaluated per cell at that cell's own dx and v_A. `N_cell` = cells per wavelength of the
structure being resolved; N_cell = 4 is the pessimistic reading (barely-resolved features), 16 the
generous one.

| ρ range (code) | cells | median η_A | η_A/η_num (N_cell=4) | (N_cell=8) | (N_cell=16) |
|---|---:|---:|---:|---:|---:|
| 1e-2 – 1e0 | 18.0 M | 8.93e-02 | 3.5e+01 | 7.0e+01 | 1.4e+02 |
| 1e0 – 1e2 | 22.9 M | 7.75e-02 | 3.1e+02 | 6.2e+02 | 1.2e+03 |
| 1e2 – 1.8e3 | 16.5 M | 5.14e-01 | 6.4e+03 | 1.3e+04 | 2.6e+04 |
| **1.8e3 – 1.8e5 (→ρ_crit)** | 13.9 M | **4.11e+00** | **9.9e+04** | 2.0e+05 | 4.0e+05 |
| **1.8e5 – 1.8e7 (first core)** | 6.5 M | **1.53e-01** | **6.0e+04** | 1.2e+05 | 2.4e+05 |
| 1.8e7 – 1.8e9 | 839 k | 7.90e-05 | 4.3e+02 | 8.6e+02 | 1.7e+03 |
| **> 1.8e9 (second core)** | 37 k | **2.86e-08** | **2.2e+00** | 4.4e+00 | 8.8e+00 |

These ratios are **not the answer** — see §4. What survives from this table is the *shape*.

η_A is non-monotonic in ρ — rising to a peak of 4.1 just below ρ_crit, then falling by eight
decades. That is the expected AD behaviour, not a bug: η_A ∝ B²/(x_e ρ²), so it climbs as x_e drops
through the C→CO ionization minimum and collapses once thermal ionization of K and H switches on
above ~10³ K. The falsification check is that the shape follows the *ionization* curve rather than
tracking ρ or B monotonically, and it does.

**The second-core caveat is real and must be carried** — and it is the one number here that §4
does *not* move, because the ceiling suppresses nothing there (1.02×). In the deepest shell the
margin is only 2–4×, i.e. numerical and physical resistivity are comparable. No flux-retention statement should be
made at ρ > 10⁴ ρ_crit on this grid. WP-2 reached the same boundary from a completely different
direction (RSLA `v·τ ≈ 1e4 km/s` vs `chat = 300 km/s`), which is corroboration rather than
coincidence: both say *first core yes, second core no*.

## 4. The applied η_A — and why the offline estimate had to be replaced

The kernel applies

```cpp
eta_A = min(eta_chem, eta_eq, eta_ad_cap)          // diffusion.hpp:310-341
```

- `eta_ad_cap` is **inactive**: `eta_ad_cap_code` is not set in the production deck and defaults to
  `numeric_limits<Real>::max()` (`hydro.cpp:1085`).
- `eta_eq`, the equilibrium NICIL/Wardle value, is **active** as a self-calibrating ceiling, added
  because "in dense gas the chemistry x_e collapses … which makes the single-fluid eta_A blow up
  and drives the explicit dt to ~0". Dense gas is the first core, so the ceiling sits exactly where
  the claim lives.

§3's number was therefore an **upper bound**, which is the wrong direction for a domination claim,
and `eta_eq` cannot be reproduced offline without porting the whole grain + Saha charge solve.

So it was measured instead. The code already caches the applied (η_O, η_H, η_A) per cell in the
Parthenon field `nonideal_eta` (`hydro.cpp:1391-1395`) — it is simply never written. Job 2448571
restarted `prod_v9` (HELD since 2026-07-23, read-only, on its own binary `17af621a`) with that
field added to `output1`'s variable list, making the applied diffusivity a phdf field.

| ρ range (code) | cells | η **applied** | η_chem (§3) | suppression | η_A/η_num (N=4) | (N=8) |
|---|---:|---:|---:|---:|---:|---:|
| 1e-2 – 1e0 | 18.0 M | 2.12e-02 | 8.93e-02 | 4.2× | 8.4 | 16.7 |
| 1e0 – 1e2 | 22.9 M | 1.46e-03 | 7.75e-02 | 53× | **5.8** | 11.7 |
| 1e2 – 1.8e3 | 16.5 M | 1.57e-03 | 5.14e-01 | 328× | 19.6 | 39.3 |
| **1.8e3 – 1.8e5 (→ρ_crit)** | 13.9 M | 2.77e-03 | 4.11e+00 | **1486×** | **66.8** | 134 |
| **1.8e5 – 1.8e7 (first core)** | 6.5 M | 5.69e-05 | 1.52e-01 | **2678×** | **22.3** | 44.6 |
| 1.8e7 – 1.8e9 | 836 k | 2.94e-06 | 7.81e-05 | 27× | 16.0 | 31.9 |
| **> 1.8e9 (second core)** | 37 k | 2.81e-08 | 2.87e-08 | **1.02×** | **2.1** | 4.3 |

**The ceiling binds hard, and precisely where it matters** — 3.2 and 3.4 decades of suppression in
the two first-core shells, and essentially none (1.02×) in the second core where the chemistry x_e
has recovered. The doc previously reasoned that `eta_eq` "would have to suppress η_A by four to
five decades before the conclusion changed"; it suppresses by 3.4, which is why this was measured
rather than argued.

**PASS stands, with much less headroom than §3 implied.** Applied physical η_A beats the numerical
bound by 22–67× at the first core on the pessimistic N_cell = 4 reading, 45–134× on N_cell = 8 —
one to two orders of magnitude, not four to five. Since η_num is an over-estimate (it absorbs
dispersion error) the true margin is larger, and since η_num ~ dx^1.97 it improves with resolution.
But the claim is now "physical AD dominates by ~1.5 decades", not "by 5 decades", and no statement
should be made at ρ > 10⁴ ρ_crit, where the margin is ~2×.

**Ohmic cap engagement, incidentally measured:** `eta_ohm_cap_code = 0.1` is on the ceiling in
**8441 of 78.7 M cells (0.011 %)**. That numerical stabilizer is not shaping the result.

## 5. Scope

- Applies to **ambipolar** diffusion, the dominant non-ideal term for flux transport here. Ohmic
  and Hall are separately capped (`eta_ohm_cap_code = 0.1`, `hall_ohmic_floor_code = 0.05`); those
  caps are deliberate numerical stabilizers and their engagement is what `cap_diag` measures.
- η_num is measured on a **smooth** Alfvén wave. At a shock or a sharp current sheet the scheme's
  dissipation is locally larger (first order at a discontinuity, by design). The comparison is
  valid for the smooth pseudo-disc/envelope field where flux transport is decided, not inside a
  captured shock.
- The snapshot is from `prod_v9`, which predates the WP-13 gravity fix. η_A and η_num are both
  *local* functions of (ρ, B, x_e, dx, v_A), so the comparison is unaffected at the decade level;
  the flux-retention numbers from that run are a separate matter and are not used here.

> *Confidence:* η_num — **measured** (WP-14 ladder). Applied η_A — **measured**, read out of the
> code's own cache on a real snapshot, no reconstruction. The offline η_chem in §3 is retained only
> to quantify the ceiling's effect. The second-core caveat — **measured**, and independently
> corroborated by WP-2.
