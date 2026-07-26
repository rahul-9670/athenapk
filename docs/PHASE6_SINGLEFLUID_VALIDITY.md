# Phase 6 gate — single-fluid validity map (RESULT)

**Date:** 2026-07-26. Diagnostic: `src/hydro/diffusion/singlefluid_validity.py`.
Tested on the deepest available FHC snapshot: `prod_v9/parthenon.out1.00096.phdf` (t=1.0946,
ρ_max=3.64e-8 g/cm³ — spans the whole collapse, diffuse envelope 2e-20 → first-core/2nd-collapse
threshold 3.6e-8, 78.7M cells).

## Question
Is the single-fluid strong-coupling approximation (the Wardle-tensor AD+Hall+Ohm non-ideal MHD the
code evolves) valid throughout the collapse, or is a true multifluid treatment (separate ion
inertia / generalized Ohm law, Phase 6) required anywhere?

## Method (grounded in `ionization.hpp`)
Per cell, from (ρ, |B|, x_e): the ion-neutral collision rate `ν_in = (m_n/(m_i+m_n)) n_n <σv>`
(Langevin, `<σv>=1.9e-9`, `m_ion=24.3 m_H`, `μ_n=2.33` — the code's constants), the dynamical
rate `ω_ff = √(4πGρ)`, and the signed ion Hall parameter `β_i = eB/(m_i c ν_in)`.
**Criterion:** two-fluid → single-fluid when ion inertia is negligible, i.e. ions equilibrate
within a dynamical time: **`χ_ion = ν_in/ω_ff ≫ 1`**.

## Result — SINGLE-FLUID NON-IDEAL MHD ADEQUATE
- **`χ_ion = 6.6×10⁶ … 2×10¹²` across the entire collapse; min = 6.6×10⁶ ≫ 1. Single-fluid valid
  in 100% of cells.** The ions are collisionally locked to the neutrals on the dynamical time at
  every density (ν_in ∝ n_n is huge), so ion inertia is negligible — the textbook reason
  single-fluid AD-MHD holds in star formation. **Phase 6 (full multifluid) is NOT required for the
  FHC collapse.**
- The result is **x_e-independent** (χ_ion depends on ν_in ∝ n_n, not the ionization fraction), so
  it is robust to the chemistry x_e floor in the dense core.
- **Self-validation against canonical theory:** the Hall parameter β_i transitions
  **ambipolar (β_i~350, diffuse envelope) → Hall (β_i~0.5, ρ~1e-13, first-core onset) → resistive/
  Ohmic (β_i≪1, dense core)** — exactly the known non-ideal MHD regime sequence for protostellar
  collapse (Wardle 2007; Tsukamoto et al.). The diagnostic reproduces established physics.

## Falsification note (load-bearing)
A first pass used the NEUTRAL-ion rate `ν_ni ∝ ρ_i` and reported "multifluid needed in 80% of
cells" — which **contradicts the established validity of AD single-fluid MHD in star formation**,
so it was the diagnostic that was wrong, not the theory. The error: `ν_ni/ω_ff < 1` at high ρ only
means AD is *active* (neutrals not flux-frozen — the intended non-ideal regime), NOT that ion
inertia matters. The correct single-fluid criterion is the ION rate `ν_in ≫ ω_ff`. Fixed; the
script now guards against and documents this pitfall.

## Bottom line for the flagship
Phase 6's gate ("single-fluid validity map") is delivered: the map shows the single-fluid
non-ideal treatment is adequate everywhere the FHC collapse reaches (through the 2nd-collapse
threshold). The expensive multifluid build is **not on the critical path** for the FHC flux result.
(Caveat: assessed on the fiducial B0z=0.15 collapse; a strongly-ionized or much-lower-B regime
should be re-checked, but those are not the fiducial FHC problem.)
