# Item 5 — what a cross-code (AthenaPK ↔ Athena++) check can still test, and what it cannot

> **CLOSED — DROPPED BY THE USER, 2026-08-06: "forget athena++, no need for any comparison."**
> No cross-code work is planned or pending. This document is retained for one reason: §2 and §3
> are findings about the *code*, not about a proposed run, and they stay true regardless. In
> particular **§2 means no future document may claim the flagship was validated against an
> independent code** — that option does not exist in this workspace, so if such a claim is ever
> wanted it needs a different second code, not a rerun of this one.

**Assessed 2026-08-06, from the tree, not from memory.** The head-to-head was formally
**concluded 2026-07-05** by the user and its large datasets deleted. Item 5 proposes re-running it
as an independent-code check on the current production code. This document establishes what is
actually reachable before any compute is spent, because two of the three limits below are hard
and were not obvious.

## 1. What survives of the old comparison

| | state |
|---|---|
| `athena++/bin/athena` | **exists** — the STS build, 2026-06-25 |
| `athena++/runs/*` | **only `figures/`**. `ad_run2match`, `ideal_match`, `test9` — all deleted |
| AthenaPK counterparts (`run2match`, `run_ideal_match`) | deleted in the 2026-07-19 cleanup |

So there is no matched pair on disk. A cross-code check means building both sides fresh: two input
decks, two submit scripts, two runs. Nothing can be resumed.

## 2. HARD LIMIT — the flagship physics is not reproducible in Athena++

Verified against this checkout, not assumed:

| flagship ingredient | Athena++ here |
|---|---|
| Hall | `HallEMF` is a **non-functional stub** (commented out); no match in `src/field/field_diffusion/` |
| AD from chemistry-evolved x_e (`ambipolar_coeff = ionization_chem`) | **absent** — Athena++ reads a constant `eta_ad` from `<problem>`, η_A = coeff·B² only |
| ionization / Wardle tensor conductivities | **absent** |
| multigroup M1 radiation (`n_group = 3`, tabulated opacity) | in-source (`src/nr_radiation`) but **compiled OFF** |
| dust evolution, sinks | **absent** |
| tabulated protostellar EOS (H₂ dissociation) | **absent** — adiabatic EOS |

The flagship deck (`runs/root_ladder/fhc_rootladder.in`) uses `ambipolar_coeff = ionization_chem`,
tabulated EOS, `n_group = 3` + tabulated opacity, and Hall. **None of it has an Athena++
counterpart.** A cross-code check of the flagship configuration is therefore impossible in
principle here — not merely expensive. Any claim that the flagship has been "validated against an
independent code" would be false.

## 3. HARD LIMIT — the production boundary condition has no Athena++ counterpart

`BoundaryFlag` in `src/bvals/bvals_interfaces.hpp:84` is:

```
block, undef, reflect, outflow, user, periodic, polar, polar_wedge, vacuum,
shear_periodic, mg_zerograd, mg_zerofixed, mg_multipole
```

There is **no `diode`**, and `grep -rn diode src/` returns nothing. (AthenaPK's own source comment
at `bvals/boundary_conditions_apk.hpp:99` says "Athena++ calls the same thing a diode" — that
refers to Athena++ *upstream//in general*, and is **not true of this checkout**. Worth fixing the
comment.)

Since 2026-08-03 the production AthenaPK fluid BC **is** `diode`. So a matched pair must pick one:

1. **Run AthenaPK with `outflow`** for the comparison — matched, but then it is not testing the
   production configuration, and the known ~51 % spurious-inflow difference is inside the
   comparison rather than outside it.
2. **Implement a diode as an Athena++ `user` BC** — a small, well-defined piece of work (zero the
   inward face-normal momentum in the ghost zone), and the only option that makes the pair *both*
   matched and production-representative.
3. **Accept the mismatch and bound it.** Measured effect on total mass drift is
   +0.10577 % → +0.05216 % at 256³. That bound is on *mass*, not on the comparison observables,
   so this is the weakest option.

## 4. What IS testable, and what it is worth

The intersection of the two codes is: **ideal MHD + self-gravity (multigrid) + Ohmic + constant-
coefficient AD (η_A = coeff·B²), adiabatic EOS.** The historical mapping is established and
survives: code units are shared by construction, `eta_ad = Q_A` maps **1:1**, and Athena++
`mu = 31.9` ↔ AthenaPK `B0z = 0.15` ↔ 7.47 µG. The B-convention difference (Athena++ Gaussian,
AthenaPK Heaviside-Lorentz) touches only `v_A`; ratios are convention-free.

That intersection is a genuine regression check on everything that changed since 2026-07-05 —
the gravity fix, relative residual, multipole BCs, v2 tables, the whole audit batch — against an
independent implementation. It is **not** a validation of the flagship.

**Cost, from the measured history rather than a guess.** At matched grid (L=16, root 256³) the
old AD pair reached first core in **3.1 h on 2×H100**, while the 64-core CPU leg **never reached
formation in 12 h** (timed out at t = 1.040). So a first-core cross-check is GPU-cheap and
CPU-prohibitive. The affordable comparison is the **diffuse phase** (the old ideal pair used
t = 0.970 → 0.995), where both sides are hours, not days.

## 5. Recommendation

Do the **reduced-tier regression check in the diffuse phase**, with option 2 (a diode `user` BC in
Athena++) if the pair is to represent production, or option 1 with the deviation stated if not.
Do **not** describe the result as validating the flagship — write the scope of §2 into whatever
document reports it.

**Sequencing note.** The user gated items 1/5/6 on "after the D1 fix". D1 is currently *measured*,
not fixed (jobs 2468612 / 2468613 / 2468889 in flight), so this document is the assessment only;
no cross-code compute has been submitted.
