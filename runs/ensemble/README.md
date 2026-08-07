# Flagship Phase 7 — IC ensemble + UQ machinery

Replaces the single fiducial BE sphere with a **controlled ensemble** so the flagship reports
flux retention as a **predictive distribution**, not one curve (audit E/F). Machinery is built +
validated; **running the ensemble is the compute step** (gated on GPU availability + user go).

## Pieces
- **`design.py`** — Latin-hypercube design over the physical IC parameters × turbulence-seed
  realizations; templates a deck per point from the tabulated-multigroup `mg_prod_tab/fhc_mgtab.in`.
- **`submit_point.sh`** — one self-chaining GPU member (timeout-safe, matched-epoch density stop),
  same rtsafe multigroup+tabulated binary as mg_prod_tab; only the sampled IC keys differ.
- **`launch_ensemble.sh`** — submits every member of a design.
- **`uq.py`** — UQ analysis: flux-retention (mu_core=M_core/Phi_core) at matched physical state
  (matched rho_max, reusing flux_retention.measure_snapshot), reports the predictive DISTRIBUTION
  + separates the IC-variance (turbulence-seed) sub-class + first-order sensitivity (Spearman of
  mu_core vs each IC parameter → which knob drives the spread).

## Sampled IC parameters (documented ranges; fiducial = mg_prod_tab/prod_v9)
| param | range | fid | knob |
|-------|-------|-----|------|
| mass [Msun] | 1–10 | 6 | core mass |
| omegatff | 0.005–0.05 (log) | 0.02 | rotation |
| B0z [HL] | 0.05–0.5 (log) | 0.15 | magnetization |
| turb_mach | 0–1 | 0.5 | turbulence |
| zeta_cr_cgs | 1e-17–1e-16 (log) | 1e-16 | CR ionization (→ non-ideal) |
| dust_kappa0_cgs | 1e-4–4e-4 (log) | 2e-4 | metallicity proxy (dust opacity) |
× turbulence-seed realizations {42,101,777} → the **IC-variance** class.
(CR rate is mirrored to BOTH `<diffusion> ion_zeta` and `<chemistry> zeta_cr_cgs`, audit #1.)

## The uncertainty-class map (audit: report classes SEPARATELY)
| class | how it's covered |
|-------|------------------|
| **IC variance** | this ensemble — turbulence seeds at fixed LHS point (uq.py isolates it) |
| **IC parameter spread** | this ensemble — the LHS sweep (uq.py sensitivity ranking) |
| **discretization** | the **convergence ladder** (`../convergence_ladder/`, njeans 4/8/16) |
| **AMR / boundary** | box-size + refinement variants (add as targeted runs) |
| **RSLA** | creduc sweep (targeted runs; the RT signal-speed gate bounds it) |
| **microphysics model** | opacity_model {tabulated,belllin}, EOS table, non-ideal coeff sweeps |
| **analysis** | measurement-method variants (mu_core radius/plane definitions) |

## Workflow
```bash
# 1. generate a design (N LHS points x seeds)
python design.py --out design01 --n_lhs 8            # -> design01/run_matrix.json + point*/fhc_ens.in
# 2. launch (compute step; one self-chaining job per member, stops at first-core epoch 1e-13)
./launch_ensemble.sh design01
# 3. analyze -> predictive distribution + sensitivity
python uq.py design01 --target-rho 1.829e5           # rho_max target = 1e-13 cgs (first core)
```

## Validation
- `design.py`: 24-point design generated, all sampled keys correctly templated into each deck
  (verified value-by-value), LHS coverage 76–91% of each range.
- `uq.py`: logic unit-tested on a synthetic ensemble with a KNOWN driver — Spearman exact
  (±1.000), predictive distribution + percentiles correct, IC-variance seed sub-class isolated,
  and the injected dominant parameter (B0z) correctly ranked #1. (Note: stable sensitivity needs
  enough LHS points — small designs show spurious rank correlations.)

## Known gaps (would need pgen work, documented not hidden)
- **Alignment** (field–rotation angle) and **external pressure** are not sampled — `collapse_be`
  fixes B along z. Add a field-tilt / P_ext knob to the pgen to include them.
- Metallicity is a dust-opacity proxy (`dust_kappa0_cgs`), not a full composition change.

## RESOLVED 2026-08-06 — the 24 decks carried a SUPERSEDED boundary condition; switched to diode

**Decision: switch. Applied to all 24 decks; pilot released.** The history below is kept because
the reasoning is what justifies the switch, and because it names the trap that produced it.

Post-switch verification (all four checks passed before any GPU job was released):
- `grep` audit: **24/24 decks now on `diode`, 0 on `outflow`**, six faces each.
- Re-`diff` against the template: the ONLY differences are the BC, the explanatory comment, and
  the 6 sampled IC keys. No collateral edit.
- **Runtime init check** (`build_cpu`, point000, 32³/numlevel=2, `nlim=0`): `exit=0`,
  `Driver completed`, 64 MeshBlocks. `diode` resolved with no unknown-BC abort — the source-level
  registration check was confirmed against an actual run, not just `git show`.
- Same run positively confirms two audit fixes on the ensemble path: the **v2 EOS table** loads
  (`eos_table_hires_v2.bin`, grid 400×1000 / 400×920), and the **v2 opacity table** loads with
  **zero `legacy` warnings** — i.e. N13's self-describing-cgs path is live here, not the v1
  verbatim-code-units path. Multigroup `n_group=3` fields (`rad.Er_g1`, `rad.*_g2`) allocated.

Two notes for whoever reads the init log next:
- The `### WARNING Radiation: <radiation> rho_unit_cgs / length_unit_cgs ... IGNORED (dead key)`
  lines are **expected and harmless**. Each ensemble point has its own BE normalization (point000
  is `mass=3.658`, not the fiducial 6.0), so the deck's hardcoded `<radiation>` unit values no
  longer match what `PhysicalUnits` derives — reported as −62.8 % / +64.2 %. They are *dead keys*:
  ignored, loudly, by design. Results are unaffected. They were deliberately NOT removed, because
  deleting them buys no physics and widens the divergence from the template.
- `### Warning in Mesh::Initialize / number of MeshBlocks increased more than twice` is an
  artifact of the deliberately tiny check mesh, not a property of the production deck.

Front-end gotcha, recorded so it is not re-derived: running the binary bare on the login node dies
with `Bus error` in `MPI_Init -> psm2_ep_open`. `OMPI_MCA_mtl=^psm2` does **not** fix it — the MTL
selects `ofi`, which loads psm2 underneath. `OMPI_MCA_pml=ob1` does (it skips MTL selection
entirely). `/tmp` on the front-end is a 50 MB tmpfs and was 100 % full; set `TMPDIR`.

### The original finding (kept for the record)

All 24 decks set `ix{1,2,3}_bc = outflow`, inherited from the `mg_prod_tab/fhc_mgtab.in` template
they were generated from on 2026-07-26. The **2026-08-03 production switch**
(`docs/validation/PRODUCTION_SWITCH_2026-08-03.md`, item B1) made **`diode`** the production fluid
BC — and it was applied to **exactly one deck in the whole tree**: `root_ladder/fhc_rootladder.in`,
the flagship. Measured across `runs/`: 1 deck on `diode`, 170 on `outflow`.

Why it is not cosmetic. Parthenon's `outflow` is a plain zero-gradient ghost copy with no inflow
suppression, so mass flows *in* through nominally-outflow faces. The diode zeroes the ghost
face-normal momentum. Measured total-mass drift over a full run: **+0.10577 % → +0.05216 % at
256³**, a **~51 % cut that is resolution-independent** — the signature of a genuine boundary
effect, and independently reproduced by `cons-Mout-solver` on a 32³ smoke deck (53 %).

Why it matters *for this ensemble specifically*: the ensemble exists to put a predictive
**distribution** around the flagship's flux-retention curve. If the members run on `outflow` and
the flagship runs on `diode`, the distribution and the curve are not like-for-like — and with
`MAX_CHAIN=30` across 24 members this would only surface after up to **720 chained jobs**.

**Honest scope, do not overstate this.** The 0.05 pp figure is *total mass drift*, not flux
retention. The effect of the BC on `mu_core` has **not** been measured. The argument for switching
is "match production", not "the bias is known and large". Correspondingly, the argument for *not*
switching is that the BC is common-mode across all 24 members and so largely cancels in the
relative *spread*, which is the ensemble's primary deliverable.

Verified before proposing the change:
- `diode` is registered for all six faces in `src/main.cpp:171-189` via `DiodeBC`, and is present
  **at the exact commit the pinned binary was built from** (`4f9adff`, checked with `git show`,
  not by reading the working tree).
- The 24 decks differ from the template *only* in the 6 sampled IC keys (`diff` clean otherwise),
  so switching the BC is a one-line change per deck with no other coupling.
- `mg_prod_tab`'s own deck must be left on `outflow` regardless: that run **has data** (17 outputs,
  2026-07-26, `STOP_CHAIN` set) and its deck has to keep describing what actually produced it.

Note `src/main.cpp:174-175` states the intent explicitly: switching production to the diode "is
result-changing and therefore an explicit, deliberate deck edit." Hence: decision required, not a
silent fix.

## OPEN DECISION 2026-08-06/07 — the matched epoch (1e-13) IS rho_crit, where r_core is degenerate

**Nothing has been changed; this is evidence for a decision, not a fix.**

`STOP_CGS` and `uq.py --target-rho` both sit at **1e-13 g/cm³**, and `RHOCRIT_CODE` in
`flux_retention.py` is **1e-13/5.467e-19 = 1.829e5** — *the same number*. But `r_core` is defined
as the radius where the density profile drops below rho_crit. So at the instant rho_max = rho_crit
the super-critical region is a single point and **r_core → 0 by construction**. The campaign
measures its headline observable at a definitional singularity.

This is not hypothetical. Measured on this campaign's own snapshots:

| measurement target | degenerate cores | matched | match quality |
|---|---|---|---|
| 1e-13 (= rho_crit) | **4 of 10** | 6 | 0.06–0.17 decades |
| 1e-12 (10× rho_crit) | **0** | 4 | **0.015–0.106 decades** |

and the degeneracy tracks how close the member landed to rho_crit:

| member | matched rho / rho_crit | outcome |
|---|---|---|
| point000 | 1.43× | r_core = 0.00696, fine |
| point003 | 1.10× | **degenerate** (r_core = 2.044e-4 = innermost profile bin) |
| point002/005/006 | ~1.0× | **degenerate** |

At the 1e-12 target every remaining exclusion is merely "didn't get there" (0.42–0.90 decades
short) because the member was stopped at 1.3–2× rho_crit — a stopping-point problem, not a data
problem. The snapshots are fine.

**The exclusions are therefore NOT random.** They correlate with proximity to rho_crit, so the
surviving distribution is drawn preferentially from members that happened to overshoot. That is a
biased sample, not a predictive distribution.

**Cost of fixing it, measured on point000, which actually made the trip:** reaching rho_crit took
80.4 min; going on to 10× rho_crit took a further **92.9 min = +116 % wall time per member**. It is
expensive precisely because it is the post-first-core regime where dt collapses.

**Recommendation (user's call — it roughly doubles the campaign):** raise `STOP_CGS` to ~2e-12 and
measure `mu_core` at 1e-12.

**NOTE THE SELF-INFLICTED PART.** The epoch-stop watchdog (`watchdog_epoch_stop.py`) tightened
overshoot from 256× down to 1.1–2×, which is excellent for compute and lands members squarely in
the degenerate regime. The members that measure *cleanly* today are the ones that ran *before* the
watchdog. Compute efficiency and measurability pull in opposite directions here; if `STOP_CGS` is
raised, the watchdog needs raising with it or it will keep trimming exactly the margin the
measurement requires.

## Status
Machinery COMPLETE + validated. Launch is **held** on the boundary-condition decision above.
Otherwise: running the ensemble (+ the convergence ladder for the discretization class) is the
remaining COMPUTE, gated on GPU availability + user go. This is the LAST code-development piece of
the flagship program (per the user's plan: ensemble is last, then resume prod_v9 once the flagship
is fully done).

Also done 2026-08-06: all 24 decks repointed to the v2 EOS/opacity tables (audit N13/N14), and
`submit_point.sh` pinned to `athenaPK_PRESERVED_84a6d248` instead of the mutable `build_gpu/bin/`
scratch slot — a rebuild mid-campaign would otherwise have silently changed the binary underneath
half the ensemble.
