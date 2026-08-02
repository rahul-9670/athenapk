# WP-13 — Restart bit-reproducibility on the current binary

**Status: ROOT-CAUSED AND FIXED 2026-08-02. The restart blow-up was caused by
`SelfGravity::FillPoissonRHS` reading `prim` — another package's `FillDerived` output — with no
guaranteed inter-package ordering. Fixed by reading `cons` (floored to match `prim`). The fix
also corrects a previously-unknown STALE-DENSITY bug that degraded gravity on EVERY step of
every self-gravity run, not just restarts: measured Jeans growth-rate error improves 2.870% ->
2.242% against the analytic rate. IT IS THEREFORE RESULT-CHANGING; see "Adoption" below.**

> ### Two earlier root causes in this document were WRONG and are retracted
> 1. ~~"Restart drops `collapse_be_rhocrit`, silently removing the cooling task."~~ The Param
>    loss is REAL and its fix is kept, but it did NOT cause the blow-up: with the Params
>    restored and the task added, the restart still floored identically.
> 2. ~~"The radiation matter coupling fails."~~ Refuted — with `matter_coupling=false` the
>    restart still floored, with ZERO radiation warnings. It was a downstream reporter.
>
> Also falsified along the way: non-ideal diffusion (all three terms off -> still fails),
> `ion_zeta` (both runs `zeta=1e-16`), FAS solver state (`do_FAS=false` -> still fails), and
> two of my own ghost-cell theories.

Binary `09e68f75f776d5c12dfb9374bb5a5059` (`build_cpu`, GCC 13.3, `OMP_NUM_THREADS=1`).
Deck: `runs/wp13_restart/straight/fhc.in` (32³, `glmmhd`, `problem_id = collapse_be`,
`eos = hydrogen` tabulated, radiation ON, self-gravity ON) — 12 cycles straight versus
6 + restart + 6, same binary, same deck.

## Result

| | straight | after restart |
|---|---|---|
| `mass` at cycle 7 | 5.16845e+04 | **1.40608e-01** |
| `KE` at cycle 7 | finite | **nan** |
| `dt` over cycles 6→8 | 8.56e−2 → 7.93e−2 | 8.56e−2 → 1.712e−1 → 3.425e−1 (**doubling**) |
| radiation coupling failures | **0** for the whole run | **262144 cells — every cell**, first step |

The final field is density `1.0e-6` in **every cell** — exactly `dfloor`. Arithmetic check:
262144 cells × 1e−6 × dV 0.5364 = 0.1406, which is the reported `mass` exactly.

## The load is bit-exact — the failure is in the first step

Isolated before diagnosing anything else. Restarting with **zero** extra cycles
(`-r <rst> parthenon/time/nlim=6`) and dumping a snapshot reproduces the restart file
**bit-exactly**:

```
after-restart prim density : min/max = 2.554041e-01 / 8.928463e+00  sum = 9.623700e+04
restart file cons density  : min/max = 2.554041e-01 / 8.928463e+00  sum = 9.623700e+04
EXACT match: True     max |diff| = 0.0
rad.Er restored exactly: True     max |diff| = 0.0
```

The two cycle-6 restart dumps (straight vs split) are themselves `np.array_equal` identical.
So **every `Independent` variable round-trips perfectly.** Nothing is wrong with the I/O.

## Root cause — an inter-package `FillDerived` ordering dependency

`SelfGravity::FillPoissonRHS` is registered as **self-gravity's** `FillDerivedMesh`
(`self_gravity.cpp:213`) and assembled the Poisson RHS from `prim`:

```cpp
rhs = four_pi_G * (prim(IDN) - grav_mean_rho);   // over IndexDomain::entire
```

But `prim` is produced by **Hydro's** `FillDerivedMesh` (`ConsToPrim`), and Parthenon does not
guarantee ordering between two packages' `FillDerived`. Measured directly, by instrumenting the
density at the point of use:

```
restart init : prim rho entire[0, 0]  interior[0, 0]        <- prim not yet filled
fresh  init  : [0,0] -> [0,0] -> [0.357928, 3.67197]        <- good by the LAST of 3 passes
```

A fresh start survives because initialization runs several `FillDerived` passes and the last one
sees a populated `prim`. **A restart runs exactly one**, and self-gravity's ran first — so the
RHS was assembled from an all-zero density.

Zero RHS ⇒ BiCGSTAB divides by a zero residual norm ⇒ **phi = NaN**, confirmed at the consumer:

```
## WP13DBG kick sees phi[1.80e+308, -1.80e+308] nNaN=5.12e+05   <- all 512000 cells
```

NaN phi ⇒ NaN gravity kick ⇒ NaN momentum and energy ⇒ every cell driven to `dfloor`.

### The bigger finding: gravity was running on a STALE density on every step

The same ordering means that in *normal* operation `FillPoissonRHS` read the **previous** pass's
`prim` while `cons` was already current. The gravity source has therefore been solved from a
one-stage-lagged density on **every step of every self-gravity run**, not only restarts. This
silently defeated audit fix #2, whose own comment in `hydro_driver.cpp` states the intent:
*"self-gravity is applied as a STAGE-CONSISTENT source. Poisson is solved from each stage's
updated density."* It was not.

## The fix

`FillPoissonRHS` now reads density from `cons`, floored to reproduce `prim` exactly:

```cpp
const Real rho_floor = pm->packages.Get("Hydro")->Param<Real>("grav_rho_floor");
const Real rho = (cons(IDN) > rho_floor) ? cons(IDN) : rho_floor;
```

`cons` is `Independent`, restored directly from the restart file, and current at every stage, so
the ordering dependency disappears entirely. The floor is required because `ConsToPrim` applies
one in place (`adiabatic_glmmhd.hpp:152`); without it the RHS would use unfloored densities.
`grav_rho_floor` is registered unconditionally in `hydro.cpp` (an earlier attempt reached for
`cons_diag_dfloor`, which exists only when that diagnostic gate is on).

## Validation

**1. Restart vs straight — PASS.** 12 cycles straight vs 6 + restart + 6, binary `f97d51a5`:
mass and KE **identical at every cycle 6→12**.

**2. Fresh runs are NOT bit-identical — by design.** The fix corrects the staleness, so it must
change results. Measured over 12 cycles: KE 3.251e-02, time 5.744e-03, tot-E 1.339e-03,
mass 1.599e-04 (max relative). Note the density floor is *not* the cause of this difference —
adding it changed nothing on this deck; the difference is purely stale-vs-current density.

**3. Physics improves — Jeans instability against the analytic growth rate** (0.777956;
`inputs/jeans_unstable.in`, 128³, `four_pi_G=1`, `cs=0.1`, `nwave=1`). The fixed binary is closer
in **every** linear-regime fit window:

| window | old error | new error |
|---|---|---|
| amp < 1.0e−03 | 6.380 % | **5.911 %** |
| amp < 3.0e−03 | 2.870 % | **2.310 %** |
| amp < 1.0e−02 | 2.870 % | **2.242 %** |

A 22 % relative error reduction in the best-resolved window. This is the evidence that the new
behaviour is *correct*, not merely different.

## Adoption — this is RESULT-CHANGING

All existing self-gravity results were produced with the stale-density RHS. They are
self-consistent with each other but not with this binary. Anything to be compared against new
runs (njeans ladder, root ladder, WP-18 σ) must be re-baselined.

## Falsification attempts

1. **"The restart loaded garbage."** Refuted — the load is bit-exact (above). The appearance
   came from a harness bug of mine: `ls -1 *.rhdf | tail -1` selected
   `parthenon.out2.final.rhdf`, which leg 2 later **overwrote** with its own cycle-10 output.
   The file actually read at restart time was leg 1's cycle-6 dump, which is intact.
2. **"Density went negative after restart."** My error, not the code's: I indexed `prim`
   component 1 (velocity_1) as density. In the `prim` dataset index 0 is density —
   `Info/ComponentNames` concatenates *all* datasets with `grav.phi` first, exactly the trap
   CLAUDE.md documents. Re-read at index 0: bit-exact.
3. **"The three `radiation/*` params differ between logs."** Red herring — they appear inside
   a PARTHENON "set but **unused**" warning (superseded by `<units>`), in the straight run only
   because that warning is emitted at problem init.
4. **"`units.json` is missing in the restart dir."** Refuted — `runs/wp13_restart/split/` **had**
   `units.json` and failed identically (0 normalization banners, same first-step blow-up).
5. **Radiation-off control** (`runs/wp13_restart/falsify/`) — **could not be run**: the binary
   correctly aborts with `eos=hydrogen requires <physics> radiation=true (else the barotropic
   cooling overwrites the table EOS thermal energy)`. Both legs exit 134. This test needs a
   `eos = adiabatic` deck instead; it is **not yet done**, so "the momentum BC is the proximate
   cause" is *inferred from source*, not directly measured.

> ## ⚠️ SECTIONS BELOW THIS LINE WERE SUPERSEDED — REWRITTEN 2026-08-02
>
> The material that used to sit here (a "Confidence" and a "Remaining work" list) belonged to
> the **first** root cause — "restart drops `collapse_be_rhocrit`, so the barotropic-cooling task
> is silently skipped". That diagnosis was **retracted**: with the Params restored and the task
> present, the restart *still* floored. Its item 4 said *"The fix (not applied — it changes
> production physics, and is the user's call)"*, which directly contradicts the headline of this
> document. **The fix IS applied, validated, and compiled into the production binary.** The stale
> text is removed rather than annotated, because a reader who stopped at that line would conclude
> the opposite of the truth.

## Confidence (current, against the ACTUAL root cause)

*Verified*: the bit-exact restart load (all Independent vars, max|diff| = 0.0); the identical
cycle-6 dumps; that `ProblemGenerator` does not run on restart (banner 2× fresh vs 0× on restart);
that `FillPoissonRHS` read `prim` while self-gravity's `FillDerived` runs **before** Hydro's, with
the density printed at the point of use (`prim rho entire[0,0] interior[0,0]` on restart init);
that this yields `phi = NaN` in all 512000 cells at the consumer; and that reading `cons` (floored
via `grav_rho_floor`) removes the ordering dependency.

*Verified by acceptance test*: restart == straight, mass and KE identical every cycle 6→12.

*Verified against an independent standard*: Jeans-instability growth rate vs the analytic
0.777956 — the fix is closer in **every** window (amp<1e-3 6.380→5.911 %, amp<3e-3 2.870→2.310 %,
amp<1e-2 2.870→**2.242 %**). This is why the change is *correct* and not merely *different*.

*Measured on the production operator (2026-08-02, re-baseline)*: at 256³ uniform, matched state,
the fix moves `MEtor/MEpol` by **+0.032 %** and mass/MEpol by <0.001 % — ~500× below WP-18's
σ ≈ 16 %. At 128³ the shift is −0.259 %, i.e. it *converges away* with resolution.

*Not claimed*: that the two binaries compared in that re-baseline differ **only** by this fix
(`49d9c257` also carries the diagnostic packages, gated default-OFF).

## Status of the follow-ups

1. **Reproduce on the production deck** — ✅ done. GPU prod-resolution test (jobs 2446344 build /
   2446345 run) and the full re-baseline above.
2. **Radiation-off control on an `eos = adiabatic` deck** — ☐ still not done. The claim that the
   lost outside-sphere momentum BC is the proximate blow-up path for the *first* (retracted)
   root cause remains *inferred from source*. Low value now that the real cause is fixed and
   independently validated, but it is not closed.
3. **Audit every other `ProblemGenerator`-registered Param** — ✅ **done 2026-08-02, CLEAN.**
   Every `AddParam` in every pgen is reached from a restart-safe hook: `collapse_be.cpp` via
   `RegisterCollapseBeParams` ← `ProblemInitPackageData` (idempotent, key-guarded);
   `cluster.cpp` via `ProblemInitPackageData`; `shattering.cpp` via `InitUserMeshData`. **No
   remaining instances of the WP-13 pattern.**
4. **The fix** — ✅ **APPLIED.** Param registration moved out of `ProblemGenerator` into
   `ProblemInitPackageData` (runs on fresh start *and* restart, idempotent), and
   `FillPoissonRHS` now reads `cons` rather than `prim`. The interim loud-failure guard is also
   in place: `hydro_driver.cpp:68` `PARTHENON_REQUIRE`s the key when `problem_id = collapse_be`,
   so the task can never again be dropped silently. User authorized adoption ("I want good
   physics") knowing it is result-changing.

## Reproducing

```bash
source ~/athenapk_env.sh
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false PMIX_MCA_gds=hash
runs/wp13_restart/go.sh          # straight vs 6+restart+6
```

Run dirs: `runs/wp13_restart/{straight,split}`, probes in `runs/wp13_restart/probe/{load,step1}`
(`load` = restart with zero extra cycles, proving the load is bit-exact; `step1` = exactly one
step after restart, isolating the failure), failed control in `runs/wp13_restart/falsify/`.
