# Production configuration change, 2026-08-03 — binary + two physics switches

Three things changed together, on the user's instruction to "adopt both switches now, re-baseline
once". They are recorded together because they were validated together and share one re-baseline.

| | before | after |
|---|---|---|
| GPU binary | `49d9c257` (`build_gpu`) | **`f181c0a1`** (`build_gpu_v2`) |
| fluid BCs | `outflow` | **`diode`** (B1) |
| gravity convergence test | absolute `1e-6` | **relative `1e-6`** (B3) |

Deck: `runs/root_ladder/fhc_rootladder.in`, md5 **`c67eb458`**.
Re-baseline legs: `runs/root_ladder/{r128_sw, r256_sw, r512_sw}` via `stage_switches.sh`.
The Gate-A-validated `r*_gfix` legs are **kept** — they are the `49d9c257` + outflow + absolute
reference the new runs are read against.

---

## 1. The binary alone changes nothing — gated byte-identical

`f181c0a1` = `49d9c257` + the diode BC + the B2 solver-converged flag and its now-unconditional
warning + the B3 startup notice + the B4 coarse-EOS-table warning. Every one of those is opt-in or
diagnostic-only, so with the deck unchanged the two binaries must agree exactly.

Job 2448500 ran both on the **production deck** (multipole self-gravity + ambipolar + RKL2 + AMR,
128³, `nlim = 40`, 2 GPUs) and byte-compared the history:

```
PASS: history files are BYTE-IDENTICAL -> f181c0a1 is an OFF-state no-op, safe to swap.
```

That is what licenses treating everything below as an effect of the two *switches*, not of the
binary.

## 2. Verifying the switches are actually live

Checked positively rather than assumed, because a deck edit that silently fails to take is exactly
the failure mode that produced a superseded IC earlier in this campaign.

- **Deck md5** echoed by the job itself is `c67eb458`, matching the edited file, and
  `submit_root.sh` passes no BC or residual overrides that could shadow it.
- **B3 notice, same-binary A/B.** The `## NOTE [self-gravity] convergence uses an ABSOLUTE residual
  tolerance` banner is printed *only* when `relative_residual` is false. On `f181c0a1`:

  | run | deck | notice printed |
  |---|---|---|
  | `gate_v2/new` | old (absolute) | **1** |
  | `r128_sw` | new (relative) | **0** |

  Same binary, opposite decks, opposite result ⇒ the switch is live. (Comparing against
  `r128_gfix` would *not* have worked: it ran on `49d9c257`, which predates the notice, so its 0
  is uninformative.)
- **Diode**, checked by its physical signature — see §3.

## 3. Result — matched state, t = 0.90, at two resolutions

Read at t = 0.90, not the t = 1.0 endpoint (WP-7: the last 0.01 t₀ is a collapse singularity on
these uniform grids).

| quantity | `r128_gfix` (old) | `r128_sw` (new) | change |
|---|---|---|---|
| `MEtor/MEpol` | 7.459354e-03 | 7.471938e-03 | **+0.169 %** |
| `ME` | 5.014867e+01 | 5.006499e+01 | −0.167 % |
| `KE` | 2.425442e+03 | 2.425459e+03 | +0.001 % |
| `mass` | 2.055642e+03 | 2.053581e+03 | −0.100 % |
| `Jsq` | 2.704267e+01 | 2.704604e+01 | +0.012 % |
| `maxRelDivB` | 5.341537e-02 | 5.343959e-02 | +0.045 % |

At 256³ (`r256_sw` vs `r256_gfix`, job 2450079):

| quantity | `r256_gfix` (old) | `r256_sw` (new) | change |
|---|---|---|---|
| `MEtor/MEpol` | 7.634908e-03 | 7.641359e-03 | **+0.085 %** |
| `ME` | 5.013619e+01 | 5.009452e+01 | −0.083 % |
| `KE` | 2.436194e+03 | 2.436213e+03 | +0.001 % |
| `mass` | 2.053600e+03 | 2.052570e+03 | −0.050 % |
| `Jsq` | 3.021358e+01 | 3.021522e+01 | +0.005 % |

**The combined switch is +0.169 % at 128³ and +0.085 % at 256³ on the flux-retention observable —
~95× and ~190× below σ ≈ 16 %, and it HALVES with resolution.** Nothing concluded on the old
configuration is invalidated.

The halving matters more than the magnitude. Both switches are boundary/tolerance effects whose
influence on the interior should weaken as the grid refines, and the measured 0.169 → 0.085 is a
factor 1.99 against a factor 2 in resolution. A change that merely *happened* to be small would
not converge away this cleanly — this is the same signature that made Gate A's gravity-fix result
credible, arrived at independently.

### The diode is doing what it was added to do

Mass drift over the full run:

| resolution | outflow (`*_gfix`) | diode (`*_sw`) | cut |
|---|---|---|---|
| 128³ | +0.21203 % | +0.10431 % | **50.8 %** |
| 256³ | +0.10577 % | +0.05216 % | **50.7 %** |

**The cut is ~51 % and is resolution-independent**, which is what a boundary-condition effect must
be — the diode removes a fixed *fraction* of the spurious inflow (the ghost-momentum half),
regardless of how well the interior is resolved. It also independently reproduces the 53 % cut
measured on the 32³ smoke deck via `cons-Mout-solver`: three configurations, two different
instruments, same number. That is what makes the mechanism credible rather than merely consistent.

(The drift itself halves with resolution, 0.212 → 0.106 %, on both branches — the boundary error
converges away even though the *fraction* the diode removes does not.)

It does **not** eliminate the inflow, and it was never expected to: zeroing the ghost-cell
face-normal momentum cannot stop a gravity-accelerated interior cell from advecting ghost density
across the face. The remaining drift (+0.104 % at 128³, +0.052 % at 256³) needs a boundary **face-flux** clamp
applied after the Riemann solve — the same surgery B6 requires, and the two should be done together.

## 4. Why these two switches, in one line each

- **B1 / diode.** Parthenon's `outflow` is a plain zero-gradient ghost copy with no inflow
  suppression, so mass flowed *in* through nominally-outflow faces. `diode` is the standard fix.
  Full analysis: the B1 entries in `CAMPAIGN_STATUS_2026-08-02.md`.
- **B3 / relative residual.** The absolute criterion does not fail — but its *fractional* accuracy
  is an artifact of units (measured: 4.4e-11 at ρ = 1, 3.1e-19 at ρ = 1e10), while the relative
  criterion holds `res/rms(rhs) = 4.4257e-14` constant to five figures over ten decades of source.
  φ's gradient drives the collapse, so a fractional demand is the one that means the same thing at
  every depth. Full analysis: `B3_gravity_tolerance.md`.

## 5. What is deliberately NOT changed

- **`runs/prod_v9/fhc.in` is untouched.** That run is HELD, and its own header records that a
  restarted run takes parameters from the rhdf + CLI, not from the deck — so editing it would be
  both out of scope and ineffective.
- **B10's chemistry instrument and B11's per-cell Hall floor are not in `f181c0a1`.** They are in
  `build_gpu_v3` = `6b1fe753`, built into a NEW directory precisely so the live re-baseline chain
  (which re-execs `build_gpu_v2` on restart) was never disturbed. Swapping production to v3 is a
  separate, gated decision.
- **`chemistry/nsub_max` is still 400**, despite WP-10 showing it binds in 100 % of cells at zero
  wall cost. Raising it is a *result-changing* physics change and would need its own re-baseline;
  it was not part of the authorised single re-baseline. Recorded as a pending decision.
