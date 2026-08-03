# WP-2 — what `creduc` actually scales, the RSLA validity estimate, and the sweep

**Status: PASS 2026-08-03.** Semantics source-verified (part 1, below); the three-leg sweep is
**run and analysed** (part 3, at the end). Production `creduc = 1000` costs **+0.003 %** on
`MEtor/MEpol` relative to the most faithful leg, read at t = 0.90 — four decades below WP-18's
σ ≈ 16 %. (Read instead at the t = 1.0 endpoint the figure is +0.49 %, but that endpoint sits on
the collapse singularity and is not the right comparison state — see the end of part 3.)

---

# Part 1 — what `creduc` actually scales, and the RSLA validity estimate

**Semantics RESOLVED (source-verified) 2026-08-02.**
No simulation, no GPU. Source of truth: `src/radiation/radiation.cpp:92-102`,
`src/units/physical_units.hpp:80`.

`VALIDATION_PLAN.md` WP-2 requires this question be settled *before* the sweep is interpreted,
because the direction of convergence depends on it. It is settled here.

## Answer: `creduc` is a DIVISOR, not the reduced light speed

```cpp
// src/radiation/radiation.cpp:95-102
const Real c_code  = pin->GetOrAddReal(bn, "c_code", c_default);   // c_default = U.c_code()
const Real creduc  = pin->GetOrAddReal(bn, "creduc", 1.0);
PARTHENON_REQUIRE(creduc >= 1.0, "radiation/creduc must be >= 1 (chat = c/creduc <= c).");
pkg->AddParam("c",      c_code);
pkg->AddParam("chat",   c_code / creduc);      // <-- the transport signal speed
pkg->AddParam("creduc", creduc);
```

with `c_code = c_light / v_unit` (`physical_units.hpp:80`). So:

| quantity | value |
|---|---|
| `c_code` | **1.5779e+06** code v |

| `creduc` | `chat` (code v) | `chat` (cm/s) | `chat` (km/s) |
|---|---|---|---|
| 1 | 1.578e+06 | 2.998e+10 | 299 792 |
| **300** | 5.259e+03 | 9.993e+07 | **999** |
| **1000 (production)** | 1.578e+03 | 2.998e+07 | **300** |
| **3000** | 5.259e+02 | 9.993e+06 | **100** |

**Larger `creduc` ⇒ slower `chat` ⇒ a MORE aggressive approximation.** Convergence is therefore
toward *smaller* `creduc`. The plan's sweep 300 / 1000 / 3000 correctly brackets production, but
it must be read as: 300 is the more-accurate leg, 3000 the more-approximate one. **The acceptance
test is whether creduc = 300 differs from creduc = 1000 by less than WP-18's σ ≈ 16 %** — not
whether the three legs agree with one another symmetrically.

`chat` is what the M1 transport actually uses: it is the wave speed in `M1FaceFlux`
(`radiation_moments.cpp:42-73`, `smin/smax = chat*lambda`) and the flux-limiter denominator
(`:211`). The matter coupling carries an explicit compensating factor `c/chat = creduc`
(`:248`, `:362`) so the equilibrium state stays correct under the RSLA — the standard
Skinner & Ostriker construction. That factor is present, which is a point in the implementation's
favour and was checked because getting it wrong is the classic RSLA bug.

## Validity estimate — and where it is expected to fail

The RSLA needs `chat` to exceed the fastest signal the radiation must actually carry. Two
conditions, with order-of-magnitude numbers for this problem (`v0 = 0.19 km/s`):

1. **Streaming / advection:** `chat ≫ v_flow`. Infall reaches a few × v0, call it ~2 km/s
   (10 v0). At production `chat = 300 km/s` the margin is ~150×. **Comfortably satisfied.**
2. **Dynamic diffusion:** `chat ≳ v · τ` across the optically thick region.

| region | ρ (g/cm³) | L (cm) | κ (cm²/g) | τ | v·τ at 10 km/s |
|---|---|---|---|---|---|
| first core | 1e−13 | 7.5e13 | ~1 | 7.5 | 75 km/s |
| second core | 1e−8 | 1e11 | ~1 | 1.0e3 | 1.0e4 km/s |

At the **first core** — the FHC target — `v·τ ≈ 75 km/s` against `chat = 300 km/s`: satisfied,
margin ~4×. Not lavish, which is exactly why WP-2 needs to be run rather than argued.

At the **second core**, `v·τ ≈ 1e4 km/s` exceeds `chat = 300 km/s` by ~30×. **The RSLA at
`creduc = 1000` is expected to be invalid in the second core.** If production is ever pushed to
second-core densities, `creduc` must come down (or the RSLA be dropped) and this is the number
that says so.

> *Confidence:* the `creduc` semantics and the `c/chat` coupling factor are **verified** — read
> directly from source. The τ estimates are **inferred** order-of-magnitude figures using κ ≈ 1
> cm²/g and eyeballed core sizes; they are not measured from a snapshot. They are good enough to
> predict *where* the sweep should show a signal, not to replace it.

## Part 2 — the sweep, as designed

`runs/wp2_creduc/` — three legs, 256³ uniform, t → 1.0, binary `49d9c257`, identical to the WP-7
root-ladder configuration but for `<radiation> creduc`. Same seed in all three, so this is a
**paired** comparison and WP-18's σ does not apply *between* the legs — but any statement about
the physical system drawn from them does carry it.

Judge: flux retention, `MEtor/MEpol`, `mag-Jsq`, at matched state (never at matched time — WP-8
and WP-18 both established that independently).

**Cost warning:** the legs are not equal. `chat` sets the radiation CFL, so `creduc = 300` takes
~3.3× more radiation substeps per hydro step than production. Budget the 300 leg at several times
the 1000 leg, and expect `creduc = 3000` to be the cheapest.

---

# Part 3 — the sweep, run (jobs 2446550 / 2446551 / 2448323)

All three legs reached `t = 1.000000` at cycle 179, so the comparison is at matched state by
construction, not by interpolation onto a common time.

| leg      | `creduc` | `chat` (km/s) | `MEtor/MEpol` | wall (s) |
|----------|---------:|--------------:|--------------:|---------:|
| `cr300`  |      300 |           999 |   0.0107850   |     6010 |
| `cr1000` |     1000 |           300 |   0.0108383   |     2090 |
| `cr3000` |     3000 |           100 |   0.0108585   |      977 |

Read against `cr300`, the most faithful leg (`chat = 999 km/s`, the largest `chat` affordable):

| vs `cr300`             | Δ(`MEtor/MEpol`) | Δ`ME`   | Δ`KE`   | Δ`Jsq`  |
|------------------------|-----------------:|--------:|--------:|--------:|
| `cr1000` (production)  |        **+0.49 %** | +0.07 % | +0.64 % |  +9.2 % |
| `cr3000`               |        **+0.68 %** | +0.10 % | +0.91 % | +13.5 % |

### The t = 1.0 endpoint is a singular stall — read the comparison at t = 0.90

The WP-7 root-grid ladder (same deck lineage, same 256³ uniform grid) established that **the last
0.01 t₀ of these runs sits on the collapse singularity**: `dt` falls by one to two decades between
t = 0.99 and t = 1.0 on every leg, and a uniform grid with no AMR cannot represent the forming
core. A comparison taken at t = 1.0 is therefore dominated by *where each leg happens to stall*,
not by the parameter under test. These three legs show it directly — `dt` agrees to five figures
across all three until t = 0.99, then splits 3.59e-4 / 3.44e-4 / 4.26e-4.

Re-read at **t = 0.90**, safely before that:

| leg | `chat` (km/s) | `MEtor/MEpol` | vs `cr300` |
|---|---:|---:|---:|
| `cr300` | 999 | 0.007635 | — |
| `cr1000` (production) | 300 | 0.007635 | **+0.003 %** |
| `cr3000` | 100 | 0.007635 | **+0.004 %** |

**PASS, decisively.** Through the entire magnetically-braked envelope phase the RSLA choice is
worth **0.003 %** on the flux-retention observable — four decades below WP-18's σ ≈ 16 % — and a
full decade of `chat` (999 → 100 km/s) changes nothing. Even read at the singular endpoint the
spread is only +0.49 % / +0.68 %, still 33× below σ. The RSLA is not a limiting approximation for
the first-core result on either reading.

Two things this does *not* say:

- **`mag-Jsq` moves 9–13 %** across the sweep. That is consistent with every other study in this
  campaign finding `Jsq` sensitive to numerics, and it remains **not quotable** as a physical
  result.
- **Scope is the first core only.** Part 1's estimate stands: at second-core densities
  `v·τ ≈ 1e4 km/s` against `chat = 300 km/s`, so the RSLA at `creduc = 1000` is expected to be
  **invalid** there. This sweep does not and cannot license pushing production to second-core
  densities at the current `creduc`.

Cost, for planning: `chat = 999 km/s` costs **6.15×** the wall time of `chat = 100 km/s` and
**2.88×** production's. Production's setting buys a 2.9× speedup for a 0.49 % bias — the right
trade at first-core densities.
