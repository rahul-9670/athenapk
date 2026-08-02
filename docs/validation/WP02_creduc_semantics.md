# WP-2 part 1 — what `creduc` actually scales, and the RSLA validity estimate

**Status: SEMANTICS RESOLVED (source-verified) 2026-08-02. The compute sweep is STAGED, not run.**
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

## The sweep (staged, not run)

`runs/wp2_creduc/` — three legs, 256³ uniform, t → 1.0, 8 GPUs, binary `49d9c257`, identical to
the WP-7 root-ladder configuration but for `<radiation> creduc`. Same seed in all three, so this
is a **paired** comparison and WP-18's σ does not apply *between* the legs — but any statement
about the physical system drawn from them does carry it.

Judge: flux retention, `MEtor/MEpol`, `mag-Jsq`, and core temperature at matched state (never at
matched time — WP-8 and WP-18 both established that independently).

**Cost warning:** the legs are not equal. `chat` sets the radiation CFL, so `creduc = 300` takes
~3.3× more radiation substeps per hydro step than production. Budget the 300 leg at several times
the 1000 leg, and expect `creduc = 3000` to be the cheapest.
