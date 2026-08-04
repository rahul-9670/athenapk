# WP-19 — the Ω·B polarity axis

**Status: PASS.** Flipping the sign of Ω relative to B changes the flux-retention observable by
**0.71 %**, against WP-18's seed-to-seed σ ≈ 16 % — **23× below the scatter**. The
single-polarity campaign is representative, and flux-retention numbers do **not** need a polarity
caveat.

Job 2450424, `build_gpu_v3` `6b1fe753`, deck `fhc_rootladder.in` `e4a7e17c`, 128³, two legs of
90 cycles each, sequential **on the same GPU** so no node-to-node confound. Read at **t = 0.90**
(WP-7: the t = 1.0 endpoint is a singular stall).

## Why this axis is not a formality

Ideal and ambipolar MHD are invariant under Ω → −Ω — the solution maps onto its mirror image, so
the polarities *must* agree. **The Hall term breaks that symmetry**: its EMF η_H (J × B)/|B| is
*linear* in B, so it changes sign with B but not with Ω. This is the origin of the aligned /
anti-aligned bimodality reported in the disc-formation literature. Production runs Hall **on**
(`hall = hall`, `hall_coeff = ionization`), so the axis is live, and nothing in this campaign had
measured its size.

`collapse_be.cpp:206-208` applies `omega_code = omegatff / tff_code` as a solid-body rotation
about z while `B0z = 0.15` lies along +z, so the sign of `omegatff` **is** the sign of Ω·B.

## Result at t = 0.90

| column | aligned (Ω·t_ff = +0.02) | anti (−0.02) | rel diff |
|---|---:|---:|---:|
| **MEtor/MEpol** | 7.471938e-03 | 7.524872e-03 | **+0.7084 %** |
| mag-MEtor | 3.713086e-01 | 3.738832e-01 | +0.6934 % |
| mag-MEpol | 4.969374e+01 | 4.968633e+01 | −0.0149 % |
| ME | 5.006499e+01 | 5.006024e+01 | −0.0095 % |
| KE | 2.425459e+03 | 2.427224e+03 | +0.0728 % |
| mag-Jsq | 2.704604e+01 | 2.704441e+01 | −0.0060 % |
| mag-dissA | 2.706161e-01 | 2.681004e-01 | −0.9296 % |
| mass | 2.053581e+03 | 2.053572e+03 | −0.0004 % |

The whole polarity effect lives in the **toroidal** field (+0.69 %) — exactly where Hall physics
should put it — while the poloidal field, total ME and mass are unchanged at the 1e-4 level. That
is the signature of a real but small Hall asymmetry, not of numerical noise.

Consistent with WP-22, which measured |η_H| < η_O in 94.8 % of cells on prod_v9.

## The helicity does not flip sign — and that is a positive control, once decomposed

`mag-Hc` is negative in **both** legs (−6.518e-02 aligned, −4.474e-02 anti), which at first
reading contradicts "the winding tracks the rotation sense". It does not, and the resolution is
quantitative rather than hand-waved.

Both legs share **one turbulent seed** (k2 sampler, seed 42). The seed's handedness is a property
of the IC and is therefore *identical* in the two legs; only the rotation-induced winding
reverses. So

```
Hc(aligned) = H_seed + H_rot        Hc(anti) = H_seed - H_rot
```

Solving:

| part | value | behaviour |
|---|---:|---|
| symmetric, seed-imposed | **−5.496e-02** | identical by construction |
| antisymmetric, rotation-induced | **∓1.022e-02** | **flips with Ω, as required** |

`|H_seed / H_rot| = 5.38` — the seed's helicity dominates the rotation-induced part by 5.4×, so
the total keeps its sign while the antisymmetric part flips exactly as the symmetry argument
demands. **The polarity flip demonstrably took effect**, which is the positive control this
campaign has learned to demand (the hst files also differ, so this is not a null test).

Checking only symmetric scalars (ME, KE) would have missed this: they agree to 0.01 % whether or
not the field geometry is mirrored.

## Incidental finding

`cap-VH`, `cap-MH`, `cap-DH` are **exactly zero in both legs** — the Hall cap never engaged
anywhere in either run at these conditions. The Hall term is active but never hit its limiter.

## Scope

128³ uniform, pre-collapse phase to t = 0.90, one seed, |Ω·t_ff| = 0.02. It does **not** license a
claim about polarity at first- or second-core densities, where η_H grows relative to η_A and the
bimodality in the literature is strongest. If production is pushed past the first core, this axis
should be re-measured there.
