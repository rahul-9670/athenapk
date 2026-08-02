# WP-21 — BEProfile analytic-fit error

**Status: RESOLVED — PASS. 2026-07-31.**
Binary/source pair: `docs/provenance/binary_5ebddce0/`.
Reproduce: `venvs/analysis_env/bin/python docs/validation/scripts/wp21_be_profile.py`
(pure analysis — no simulation, no GPU).

## Question

`src/pgen/collapse_be.cpp:66-68` initializes the sphere from a closed-form approximation
("Tomida 2011 PhD thesis"), not from an integrated Bonnor–Ebert solution:

```cpp
Real BEProfile(Real r) { return std::pow(1.0 + r*r/rcsq_code, -1.5); }   // rcsq = 26/3
```
truncated at `rc_code = 6.45` (line 40), with `bemass_code = 197.561` (line 42).

WP-21 asks how far that is from the truth, and therefore **whether t = 0 is where we think it is.**

## Method

Integrated the isothermal Lane–Emden equation directly,

> ψ'' + (2/ξ)ψ' = e^(−ψ),  ψ(0) = ψ'(0) = 0,  ρ/ρ_c = e^(−ψ)

with DOP853 at rtol 1e-12 / atol 1e-14 (series expansion near the 2/ξ singularity). Because
`four_pi_G_code = 1` and c_s = 1, the code length unit **is** the isothermal scale height
a = c_s/√(4πGρ_c), so code radius r is exactly the dimensionless ξ — no conversion needed.

**Integration validated independently:** the enclosed mass computed by quadrature and by the
Lane–Emden first integral M = 4πξ²ψ' agree to 7 significant figures (197.3190 both ways), and
the recovered critical radius is ξ_crit = 6.4508 against the textbook 6.451.

## Results

| quantity | exact (integrated) | code approximation | error |
|---|---|---|---|
| critical radius ξ_crit | 6.4508 | 6.45 (`rc_code`) | **−0.012 %** |
| ρ(edge)/ρ_c | 0.071235 | 0.071586 | +0.49 % |
| centre-to-edge contrast | 14.038 | 13.969 | −0.49 % (textbook 14.04) |
| enclosed mass (code units) | 197.3190 | 197.5614 | **+0.123 %** |
| mean density | 0.175550 | 0.175766 | +0.123 % |

**Density profile error over 0 ≤ ξ ≤ 6.45:**

- **max |rel. error| = 0.937 %**, at ξ = 1.99
- RMS rel. error = 0.653 %
- exact at the centre (both normalized to ρ_c)

| ξ | ρ_exact | ρ_approx | rel. err |
|---|---|---|---|
| 0.00 | 1.000000 | 1.000000 | +0.000 % |
| 1.00 | 0.853143 | 0.848913 | −0.496 % |
| 2.00 | 0.571310 | 0.565959 | **−0.937 %** |
| 3.00 | 0.345302 | 0.343595 | −0.495 % |
| 4.00 | 0.207581 | 0.208263 | +0.329 % |
| 5.00 | 0.129498 | 0.130610 | **+0.859 %** |
| 6.45 | 0.071235 | 0.071586 | +0.492 % |

The error is **sign-structured**, not random: the approximation is slightly *under*-dense at
intermediate radii (ξ ≈ 1–3.5) and slightly *over*-dense in the envelope (ξ ≳ 4). Net effect is
the +0.123 % mass excess.

## Does t = 0 move?

Since t_ff ∝ ρ̄^(−1/2), the +0.123 % mean-density excess shifts the free-fall clock by

> **Δt_ff / t_ff = −0.061 %**

i.e. the sphere starts collapsing 0.06 % *early* relative to a true critical BE sphere of the
same ρ_c and radius. **t = 0 is where we think it is.**

## Incidental finding — `bemass_code` is self-consistent

`bemass_code = 197.561` reproduces this analysis's mass of the **approximation**
(197.5614) to five significant figures, *not* the exact BE mass (197.3190). So the constant was
derived from the same approximate profile the pgen actually lays down. That is the correct,
self-consistent choice — the mass normalization and the density field agree with each other, and
there is no double-counting. The 0.12 % offset is against the *idealized* BE sphere only.

## Verdict

**PASS, with a wide margin.** Every error is ≤ 1 %, and the quantity that matters for the
collapse clock is 0.06 %. For context, the IC deliberately applies an overdensity factor
**f = 5** to drive collapse — a 5× departure from hydrostatic equilibrium by construction. A
0.9 % shape error is three orders of magnitude below the deliberate perturbation and cannot be
a limiting error anywhere in this campaign.

No action required. Cite these numbers in the paper's IC section rather than describing the
profile as "a Bonnor–Ebert sphere" without qualification — it is a 0.9 %-accurate closed-form
approximation to one, truncated at the critical radius.

## Confidence

*Verified.* The ODE integration is cross-checked two independent ways (quadrature vs. first
integral, and recovered ξ_crit vs. textbook), and the unit identification (code r = ξ) follows
from `four_pi_G_code = 1` with c_s = 1.
