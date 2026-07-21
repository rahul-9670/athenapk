Subject: Core rotation count + quantitative state of the 6 Msun collapse run (prod_t4_full)

Dear Professor [Name],

You asked how many rotations the core has completed. Short answer: **through the end of its quasi-static (first-core plateau) phase the core had completed ~0.1 rotations; integrating through our current simulation front gives ~0.24 rotations** — with the caveat, explained below, that the final runaway segment is being re-simulated after we caught a solver bug, so end-state numbers are preliminary. The core is spinning up rapidly (Omega ~ rho^0.57; rotation period ~7 yr at the current front).

**How N_rot is calculated.** At 16 snapshot epochs from first-core formation onward I compute the mass-weighted angular velocity of the core as a body, Omega = |L|/I, about the instantaneous density peak (bulk motion subtracted, I about the measured angular-momentum axis), and integrate N = (1/2pi) Int Omega dt — the true accumulated rotation angle. The commonly quoted alternative (core age divided by the *current* rotation period) assumes the core always spun at today's rate and overestimates by ~2 orders of magnitude during spin-up. Inner material naturally leads this body-average.

**Comparison with the literature.** RHD collapse calculations (Masunaga & Inutsuka 2000; Tomida et al. 2013; Vaytet & Haugbolle 2017) find first cores complete of order one to a few rotations before second collapse, the count scaling with the core lifetime (~10^2-10^3 yr, inversely proportional to accretion rate). Ours sits at the low end for two reasons: (i) the initial cloud is a slow rotator (beta_rot = 3.6e-4, low end of the Goodman et al. 1993 range), and (ii) the f = 5 overdensity drives accretion ~100x the Shu rate, so the first core lives only ~10^2 yr — less time to turn.

**State of the run** (front at cycle 77,250; "clean" = from the validated portion of the run, "prelim" = from the final runaway segment that is being re-simulated):

| Quantity | Value | Status |
|---|---|---|
| Simulated time | 50.80 kyr = 1.26 initial free-fall times | clean |
| First core formed | at ~50.4-50.6 kyr | clean |
| Quasi-static plateau | ~23 yr at rho ~ 1e-10 g/cm3 | clean |
| Peak density reached | 1.44e-7 g/cm3 (n ~ 4e16 cm-3) | prelim |
| Central temperature | expected ~1500-2200 K at the front (near-adiabatic, exact-EOS isentrope; dissociation onset imminent or under way) | prelim (see caveat) |
| Core mass | 0.0029 Msun (central object, r_rms 0.28 AU); 0.017 Msun (first core, r_rms 9 AU) | prelim |
| Mass accretion rate | 5.2e-4 Msun/yr at 1 AU; 1.2e-4 at 10 AU (~100x Shu c_s^3/G) | prelim |
| Infall speed | 2.5 km/s = 13 c_s at 0.35 AU | prelim |
| Magnetic field | 7.47 uG initial cloud -> ~8 mG core plateau | prelim (see fossil-field note) |
| Magnetic flux (equatorial) | 5.5e24 G cm2 within 1 AU; 4.6e26 within 10 AU | prelim |
| Flux accretion rate | ~3.4e24 G cm2/yr at 1 AU | prelim |
| Mass-to-flux lambda | 31.9 (initial cloud) -> ~2250 within 1 AU | prelim |
| Rotation | N = 0.09 (clean, through plateau end) / 0.24 (prelim, through front); period now ~7 yr (central object), 370 yr (first core); initial Omega t_ff = 0.02, beta_rot = 3.6e-4 | mixed |
| B-rho relation | kappa = 0.42 while flux-frozen (1e-16 to 1e-12 g/cm3), then a flat ~8 mG decoupling plateau (kappa ~ 0) | clean regime / prelim plateau |
| AMR / resolution | level 14, 2843 blocks (93M cells); finest cell 0.023 AU = 5.0 Rsun; 256^3 root on (0.47 pc)^3 | — |
| Current timestep | 6.8e-9 code = 1.0e4 s ~ 2.8 h physical per step (18-19 s wall on 5x H100) | — |

**The caveat (and why I trust the pipeline more, not less, after it).** While auditing the run we discovered that a super-time-stepping refactor (in production since sim cycle 71,000) was silently destroying the stored radiation field once per cycle; the matter-radiation coupling then rebuilt it every cycle, and in the reduced-speed-of-light bookkeeping that rebuild overcharged the gas by c/c_hat = 1000x — an artificial cooling channel worth ~30-70% of the compression heating during the runaway. We proved the thermal impact three independent ways (exact-EOS isentrope from the last clean state reaches ~2200 K at the front vs the contaminated readout of ~590 K; the restart transient magnitude matches the predicted 1000x aT^4 charge; the pre-bug radiation field sits exactly at aT^4). The bug is fixed and validated (gas evolution bit-identical in the fixed solver where radiation is absent; radiation now restart-continuous), everything through cycle 71,000 — including the plateau physics — is uncontaminated, and we are re-running the ~3-decade density runaway from the last clean checkpoint with the corrected binary (~2-4 days on 5x H100). Numbers marked "prelim" above will be finalized then; the corrected core will be hotter, so the second collapse should proceed closer to the textbook H2-dissociation channel, and the hotter core may partially re-couple the field via thermal (potassium) ionization — which directly matters for the fossil-flux budget.

**Fossil-field reading (preliminary).** The run measures the flux problem directly: B ~ rho^0.42 while coupled, then a ~8 mG "magnetic wall" once the ionization fraction floors — the material forming the star has already shed a factor ~70 in mass-to-flux relative to the cloud (lambda 31.9 -> ~2250 within 1 AU). For scale, a kG Ap-star field on 2 Rsun is ~6e25 G cm2, so the fossil budget question is quantitatively in play; the corrected (hotter) re-run will refine the plateau value and the degree of thermal re-coupling.

**How far we can push.** Resolution is provisioned through second-core formation (20 AMR levels -> 0.16 Rsun finest cell, Truelove maintained to ~1e-2 g/cm3): ~11 density e-folds beyond the front, ~6 GPU-days on 5x H100 at current throughput — less once the Hall-timestep optimization now being priced is enabled (removing Hall entirely measures 6.1x; a capped-Hall compromise is being evaluated). Beyond the second core the timestep collapses permanently and we hand over to the validated sink-particle module for the main accretion phase.

Happy to send the full flow-diagram/results document or any underlying plots.

Best regards,
Rahul
