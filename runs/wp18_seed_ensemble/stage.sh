#!/bin/bash
# WP-18 -- seed-variance ensemble. Stages one run dir per turbulence seed, varying ONLY
# problem/collapse_be/turb_seed. Everything else -- physics, resolution, and the CORRECTED
# turbulence IC (turb_ksample=k2, turb_nmodes=2048, turb_zeta=0.6667) -- is held fixed.
#
# WHY THIS IS THE GATING WP: every convergence work package (WP-1,2,3,7,8,9) compares two
# numbers and asks "did it move?". Without sigma there is no threshold, so those tests return
# numbers but not verdicts. WP-18 supplies the threshold.
#
# WHY IT MATTERS MORE AFTER WP-20: the corrected k^2 sampler concentrates the energy-bearing
# modes into the sparsely-populated low-k shells. Measured by Monte Carlo on vrms_analytic,
# realization scatter at the OLD nmodes=128 would have been 23.3% (vs 11.7% under the old flat
# sampler). nmodes=2048 brings it to 6.0%, but 6% is still the floor on every comparison.
#
# VEHICLE: the pre-collapse uniform leg (fhc_rootladder.in: refinement=none, tlim=1.0). That
# phase is 92% of the evolved time and only ~6% of the wall cost, and it is where the magnetic
# braking that sets flux retention actually happens. Measured cost: 2040 s at 256^3 on 4 GPUs
# (job 2433514), so ~1100 s on 8, plus ~3 min one-off IC generation at 2048 modes.
# A DEEP sub-ensemble (3-4 seeds to rho=1e-13) is a separate, later step -- sigma on the
# envelope first, because that is what the convergence WPs need.
#
# REQUIRES a binary that knows turb_ksample (post-5ebddce0). On an older binary the key is
# silently ignored and every seed runs the OLD IC -- the run.log banner is the check:
#   "k sampling          : k2  (dN/dk ~ k^2; alpha is the 3D PSD slope)"
set -e
BASE=/beegfs/u/bbg6470/athenapk/runs/wp18_seed_ensemble
DECK=/beegfs/u/bbg6470/athenapk/runs/root_ladder/fhc_rootladder.in
WRAP=/beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh          # WP-12-fixed wrapper
SEEDS="${SEEDS:-42 101 202 303 404 505 606 707 808 909 1010 1111}"

for S in $SEEDS; do
  D=$BASE/seed_$S
  mkdir -p $D
  sed "s|^turb_seed = .*|turb_seed = $S|" $DECK > $D/fhc_seed.in
  cp $WRAP $D/wrap_mod.sh; chmod +x $D/wrap_mod.sh
  grep -q "^turb_seed = $S\$" $D/fhc_seed.in || { echo "FAILED to set seed $S"; exit 1; }
done
echo "staged seeds: $SEEDS"
echo "verify only turb_seed differs between any two:"
echo "  diff $BASE/seed_42/fhc_seed.in $BASE/seed_101/fhc_seed.in"
