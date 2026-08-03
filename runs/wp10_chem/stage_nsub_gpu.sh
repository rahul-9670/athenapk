#!/bin/bash
# WP-10 follow-up — is raising `chemistry/nsub_max` FREE on the PRODUCTION configuration?
#
# WHY THE EXISTING MEASUREMENT IS NOT ENOUGH. runs/wp10_chem/ measured nsub_max 400 vs 3200 vs
# 100000 at identical wall time (380 s every leg) and concluded "free". That was **32^3 on CPU**,
# where the chemistry is a small fraction of a step dominated by M1 radiation, multigrid gravity
# and the non-ideal RKL2 stack. Two things make it non-transferable:
#
#   1. COST MODE. The sub-cycler is a SERIAL per-cell loop inside a Kokkos kernel. On a GPU the
#      cells in a warp advance in lockstep, so a warp costs its SLOWEST cell: if one cell wants
#      10^4 sub-steps and its neighbours want 10, all of them pay 10^4. On CPU (one cell per
#      thread, no lockstep) that penalty simply does not exist, so the CPU test cannot see it.
#   2. REGIME. Production reaches rho/rho_0 ~ 1e10 with a chemistry x_e that collapses through
#      the C->CO ionization minimum. The number of sub-steps the accuracy criterion asks for is a
#      function of the local state, and the smoke deck never gets there.
#
# WHAT COULD ALSO CHANGE PHYSICALLY. x_e feeds eta_A through `ambipolar_coeff = ionization_chem`.
# WP-22 measured the equilibrium ceiling suppressing the chemistry branch by 1486-2678x IN THE
# FIRST CORE, so x_e is largely masked there -- but in the ENVELOPE, which is where the magnetic
# braking that sets flux retention actually happens, the chemistry branch may well be the binding
# one. So this is not only a cost question.
#
# DESIGN. Production deck, production binary, 128^3 (the rung r128_sw runs in ~5 min on 3 GPUs),
# three legs differing ONLY in nsub_max:
#     400    production
#     4000   10x
#     40000  100x
# Judge (a) wsec_step -- the cost answer; (b) MEtor/MEpol at t = 0.90 -- the physics answer, read
# BEFORE the singular endpoint per WP-7; (c) the B10 floored-cell warning, which build_gpu_v3
# carries, so we learn how much of the domain is still cap-limited at each setting.
#
# ACCEPTANCE. Raise it if the cost is small (say < 10 % on wsec_step) AND the physics shift is
# well below sigma = 16 %. If cost is large, the answer is a middling nsub_max, not 400 and not
# 100000. If the physics shift is LARGE, that is a finding in its own right -- it would mean
# production's chemistry has never been converged and every flux-retention number inherits it.
#
# Usage: ./stage_nsub_gpu.sh        -> print
#        ./stage_nsub_gpu.sh go     -> submit
set -o pipefail
LAD=/beegfs/u/bbg6470/athenapk/runs/root_ladder
HERE=/beegfs/u/bbg6470/athenapk/runs/wp10_chem
WRAPSRC=/beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh
BIN="${BIN:-/beegfs/u/bbg6470/athenapk/build_gpu_v3/bin/athenaPK}"
[ -r "$WRAPSRC" ] || { echo "FATAL: no GPU-pinning wrapper at $WRAPSRC"; exit 1; }
[ -x "$BIN" ]     || { echo "FATAL: no binary at $BIN (build_gpu_v3 not built yet?)"; exit 1; }
echo "binary: $(md5sum $BIN)"
echo "deck:   $(md5sum $LAD/fhc_rootladder.in)"

for NS in 400 4000 40000; do
  DIR=$HERE/gpu_ns$NS
  mkdir -p $DIR
  install -m 755 $WRAPSRC $DIR/wrap_mod.sh || { echo "FATAL: cannot install wrapper"; exit 1; }
  CMD="sbatch --job-name=ns$NS --nodes=1 --ntasks=3 --gres=gpu:h100:3 \
--export=ALL,NX=128,RUNDIR=$DIR,NRANK=3,BIN=$BIN,OV=chemistry/nsub_max=$NS \
$LAD/submit_root.sh"
  if [ "${1:-}" = "go" ]; then echo "+ $CMD"; eval $CMD; else echo "would run: $CMD"; fi
done
