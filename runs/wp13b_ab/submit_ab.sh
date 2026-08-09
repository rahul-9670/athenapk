#!/bin/bash
#SBATCH --job-name=wp13ab
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:h100:4
#SBATCH --cpus-per-task=8
#SBATCH --time=02:30:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp13b_ab/%x_%j.out
set -o pipefail
#
# WP-13b LONG-BASELINE A/B — does the DiodeBC radiation-ghost fix change the SCIENCE?
#
# THE QUESTION. WP-13b proved `DiodeBC` packed only "cons", so the M1 radiation moments had no
# domain boundary condition; it is fixed in 967fced6. The fix is RESULT-CHANGING even for a fresh
# run, but the only measurement so far is 31 cycles deep: Er_tot +5.5e-05, mass -7.0e-10,
# ME -1.6e-08, KE -1.5e-06, rho_max -5.0e-07. Those are tiny. They are ALSO meaningless as a
# bound, because the defect is a BOUNDARY FLUX error: it injects/removes radiation energy through
# the domain face every step, so its effect ACCUMULATES. A production member runs 1e4-1e5 cycles.
# Nothing measured so far says whether 5.5e-05 at cycle 31 becomes 5e-05 or 5e-01 at cycle 1e4.
#
# THIS TEST ANSWERS THAT AND NOTHING ELSE. It does not re-validate the fix (job 2491892 did that).
#
# THE DESIGN — THREE LEGS, and the third is the whole point.
#   old_a : binary 84a6d248 (what the 24-member ensemble actually ran)
#   old_b : binary 84a6d248, IDENTICAL to old_a
#   new   : binary 967fced6 (the fix)
# old_a vs old_b measures the GPU/MPI NON-DETERMINISM FLOOR at this depth. Without it the test is
# uninterpretable, because 4-rank GPU reductions are not order-deterministic and ~1e4 cycles of a
# collapsing turbulent flow amplify a last-bit difference into a visible one. The claim "the fix
# does not move the science" is only meaningful as "old-vs-new sits inside old-vs-old".
# This is the same discriminating design that made the original WP-13b readable; see
# docs/validation/WP13b_restart_gpu_amr.md.
#
# ALL THREE LEGS ARE FRESH FROM t=0 AND RESTART-FREE (single slot, no chaining). That is
# deliberate. A restart inside the old-binary leg would apply the zeroed-ghost corruption on top
# of the fresh-run difference and CONFOUND the two effects, and separating them is exactly what
# this test exists to do.
#
# DECK = ensemble point000's fhc_ens.in, unmodified, with three CLI overrides:
#   single_precision_output=0  -- the deck ships =1 (~1e-7 relative). The 31-cycle bulk signal is
#                                ~1e-6, so single precision would sit right on top of the effect
#                                being measured and could hide it entirely.
#   variables=...              -- grav.phi / rad.Er_g1 / rad.Er_g2 dropped. Doubling the precision
#                                doubles the file size, and 3 legs x 17 snapshots x 4 GB would be
#                                ~400 GB against ~2.8 TB of free quota. gravity is exonerated on
#                                this deck (restart divergence 0.88x its OWN floor) and the
#                                per-group energies are not science observables.
#   hst dt=0.001               -- 5x finer than the deck's 0.005. The .hst is FULL PRECISION and
#                                tiny, so it, not the phdf, is the primary comparison vehicle for
#                                mass / KE / ME / tot-E. The phdf is only for rho_max, mu_core and
#                                the field structure at the matched epoch.
#
# WALLTIME. point000 reached rho = 2.00e-12 cgs (20x rhocrit, PAST the 1e-12 measurement epoch)
# in one 1h45m 4-GPU slot. Internal -t 02:15:00 therefore passes the epoch with ~30 % margin while
# still CAPPING how deep the run goes -- depth is what drives block count, and the measured memory
# law (D1: 0.159 GiB x blocks/rank) puts the 1373-block endpoint at 54.6 GiB of 79.2 on 4 GPUs.
# Running longer buys nothing here and risks OOM.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

: "${AB_LEG:?set AB_LEG (old_a|old_b|new)}"
R=/beegfs/u/bbg6470/athenapk/runs/wp13b_ab
D=$R/$AB_LEG
case "$AB_LEG" in
  old_a|old_b) BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_PRESERVED_84a6d248 ;;
  new)         BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_PRESERVED_967fced6 ;;
  *)           echo "unknown AB_LEG=$AB_LEG"; exit 2 ;;
esac

mkdir -p $D
SRC=/beegfs/u/bbg6470/athenapk/runs/ensemble/design01/point000
cp -n $SRC/fhc_ens.in $D/ 2>/dev/null
cp -n $SRC/units.json $D/ 2>/dev/null
cp -n $SRC/wrap_mod.sh $D/ 2>/dev/null
WRAP=$D/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

echo "=== WP-13b A/B leg=$AB_LEG job=$SLURM_JOB_ID $(date) host=$(hostname) ==="
echo "binary md5 $(md5sum $BIN | cut -c1-8)  path $BIN"

cd $D
# PROD flags verbatim from ensemble/submit_point.sh -- the physics must be the production physics,
# or this measures something other than what production does.
stdbuf -oL -eL mpirun -n 4 $MCA $WRAP $BIN -i $D/fhc_ens.in -t 02:15:00 \
  parthenon/mesh/do_coalesced_comms=true diffusion/integrator=rkl2 \
  diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=1000 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 \
  parthenon/output2/dn=100000 parthenon/output1/dn=50 \
  parthenon/output1/single_precision_output=0 \
  "parthenon/output1/variables=prim, rad.Er, rad.Fr1, rad.Fr2, rad.Fr3" \
  parthenon/output0/dt=0.001 \
  >> $D/run.log 2>&1
echo "RUN_EXIT $? $(date)" | tee -a $D/run.log
# output2/dn=100000 suppresses restart dumps on purpose: these legs must never restart, and a
# 4 GB rhdf every 250 cycles is pure quota cost for a file this test must not use.
grep -c "^cycle=" $D/run.log | sed 's/^/cycles_done=/'
