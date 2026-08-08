#!/bin/bash
#SBATCH --job-name=wp13b
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:h100:4
#SBATCH --cpus-per-task=8
#SBATCH --time=00:40:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp13b_%j.out
set -o pipefail
#
# WP-13b — restart reproducibility in the PRODUCTION configuration (GPU, 4 ranks, AMR, rkl2).
#
# WHY THIS EXISTS. WP-13 is closed and its fix is real, but it was validated on `build_cpu`, one
# rank, `diffusion/integrator = unsplit`, 12 cycles. Production is GPU + 4 MPI ranks + RKL2 STS.
# Nothing had tested restart in THAT configuration. The trigger was an observation on ensemble
# member point012 (2026-08-08): a fresh leg and a leg restarted from the t=0 restart file reached
# cycle 200 with t = 1.0262248 vs 0.9922842 and dt = 1.940e-4 vs 9.970e-5 -- a factor-2 dt
# difference. That is NOT by itself evidence of a restart bug, because GPU reductions are not
# order-deterministic and 200 cycles of a collapsing turbulent flow amplify a last-bit difference.
#
# THE DISCRIMINATING DESIGN. Three legs, same deck, same binary, same rank count:
#   fresh_a : -i fhc.in, nlim=60                       reference
#   fresh_b : -i fhc.in, nlim=60                       IDENTICAL to fresh_a
#   split   : -i fhc.in nlim=30, then -r <rst> nlim=60  restart at the midpoint
# fresh_a vs fresh_b measures the GPU/MPI non-determinism FLOOR. fresh_a vs split measures
# restart divergence. The claim "restart is as reproducible as this code is" holds iff the
# second difference is within the first. Without leg B the test cannot distinguish a restart
# defect from ordinary non-determinism -- which is exactly the ambiguity point012 left open.
#
# 60 cycles rather than WP-13's 12: numlevel=2 AMR must actually regrid inside the window, and
# the restart must land after at least one regrid, or the AMR path is untested.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

# PARAMETERIZED so the isolation experiments reuse this exact harness rather than a copy that
# can drift. Defaults reproduce the production configuration.
#   WP13B_DIR   = run directory (each experiment gets its own, so evidence is never clobbered)
#   WP13B_DIFF  = diffusion/integrator   (rkl2 = production; unsplit = what WP-13 actually tested)
#   WP13B_NLIM  = final cycle            (60 = default; 31 = one step past the restart)
#   WP13B_SPLIT = restart cycle          (30)
#   WP13B_DECK  = the input deck to test. DEFAULT (wp13_restart/straight/fhc.in) is the GRAY
#                 radiation path -- it sets no `n_group`, so n_group=1 and MatterCoupling's
#                 plain linearized Newton runs. The FLAGSHIP and all 69 production decks set
#                 `n_group = 3`, which dispatches to MatterCouplingMultigroup and its SAFEGUARDED
#                 rtsafe solve instead. Those are DIFFERENT SOLVERS, so a result on one says
#                 nothing about the other; use ensemble/design01/point000/fhc_ens.in to test the
#                 solver production actually runs.
R=${WP13B_DIR:-/beegfs/u/bbg6470/athenapk/runs/wp13b_gpu_restart}
DECK=${WP13B_DECK:-/beegfs/u/bbg6470/athenapk/runs/wp13_restart/straight/fhc.in}
DIFFINT=${WP13B_DIFF:-rkl2}
NLIM=${WP13B_NLIM:-60}
SPLIT=${WP13B_SPLIT:-30}
# rkl2_freeze_eta is GATED: the code hard-fails with
#   "diffusion/rkl2_freeze_eta=true requires diffusion/eta_cache=true and diffusion/integrator=rkl2"
# so an unsplit leg MUST also pass WP13B_FREEZE=false. Experiment A (job 2490926) aborted 134 on
# all three legs for exactly this reason.
FREEZE=${WP13B_FREEZE:-true}
EXTRA=${WP13B_EXTRA:-}
mkdir -p $R/fresh_a $R/fresh_b $R/split
DECKDIR=$(dirname $DECK)
for d in fresh_a fresh_b split; do
  cp -n $DECK $R/$d/fhc.in 2>/dev/null
  cp -n $DECKDIR/units.json $R/$d/ 2>/dev/null
done
cp -n /beegfs/u/bbg6470/athenapk/runs/wp13b_gpu_restart/wrap_mod.sh $R/ 2>/dev/null
cp -n /beegfs/u/bbg6470/athenapk/runs/wp13b_gpu_restart/compare_wp13b.py $R/ 2>/dev/null
# Same pinned binary the ensemble ran, so this statement covers the ensemble's own restarts.
# WP13B_BIN lets a leg test a CANDIDATE binary. Default stays the binary the ensemble ran, so
# every result already recorded remains reproducible with the same command.
BIN=${WP13B_BIN:-/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_PRESERVED_84a6d248}
WRAP=$R/wrap_mod.sh
# RAD_PRINT_NSUB / RAD_DISABLE_TRANSPORT are radiation_moments.cpp's own getenv diagnostics
# (nsub = ceil(dt/dt_rad) is an INTEGER, so it is a candidate discrete branch across a restart).
# They must be forwarded explicitly: mpirun does not inherit the submitter's environment.
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
[ -n "$RAD_PRINT_NSUB" ] && MCA="$MCA -x RAD_PRINT_NSUB"
[ -n "$RAD_DISABLE_TRANSPORT" ] && MCA="$MCA -x RAD_DISABLE_TRANSPORT"
# Production physics flags, verbatim from runs/ensemble/submit_point.sh. output2/dn=30 puts a
# restart file exactly at the split point; output1/dn=100000 leaves only the final phdf, which is
# the artifact being compared.
#
# tlim=1000 IS LOAD-BEARING. The deck's own `tlim = 1.5` is reached at cycle 22, BEFORE the
# nlim=30 split point. The first version of this script omitted it and the test was VACUOUS:
# split leg 1 ran to tlim at cycle 22, and the "restart" then started at cycle 22 already at
# tlim and completed in 0.108 s having advanced ZERO cycles. All three legs were therefore just
# the same straight run, and they compared bit-identical for a trivial reason. Only nlim may
# terminate these legs.
PROD="parthenon/mesh/do_coalesced_comms=true diffusion/integrator=$DIFFINT \
diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=1000 \
diffusion/rkl2_freeze_eta=$FREEZE diffusion/eta_ohm_cap_code=0.1 parthenon/output2/dn=$SPLIT \
parthenon/time/tlim=1000.0 $EXTRA"

echo "=== WP-13b dir=$R diff=$DIFFINT nlim=$NLIM split=$SPLIT $(date) host=$(hostname) job=$SLURM_JOB_ID ==="
echo "binary md5 $(md5sum $BIN | cut -c1-8)"

run() {  # run <dir> <args...>
  local d=$1; shift
  ( cd $R/$d && stdbuf -oL -eL mpirun -n 4 $MCA $WRAP $BIN "$@" $PROD ) >> $R/$d/run.log 2>&1
  echo "$d exit=$?" >> $R/status
}

: > $R/status
run fresh_a -i $R/fresh_a/fhc.in parthenon/time/nlim=$NLIM
run fresh_b -i $R/fresh_b/fhc.in parthenon/time/nlim=$NLIM
run split   -i $R/split/fhc.in   parthenon/time/nlim=$SPLIT

RST=$(ls -1t $R/split/parthenon.out2.*.rhdf 2>/dev/null | head -1)
echo "restart_from=$RST" >> $R/status
# A restart from the t=0 file would silently re-run from scratch and make the whole test vacuous.
# That is the exact trap that hit ensemble members point012/013/014 on 2026-08-08.
case "$RST" in
  *out2.00000.rhdf) echo "ABORT: newest restart is the t=0 file -- output2/dn did not fire" >> $R/status; exit 1 ;;
  "")               echo "ABORT: no restart file written" >> $R/status; exit 1 ;;
esac
run split -r $RST parthenon/time/nlim=$NLIM

# POSITIVE CHECK. Exit 0 is not evidence the restart leg did anything -- the first version of this
# test passed while the restart advanced ZERO cycles. Every leg must END at cycle $NLIM, and the
# restart leg must have STARTED at cycle $SPLIT.
for d in fresh_a fresh_b split; do
  echo "$d last_cycle=$(grep -oE '^cycle=[0-9]+' $R/$d/run.log | tail -1 | cut -d= -f2)" >> $R/status
done
echo "split restart_started_at_cycle=$(awk '/Starting up hydro driver/{n++} n==2' $R/split/run.log \
      | grep -m1 -oE '^cycle=[0-9]+' | cut -d= -f2)" >> $R/status

echo DONE >> $R/status
cat $R/status
