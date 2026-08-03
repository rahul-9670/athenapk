#!/bin/bash
#SBATCH --job-name=mg_fix_regress
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/mg_fix_regress/%x_%j.out
set -o pipefail
#
# BOUNDED regression test of the rtsafe multigroup-coupling fix: the EXACT production config
# (256^3 numlevel=20, full physics, tabulated multigroup) on the FIXED binary, run to tlim=1.05 --
# past t~1.01 where the un-fixed run NaN-blew-up. One-shot (no chain). Confirms the fix at
# production fidelity: crosses the blowup point with zero non-convergence + finite fields.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK   # FIXED (rtsafe) multigroup+tabulated
WDIR=/beegfs/u/bbg6470/athenapk/runs/mg_fix_regress
WRAP=$WDIR/wrap_mod.sh
DECK=/beegfs/u/bbg6470/athenapk/runs/mg_prod_tab/fhc_mgtab.in
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
cd $WDIR
echo "binary:"; md5sum $BIN
echo "=== regress start $(date) job $SLURM_JOB_ID ===" >> $WDIR/run.log
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN -i $DECK -t 00:55:00 \
  parthenon/time/tlim=1.05 parthenon/time/nlim=100000 parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
  >> $WDIR/run.log 2>&1
RC=$?; echo "RUN_EXIT $RC $(date)" >> $WDIR/run.log; tail -4 $WDIR/run.log
