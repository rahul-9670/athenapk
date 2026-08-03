#!/bin/bash
#SBATCH --job-name=mg_prod_val
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/mg_prod_validate/%x_%j.out
#
# BOUNDED multigroup-RHD validation of the production FHC config (n_group=3 + Semenov opacity).
# Fresh from t=0, nlim=15 (self-terminates). SEPARATE dir + binary from prod_v9 (uses the
# multigroup build_gpu/bin/athenaPK; prod_v9 uses athenaPK_eos_v9, untouched). No chaining.

source ~/athenapk_env.sh
export OMP_NUM_THREADS=8
export OMPI_MCA_io=romio341
export PMIX_MCA_gds=hash

BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK   # multigroup + Semenov (NOT the prod binary)
WDIR=/beegfs/u/bbg6470/athenapk/runs/mg_prod_validate
WRAP=$WDIR/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

cd $WDIR
echo "binary:"; md5sum $BIN
echo "=== mg_prod_validate start $(date) job $SLURM_JOB_ID ===" >> $WDIR/run.log
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN -i fhc_mg.in \
  parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
  >> $WDIR/run.log 2>&1
RC=$?; echo "RUN_EXIT $RC $(date)" >> $WDIR/run.log
tail -5 $WDIR/run.log
