#!/bin/bash
#SBATCH --job-name=flagship8
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gres=gpu:h100:8
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
# AMR-verification rerun: SAME flagship config (CT + multigroup tabulated + chem + non-ideal +
# gravity + jeans_nonideal + coalesced) but on 8 GPUs (60% more HBM headroom than the 5-GPU OOM run)
# and with per-CYCLE hst logging (nbtotal every cycle) + GPU-mem trace. Tests whether the
# wsec_AMR explosion + OOM is memory-pressure (goes away with headroom) or real volume-filling
# over-refinement by curr_nsheet on the turbulent field (persists -> needs curr_nsheet fix).
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WDIR=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WDIR
WRAP=$WDIR/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
echo "binary:"; md5sum $BIN
( while true; do echo "$(date +%s) $(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits|tr '\n' ' ')" >> $WDIR/gpumem8.log; sleep 15; done ) & GPULOG=$!
echo "=== flagship8 start $(date) job $SLURM_JOB_ID ===" >> $WDIR/run8.log
stdbuf -oL -eL mpirun -n 8 $MCA $WRAP $BIN -i fhc_flagship.in -t 01:20:00 \
  parthenon/time/tlim=1.07 parthenon/time/nlim=100000 parthenon/mesh/do_coalesced_comms=true \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 parthenon/output0/dt=-1 parthenon/output0/dn=1 \
  parthenon/time/ncycle_out=1 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/ion_zeta=1.0e-16 >> $WDIR/run8.log 2>&1
RC=$?; kill $GPULOG 2>/dev/null; echo "RUN_EXIT $RC $(date)" >> $WDIR/run8.log; tail -5 $WDIR/run8.log
