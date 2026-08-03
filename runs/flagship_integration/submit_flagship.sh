#!/bin/bash
#SBATCH --job-name=flagship_integ
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
#
# FULL-FLAGSHIP INTEGRATION TEST: the maximal config -- glmMHD + self-gravity multigrid + hydrogen
# hi-res EOS + gow17 chemistry + all three non-ideal MHD (RKL2 STS) + M1 RT (MULTIGROUP n_group=3
# + TABULATED opacity + rtsafe coupling) + PHYSICS-BASED AMR (jeans_nonideal current-sheet
# refinement) -- everything the flagship built, together, at production resolution (256^3
# numlevel=20). BOUNDED to tlim=1.07 (through collapse onset: AMR first refines + the coupling
# enters its stiff regime -- the maximal-stress integration point). One-shot (no chain).
# Confirms the whole stack coexists and holds. SEPARATE dir + fixed binary; prod_v9 untouched.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK   # rtsafe multigroup+tabulated (fixed)
WDIR=/beegfs/u/bbg6470/athenapk/runs/flagship_integration
WRAP=$WDIR/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
cd $WDIR
echo "binary:"; md5sum $BIN
# GPU-memory logger (jeans_nonideal at numlevel=20 could over-refine -> watch for allocator growth)
( while true; do echo "$(date +%s) $(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits | tr '\n' ' ')" >> $WDIR/gpumem.log; sleep 30; done ) &
GPULOG=$!
echo "=== flagship-integration start $(date) job $SLURM_JOB_ID ===" >> $WDIR/run.log
# do_coalesced_comms is ON: the coalesced+CT crash (asymmetric flxcor buffer set for CT's Edge
# flux field bnd_flux::Bf) is FIXED -- flxcor is now excluded from coalescing (individual path),
# validated coalesced==non-coalesced to within the GPU non-determinism floor (1.2e-10 < 5.8e-10).
# Recovers coalesced comms for the performance-critical halo + GMG exchange under CT.
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN -i fhc_flagship.in -t 01:50:00 \
  parthenon/time/tlim=1.07 parthenon/time/nlim=100000 parthenon/mesh/do_coalesced_comms=true \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=400 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 diffusion/ion_zeta=1.0e-16 \
  >> $WDIR/run.log 2>&1
RC=$?; kill $GPULOG 2>/dev/null; echo "RUN_EXIT $RC $(date)" >> $WDIR/run.log; tail -4 $WDIR/run.log
