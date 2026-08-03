#!/bin/bash
#SBATCH --job-name=ncu_v8
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/gpu_profile_v8/ncu_%j.out
set -o pipefail
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
export KOKKOS_TOOLS_LIBS=/beegfs/u/bbg6470/athenapk/runs/gpu_profile_v8/libkp_nvtx.so
export KOKKOS_PROFILE_LIBRARY=$KOKKOS_TOOLS_LIBS
# After nsys names the top kernels, set NCU_SKIP to land past the warmup step and NCU_COUNT
# to a small number (ncu replays each kernel several passes -> slow). Optionally restrict to
# the top NVTX labels by editing wrap_ncu.sh to add e.g.
#   --nvtx-include "PrecomputeNonidealEta/" --nvtx-include "Ambipolar X1 fluxes/" ...
export NCU_SKIP=${NCU_SKIP:-300}
export NCU_COUNT=${NCU_COUNT:-40}

BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_eos_v8_auditfixes
WDIR=/beegfs/u/bbg6470/athenapk/runs/gpu_profile_v8
WRAP=$WDIR/wrap_ncu.sh
RST=$WDIR/restart_c1000.rhdf
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR -x KOKKOS_TOOLS_LIBS -x KOKKOS_PROFILE_LIBRARY -x NCU_SKIP -x NCU_COUNT"
cd $WDIR
mkdir -p ncu
echo "=== ncu profile start $(date) job $SLURM_JOB_ID on $(hostname) ==="
md5sum $BIN
# nlim=1002: rank0 runs slowly under ncu replay; a couple steps is plenty (kernel metrics
# are structural). Other ranks block at MPI collectives while rank0 replays -- expected.
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN -r $RST -t 01:50:00 \
  parthenon/time/nlim=1002 \
  parthenon/output1/dn=1000000 parthenon/output2/dn=1000000 \
  parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
  2>&1 | tee $WDIR/ncu_run.log
echo "=== ncu profile done $(date) rc=$? ==="
ls -la $WDIR/ncu/
# Summarize: ncu -i ncu/ncu_rep.ncu-rep --page raw --csv > ncu/summary.csv
