#!/bin/bash
#SBATCH --job-name=ws1_gpusmoke
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --output=%x_%j.out
# WS-1 GPU sink smoke test: run twobody.in (2 sinks, N-body + MPI gather + accretion-off) on
# 2 H100s for a few orbits to confirm the sink host-side deep_copy/MPI path runs on GPU.
set -o pipefail
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WRAP=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io"
WDIR=/beegfs/u/bbg6470/athenapk/runs/validation_ws1/gpu_smoke
rm -rf "$WDIR" && mkdir -p "$WDIR" && cd "$WDIR"
echo "binary:"; md5sum $BIN
# 3 orbits, coarse, to check the sink N-body runs on GPU without crash and energy is sane.
stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN -i ../twobody.in \
  parthenon/time/tlim=37.7 parthenon/output0/dt=6.28 >run.log 2>&1
echo "EXIT $? $(date)"; grep -iE "Driver completed|error|abort|nan|PARTHENON ERROR" run.log | grep -ivE "residual" | tail -4
