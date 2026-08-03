#!/bin/bash
#SBATCH --job-name=nsys_v8
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/gpu_profile_v8/nsys_%j.out
set -o pipefail
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
# Kokkos->NVTX connector so nsys shows physics kernel LABELS (freeze eta, rad flux, MG, ...)
export KOKKOS_TOOLS_LIBS=/beegfs/u/bbg6470/athenapk/runs/gpu_profile_v8/libkp_nvtx.so
export KOKKOS_PROFILE_LIBRARY=$KOKKOS_TOOLS_LIBS

BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_eos_v8_auditfixes
WDIR=/beegfs/u/bbg6470/athenapk/runs/gpu_profile_v8
WRAP=$WDIR/wrap_nsys.sh
RST=$WDIR/restart_c1000.rhdf
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR -x KOKKOS_TOOLS_LIBS -x KOKKOS_PROFILE_LIBRARY"
cd $WDIR
mkdir -p nsys
echo "=== nsys profile start $(date) job $SLURM_JOB_ID on $(hostname) ==="
md5sum $BIN
# Restart at cycle 1000 (first-core, 17-19 STS substeps). Run 3 steps (1001 warmup +
# 1002/1003 representative). Outputs disabled (huge dn). Physics CLI mirrors prod submit.sh.
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN -r $RST -t 00:50:00 \
  parthenon/time/nlim=1003 \
  parthenon/output1/dn=1000000 parthenon/output2/dn=1000000 \
  parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
  2>&1 | tee $WDIR/nsys_run.log
echo "=== nsys profile done $(date) rc=$? ==="
ls -la $WDIR/nsys/
