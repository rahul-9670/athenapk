#!/bin/bash
#SBATCH --job-name=ct_gpu_on
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/ct_coalesced_fix/%x_%j.out
set -o pipefail
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/ct_coalesced_fix; cd $WD
WRAP=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
DECK=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/fhc_flagship.in
# REAL physics (non-ideal on, normal creduc) so the collapse+AMR dynamics match the flagship;
# coalesced ON so it runs FAST until the crash (~collapse onset). Full flagship stack.
COMMON="parthenon/mesh/nx1=64 parthenon/mesh/nx2=64 parthenon/mesh/nx3=64 \
  parthenon/meshblock/nx1=16 parthenon/meshblock/nx2=16 parthenon/meshblock/nx3=16 \
  parthenon/mesh/numlevel=4 parthenon/time/nlim=200 parthenon/time/ncycle_out=5 \
  parthenon/mesh/task_collection_timeout_in_seconds=1200 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/ion_zeta=1.0e-16"
echo "binary:"; md5sum $BIN
echo "=== ON (coalesced, real physics) $(date) ==="; rm -rf gpu_on; mkdir gpu_on
mpirun -n 2 $MCA $WRAP $BIN -i $DECK $COMMON parthenon/mesh/do_coalesced_comms=true -d gpu_on > gpu_on/log 2>&1
echo "ON exit=$? cyc=$(grep -cE '^cycle=' gpu_on/log) recv-buf=$(grep -c 'receiving buffer' gpu_on/log) driver=$(grep -c 'Driver completed' gpu_on/log)"
echo "=== DONE $(date) ==="
