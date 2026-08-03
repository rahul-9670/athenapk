#!/bin/bash
#SBATCH --job-name=ct_nondet
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=00:20:00
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
COMMON="parthenon/mesh/nx1=64 parthenon/mesh/nx2=64 parthenon/mesh/nx3=64 \
  parthenon/meshblock/nx1=16 parthenon/meshblock/nx2=16 parthenon/meshblock/nx3=16 \
  parthenon/mesh/numlevel=4 parthenon/time/nlim=6 parthenon/time/ncycle_out=5 \
  parthenon/output1/single_precision_output=0 parthenon/output1/dn=1 parthenon/output1/id=g \
  parthenon/mesh/task_collection_timeout_in_seconds=1200 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 diffusion/ion_zeta=1.0e-16 \
  parthenon/mesh/do_coalesced_comms=true"
echo "=== ON run A ==="; rm -rf nd_a; mkdir nd_a
mpirun -n 2 $MCA $WRAP $BIN -i $DECK $COMMON -d nd_a > nd_a/log 2>&1; echo "A exit=$? driver=$(grep -c 'Driver completed' nd_a/log)"
echo "=== ON run B (identical) ==="; rm -rf nd_b; mkdir nd_b
mpirun -n 2 $MCA $WRAP $BIN -i $DECK $COMMON -d nd_b > nd_b/log 2>&1; echo "B exit=$? driver=$(grep -c 'Driver completed' nd_b/log)"
echo "=== DONE ==="
