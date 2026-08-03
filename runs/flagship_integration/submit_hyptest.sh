#!/bin/bash
#SBATCH --job-name=hyptest
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WD
WRAP=$WD/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
# 128^3 collapse to onset (t=1.06, past the deep-core regrid), 2 GPUs. Compare block growth for
# type=jeans (mg_prod_tab, reached 1st core) vs jeans_nonideal (flagship, OOM'd). Isolates the
# current-sheet criterion's memory cost.
COMMON="parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
  parthenon/mesh/numlevel=14 parthenon/time/tlim=1.06 parthenon/time/nlim=100000 parthenon/time/ncycle_out=10 \
  parthenon/mesh/do_coalesced_comms=true parthenon/mesh/task_collection_timeout_in_seconds=1800 \
  parthenon/output0/dt=-1 parthenon/output0/dn=1 parthenon/output1/dn=100000 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 diffusion/ion_zeta=1.0e-16"
for RT in jeans jeans_nonideal; do
  echo "=== refinement/type=$RT 128^3 $(date) ==="; rm -rf ht_$RT; mkdir ht_$RT
  mpirun -n 2 $MCA $WRAP $BIN -i fhc_flagship.in $COMMON refinement/type=$RT -d ht_$RT > ht_$RT/log 2>&1
  echo "$RT: exit=$? maxblk=$(awk 'NR>2{print $4}' ht_$RT/*.hst 2>/dev/null|sort -n|tail -1) OOM=$(grep -ciE 'failed to allocate' ht_$RT/log) lastt=$(grep '^cycle=' ht_$RT/log|tail -1|grep -oE 'time=[0-9.e+-]+')"
done
echo "DONE"
