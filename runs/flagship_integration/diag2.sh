#!/bin/bash
#SBATCH --job-name=diag2
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=00:40:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WD; WRAP=$WD/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

# Shared config: restart cycle 200 (t=1.100); advance 3 REAL cycles (nlim=203); print every cycle.
# NO big science/restart dumps (dn huge). Only the DIAG_* stdout lines matter.
COMMON="parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
  parthenon/mesh/numlevel=14 refinement/curr_max_level=4 parthenon/time/tlim=100.0 parthenon/time/nlim=203 \
  parthenon/time/ncycle_out=1 parthenon/mesh/do_coalesced_comms=true \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 \
  parthenon/output0/dt=-1 parthenon/output0/dn=100000 parthenon/output1/dn=100000 parthenon/output2/dn=100000 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/ion_zeta=1.0e-16"

for CD in ct glm; do
  OUT=diag_${CD}2; rm -rf $OUT; mkdir -p $OUT
  echo "=== DIAG2 ${CD^^}: restart fc128b/out2.00002 cyc200 t=1.100, nlim=203 (3 real cyc), binary=$(md5sum $BIN|cut -c1-8) $(date) ===" >> $OUT/run.log
  stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN -r fc128b/parthenon.out2.00002.rhdf -t 00:18:00 \
    $COMMON hydro/divergence_control=${CD} hydro/ct_emf=gs05 -d $OUT >> $OUT/run.log 2>&1
  echo "RUN_EXIT $? $(date)" >> $OUT/run.log
done
echo "ALL DONE $(date)"
