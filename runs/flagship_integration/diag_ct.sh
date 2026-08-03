#!/bin/bash
#SBATCH --job-name=diagct
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
OUT=diag_ct; mkdir -p $OUT
echo "=== DIAG CT: eta_A/M(Y0) instrumentation, restart fc128b/out2.00002 t=1.100, nlim=2, binary=$(md5sum $BIN|cut -c1-8) $(date) ===" >> $OUT/run.log
stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN -r fc128b/parthenon.out2.00002.rhdf -t 00:35:00 \
  parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
  parthenon/mesh/numlevel=14 refinement/curr_max_level=4 parthenon/time/tlim=100.0 parthenon/time/nlim=2 \
  parthenon/time/ncycle_out=1 parthenon/mesh/do_coalesced_comms=true \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 parthenon/output0/dt=-1 parthenon/output0/dn=100000 \
  parthenon/output1/dn=100000 parthenon/output2/dn=100000 \
  hydro/divergence_control=ct hydro/ct_emf=gs05 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/ion_zeta=1.0e-16 -d $OUT >> $OUT/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> $OUT/run.log
