#!/bin/bash
#SBATCH --job-name=capprobe
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
# CONTROLLED one-variable AD-cap bisection: restart the SAME clean pre-blowup state
# (fc128c/out2.00010 = cyc1000, t=1.103539, 31 cycles before the cap=3 NaN at cyc1031)
# with a single AD cap value (env CAP) and run just past the blowup epoch (tlim=1.1050).
# STABLE => that cap survives cyc1031; NaN => still destabilizes. Non-chaining, short.
: "${CAP:?set CAP}" ; : "${TAG:?set TAG}"
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WD
WRAP=$WD/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
OUT=probe_$TAG; mkdir -p $OUT
RESTART=fc128c/parthenon.out2.00010.rhdf
echo "=== cap-probe TAG=$TAG CAP=$CAP $(date) restart=$RESTART tlim=1.1050 ===" >> $OUT/run.log
stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN -r $RESTART -t 01:20:00 \
  parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
  parthenon/mesh/numlevel=14 refinement/curr_max_level=4 parthenon/time/tlim=1.1050 parthenon/time/nlim=100000000 \
  parthenon/time/ncycle_out=10 parthenon/mesh/do_coalesced_comms=true \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 parthenon/output0/dt=-1 parthenon/output0/dn=1 \
  parthenon/output1/dn=100000 parthenon/output2/dn=100000 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/eta_ad_cap_code=$CAP \
  diffusion/ion_zeta=1.0e-16 -d $OUT >> $OUT/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> $OUT/run.log
