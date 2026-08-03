#!/bin/bash
#SBATCH --job-name=stsopt_ctrl
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/stsopt_ab/%x_%j.out
set -o pipefail
# Correctness with a NON-DETERMINISM CONTROL: run base TWICE (A,B) + cand once (C), each
# 1 cycle from the same seed. A-vs-B = the GPU non-determinism floor (Poisson SUM-reductions
# / atomics); A-vs-C = base vs STS-opt. Result-preserving <=> A-vs-C ~ A-vs-B.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
W=/beegfs/u/bbg6470/athenapk/runs/stsopt_ab; WRAP=$W/wrap_mod.sh
CLI="parthenon/mesh/do_coalesced_comms=true diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 parthenon/time/nlim=251 parthenon/output1/dn=1 parthenon/output2/dn=1000000"
cd $W
run() { # $1=tag $2=binary
  local OUT=$W/ctrl_$1; rm -rf $OUT; mkdir -p $OUT
  stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $2 -r $W/small_seed.rhdf -d $OUT $CLI > $OUT/log.txt 2>&1
  echo "$1 exit=$? $(grep '^cycle=' $OUT/log.txt | tail -1)"
}
run A /beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_stsopt_base
run B /beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_stsopt_base
run C /beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_stsopt_cand
echo "DONE $(date)"
