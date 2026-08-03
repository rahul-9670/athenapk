#!/bin/bash
#SBATCH --job-name=fc128
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
# Full flagship -> FIRST CORE at 128^3 on 2 GPUs (instant slot; compute is infinite so slow is fine).
# curr_max_level=4 bounds the current-sheet refinement (the deep-core block explosion) while Jeans
# owns the deep collapse. Self-chaining + first-core density stop. Validates the fix AND reaches
# first core with the full stack. Binary 1603b221.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WD
WRAP=$WD/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python; RHO0=5.467e-19; STOP_CGS=1.0e-12
N=$(cat fc128b_n 2>/dev/null||echo 0); N=$((N+1)); echo $N > fc128b_n
[ -f STOP_128b ] && { echo "STOP -> exit"; exit 0; }
NEWEST=$(ls -t fc128b/parthenon.out1.*.phdf 2>/dev/null|head -1)
if [ -n "$NEWEST" ]; then
  RM=$($PY - "$NEWEST" <<'PYEOF'
import sys,h5py,numpy as np
try:
 with h5py.File(sys.argv[1],"r") as h: print("%.6e"%float(np.array(h["prim"][:,0,...]).max()))
except: print("0.0")
PYEOF
)
  [ "$($PY -c "print(1 if $RM*$RHO0>=$STOP_CGS else 0)")" = "1" ] && { echo "FIRST CORE rho=$RM -> STOP"; echo "first core $RM $(date)" > STOP_128b; exit 0; }
fi
LATEST=$(ls -t fc128b/parthenon.out2.*.rhdf 2>/dev/null|head -1)
[ $N -ge 2 ] && [ -z "$LATEST" ] && { echo "no restart -> STOP"; echo "no-restart" > STOP_128b; exit 0; }
[ $N -lt 40 ] && sbatch --dependency=afterany:$SLURM_JOB_ID $0 && echo "successor queued"
if [ -n "$LATEST" ]; then RA="-r $LATEST"; else RA="-i fhc_flagship.in"; fi
mkdir -p fc128b
echo "=== fc128 slot $N $(date) 128^3 curr_max_level=4 ===" >> fc128b/run.log
stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN $RA -t 07:45:00 \
  parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
  parthenon/mesh/numlevel=14 refinement/curr_max_level=4 parthenon/time/tlim=1.30 parthenon/time/nlim=100000000 \
  parthenon/time/ncycle_out=20 parthenon/mesh/do_coalesced_comms=true \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 parthenon/output0/dt=-1 parthenon/output0/dn=1 \
  parthenon/output1/dn=200 parthenon/output2/dn=100 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/ion_zeta=1.0e-16 -d fc128b >> fc128b/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> fc128b/run.log
