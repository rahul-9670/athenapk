#!/bin/bash
#SBATCH --job-name=fc128d
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
# Full flagship -> FIRST CORE at 128^3 on 2 GPUs, WITH the ambipolar eta_A cap dt-fix.
# Root cause of the 7e-8 dt-wall = UNCAPPED eta_A (~1/rho^2 runaway in low-rho refined cells),
# the one non-ideal term lacking the cap Ohm/Hall have. New key diffusion/eta_ad_cap_code=30.0
# (analytic optimum <3.33 crossover, depth-independent => ~80x dt gain to the Hall limit while
# preserving real AD to ~20x above the dense-core value). Restarts from the fc128b deep-collapse
# state (rho~6.9e-14) so the fix's dt jump is directly visible and the expensive runup is reused.
# New GPU binary carries the eta_ad_cap mechanism. Self-chaining + first-core density stop.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WD
WRAP=$WD/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python; RHO0=5.467e-19; STOP_CGS=1.0e-12
N=$(cat fc128d_n 2>/dev/null||echo 0); N=$((N+1)); echo $N > fc128d_n
[ -f STOP_128d ] && { echo "STOP -> exit"; exit 0; }
NEWEST=$(ls -t fc128d/parthenon.out1.*.phdf 2>/dev/null|head -1)
if [ -n "$NEWEST" ]; then
  RM=$($PY - "$NEWEST" <<'PYEOF'
import sys,h5py,numpy as np
try:
 with h5py.File(sys.argv[1],"r") as h: print("%.6e"%float(np.array(h["prim"][:,0,...]).max()))
except: print("0.0")
PYEOF
)
  [ "$($PY -c "print(1 if $RM*$RHO0>=$STOP_CGS else 0)")" = "1" ] && { echo "FIRST CORE rho=$RM -> STOP"; echo "first core $RM $(date)" > STOP_128d; exit 0; }
fi
# First slot restarts from the fc128b deep state; later slots chain from fc128d.
LATEST=$(ls -t fc128d/parthenon.out2.*.rhdf 2>/dev/null|head -1)
if [ -z "$LATEST" ]; then LATEST=$(ls -t fc128b/parthenon.out2.*.rhdf 2>/dev/null|head -1); fi
[ $N -ge 2 ] && [ -z "$(ls -t fc128d/parthenon.out2.*.rhdf 2>/dev/null|head -1)" ] && { echo "no fc128d restart -> STOP"; echo "no-restart" > STOP_128d; exit 0; }
[ $N -lt 40 ] && sbatch --dependency=afterany:$SLURM_JOB_ID $0 && echo "successor queued"
if [ -n "$LATEST" ]; then RA="-r $LATEST"; else RA="-i fhc_flagship.in"; fi
mkdir -p fc128d
echo "=== fc128d slot $N $(date) 128^3 curr_max_level=4 eta_ad_cap=30 restart=$LATEST ===" >> fc128d/run.log
stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN $RA -t 07:45:00 \
  parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
  parthenon/mesh/numlevel=14 refinement/curr_max_level=4 parthenon/time/tlim=1.30 parthenon/time/nlim=100000000 \
  parthenon/time/ncycle_out=20 parthenon/mesh/do_coalesced_comms=true \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 parthenon/output0/dt=-1 parthenon/output0/dn=1 \
  parthenon/output1/dn=200 parthenon/output2/dn=100 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/eta_ad_cap_code=30.0 \
  diffusion/ion_zeta=1.0e-16 -d fc128d >> fc128d/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> fc128d/run.log
