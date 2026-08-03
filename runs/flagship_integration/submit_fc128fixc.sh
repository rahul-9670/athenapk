#!/bin/bash
#SBATCH --job-name=fc128fxc
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
# CT AD FIX CONFIRMATION: forward from fc128b/out2.00002 (t=1.100) past the OLD blowup
# (t=1.1035) toward FIRST CORE, with the NEW direct-edge ambipolar CT stencil (binary
# 68497eb4). EXACT flagship config -- only the binary (AD edge-EMF) differs from the run
# that over-magnetized (ME/E->0.98) and NaN'd. Self-chaining + first-core density stop.
# New output dir fc128fix/ (fc128b restarts untouched). If it crosses t=1.1035 clean and
# reaches rho>=1e-12 cgs => CT fixed.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WD
WRAP=$WD/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python; RHO0=5.467e-19; STOP_CGS=1.0e-12
N=$(cat fc128fixc_n 2>/dev/null||echo 0); N=$((N+1)); echo $N > fc128fixc_n
[ -f STOP_128fixc ] && { echo "STOP -> exit"; exit 0; }
mkdir -p fc128fixc
# first-core density check on the newest science dump
NEWEST=$(ls -t fc128fixc/parthenon.out1.*.phdf 2>/dev/null|head -1)
if [ -n "$NEWEST" ]; then
  RM=$($PY - "$NEWEST" <<'PYEOF'
import sys,h5py,numpy as np
try:
 with h5py.File(sys.argv[1],"r") as h: print("%.6e"%float(np.array(h["prim"][:,0,...]).max()))
except: print("0.0")
PYEOF
)
  [ "$($PY -c "print(1 if $RM*$RHO0>=$STOP_CGS else 0)")" = "1" ] && { echo "FIRST CORE rho=$RM -> STOP"; echo "first core $RM $(date)" > STOP_128fixc; exit 0; }
fi
# restart: first slot from the external t=1.100 diagnostic restart; later from fc128fix's own
LATEST=$(ls -t fc128fixc/parthenon.out2.*.rhdf 2>/dev/null|head -1)
if [ -n "$LATEST" ]; then RA="-r $LATEST"; else RA="-r fc128fix/parthenon.out2.00003.rhdf"; fi
[ $N -lt 40 ] && sbatch --dependency=afterany:$SLURM_JOB_ID $0 && echo "successor queued"
echo "=== fc128fix slot $N $(date) binary=$(md5sum $BIN|cut -c1-8) restart=$RA ===" >> fc128fixc/run.log
stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN $RA -t 07:45:00 \
  parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
  parthenon/mesh/numlevel=14 refinement/curr_max_level=4 parthenon/time/tlim=1.30 parthenon/time/nlim=100000000 \
  parthenon/time/ncycle_out=20 parthenon/mesh/do_coalesced_comms=true \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 parthenon/output0/dt=-1 parthenon/output0/dn=200 \
  parthenon/output1/dn=200 parthenon/output2/dn=100 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/ion_zeta=1.0e-16 diffusion/eta_ad_cap_code=3.0 -d fc128fixc >> fc128fixc/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> fc128fixc/run.log
