#!/bin/bash
#SBATCH --job-name=flag_fc
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --gres=gpu:h100:6
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
# Full-flagship config -> FIRST CORE. Physics verified genuine (coherent collapsing core, r^-2
# profile, infall, barotropic matches mg_prod_tab). The t=1.05 OOM was a real memory limit of the
# CT + multigroup + current-sheet-AMR + coalesced footprint on 5 GPUs at deep refinement.
# Fix = (1) 6 GPUs (more HBM than the 5-GPU OOM; a genuine fixed shortfall, unlike the falsified
# envelope case; 6 not 8 so it fits partially-freed nodes and schedules sooner) +
# (2) numlevel 20->14 (first core needs only ~5 levels; caps the current-sheet cascade that was
# refining thin sheets to many levels). Self-chaining + restart so it survives the walltime and
# runs through first core. All fixes carried (binary 19be8c1b). prod_v9 untouched.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WDIR=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WDIR
WRAP=$WDIR/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python; RHO0=5.467e-19; STOP_CGS=1.0e-12

N=$(cat fc_chain_n 2>/dev/null || echo 0); N=$((N+1)); echo $N > fc_chain_n
[ -f STOP_FC ] && { echo "STOP_FC -> exit"; exit 0; }
# first-core density stop
NEWEST=$(ls -t fc/parthenon.out1.*.phdf 2>/dev/null | head -1)
if [ -n "$NEWEST" ]; then
  RM=$($PY - "$NEWEST" <<'PYEOF'
import sys,h5py,numpy as np
try:
 with h5py.File(sys.argv[1],"r") as h: print("%.6e"%float(np.array(h["prim"][:,0,...]).max()))
except: print("0.0")
PYEOF
)
  [ "$($PY -c "print(1 if $RM*$RHO0>=$STOP_CGS else 0)")" = "1" ] && {
    echo "FIRST CORE reached rho=$RM*rho0>=$STOP_CGS -> STOP"; echo "first core $RM $(date)" > STOP_FC; exit 0; }
fi
LATEST=$(ls -t fc/parthenon.out2.*.rhdf 2>/dev/null | head -1)
if [ $N -ge 2 ] && [ -z "$LATEST" ]; then echo "no restart -> STOP"; echo "no-restart $(date)" > STOP_FC; exit 0; fi
[ $N -lt 30 ] && sbatch --dependency=afterany:$SLURM_JOB_ID $0 && echo "successor queued"

if [ -n "$LATEST" ]; then RA="-r $LATEST"; else RA="-i fhc_flagship.in"; fi
mkdir -p fc
echo "=== flag_fc slot $N $(date) numlevel=14 6GPU ===" >> fc/run.log
stdbuf -oL -eL mpirun -n 6 $MCA $WRAP $BIN $RA -t 03:45:00 \
  parthenon/mesh/numlevel=14 parthenon/time/tlim=1.10 parthenon/time/nlim=1000000 \
  parthenon/time/ncycle_out=5 parthenon/mesh/do_coalesced_comms=true \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 \
  parthenon/output1/dn=100 parthenon/output2/dn=300 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/ion_zeta=1.0e-16 -d fc >> fc/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> fc/run.log
