#!/bin/bash
#SBATCH --job-name=fc128gex
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
# MATCHED-EPOCH GLM CONTROL (2026-07-29). Continues the EXISTING fc128glm run from its own
# newest restart (out2.00002, t=1.102925) to tlim=1.1040, purely to answer one question:
#   does GLM, at t~1.1036, also develop the low-density magnetically-dominated region at
#   r~0.41 that the CT runs show (rho_min(r<1) ~ 1.1e-2, ~130 cells below rho=0.1, beta=2.6e-3)?
# Without this the "GLM has zero such cells" column is comparing DIFFERENT EPOCHS (fc128glm's
# only science dump is t=1.102925, and at that time the CT runs had almost no pocket either --
# fc128b had 3 cells below rho=0.1). This run makes the comparison like-for-like and decides
# whether the residual in fc128full is physical (a magnetic tower / accretion-shock cavity) or
# a third, still-unidentified defect.
# Config is UNCHANGED from the original fc128glm chain (its own knobs, no eta_ad_cap_code), so
# its trajectory stays self-consistent with its own history. tlim=1.1040 bounds the compute to
# just past the comparison epoch; slot cap 3.
# BINARY NOTE: runs athenaPK_ctfull_faf89f87, not the chain's original 1603b221. All changes
# since are CT-gated (inside CT_* / use_ct branches) and the GLM path was verified BIT-IDENTICAL
# on CPU against the pre-change log (orszag_tang_ad_glm -> 2.6565828685855680e-01). Stated here
# because it is an assumption, not a measurement, for the pre-2026-07-29 CT stencil work.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_ctfull_faf89f87
WD=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WD
WRAP=$WD/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python; RHO0=5.467e-19; STOP_CGS=1.0e-12
N=$(cat fc128glm_ext_n 2>/dev/null||echo 0); N=$((N+1)); echo $N > fc128glm_ext_n
[ -f STOP_128glmext ] && { echo "STOP -> exit"; exit 0; }
NEWEST=$(ls -t fc128glm/parthenon.out1.*.phdf 2>/dev/null|head -1)
if [ -n "$NEWEST" ]; then
  RM=$($PY - "$NEWEST" <<'PYEOF'
import sys,h5py,numpy as np
try:
 with h5py.File(sys.argv[1],"r") as h: print("%.6e"%float(np.array(h["prim"][:,0,...]).max()))
except: print("0.0")
PYEOF
)
  [ "$($PY -c "print(1 if $RM*$RHO0>=$STOP_CGS else 0)")" = "1" ] && { echo "FIRST CORE rho=$RM -> STOP"; echo "first core $RM $(date)" > STOP_128glmext; exit 0; }
fi
LATEST=$(ls -t fc128glm/parthenon.out2.*.rhdf 2>/dev/null|head -1)
[ $N -ge 2 ] && [ -z "$LATEST" ] && { echo "no restart -> STOP"; echo "no-restart" > STOP_128glmext; exit 0; }
[ $N -lt 3 ] && sbatch --dependency=afterany:$SLURM_JOB_ID $0 && echo "successor queued"
if [ -n "$LATEST" ]; then RA="-r $LATEST"; else RA="-i fhc_flagship.in"; fi
mkdir -p fc128glm
echo "=== fc128glm-EXT slot $N $(date) 128^3 curr_max_level=4 hydro/divergence_control=glm ===" >> fc128glm/run.log
stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN $RA -t 07:45:00 \
  parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
  parthenon/mesh/numlevel=14 refinement/curr_max_level=4 hydro/divergence_control=glm parthenon/time/tlim=1.1040 parthenon/time/nlim=100000000 \
  parthenon/time/ncycle_out=20 parthenon/mesh/do_coalesced_comms=true \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 parthenon/output0/dt=-1 parthenon/output0/dn=1 \
  parthenon/output1/dn=200 parthenon/output2/dn=100 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/ion_zeta=1.0e-16 -d fc128glm >> fc128glm/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> fc128glm/run.log
