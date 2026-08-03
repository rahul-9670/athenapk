#!/bin/bash
#SBATCH --job-name=fc128dia
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00   # diagnostic: the identical window took ~40 min in fc128noamr. Do NOT inherit the production chain 8h -- it blocks backfill scheduling.
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
# DEFECT-3 DIRECT MEASUREMENT (2026-07-29). Seven ablations located the defect but produced
# no mechanism; this measures it instead of eliminating around it.
# hydro/ct_proj_diag=true records per cell
#     dEint_rel = (eint_post - eint_pre)/eint_pre = -(me_post - me_pre)/eint_pre
# i.e. the RELATIVE internal-energy transfer the face->cell projection imposes, taken BEFORE
# the eint guard acts, so it is the raw E-vs-ME bookkeeping mismatch. Negative = the projection
# is COOLING the gas (the damaging direction). Exposed as history vars ct_projEintMin /
# ct_projEintMaxAbs (per cycle) and as the derived field ct.dEint (per cell, for localization).
# Binary athenaPK_ctdiag_42f8c311. CPU-validated: OFF-state bit-identical; signal separates the
# known-broken from the known-fixed CT+AD+RKL2 config by 3.9x at t=0.1 as the shock strengthens.
# Frozen mesh (amr_check_interval=1e9) + Hall ON = exactly the fc128noamr configuration, whose
# hole is characterized on disk (t=1.103633: rho_min 6.81e-2, 5 cells <0.1, max ME/E 0.835).
# READ-OUT: compare ct_projEintMin against the hole's ME/IE ~ 150.
#   |dEint_rel| small (<<1/150) where the hole forms -> projection EXONERATED; the cause is in
#       the source-term coupling (self-gravity / radiation / chemistry) and this line closes.
#   |dEint_rel| O(1/150) or larger there              -> MECHANISM IDENTIFIED, with a number.
# Localize with the ct.dEint field in the phdf against the rho<0.1 cells. NOTE: the deck's
# output1 'variables' list is explicit, so ct.dEint MUST be added there or the field is
# computed but never written; overridden on the CLI below (rad.* dropped -- not needed here
# and it halves the dump size).
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_ctdiag_42f8c311
WD=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WD
WRAP=$WD/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python; RHO0=5.467e-19; STOP_CGS=1.0e-12
N=$(cat fc128diag_n 2>/dev/null||echo 0); N=$((N+1)); echo $N > fc128diag_n
[ -f STOP_128diag ] && { echo "STOP -> exit"; exit 0; }
mkdir -p fc128diag
# first-core density check on the newest science dump
NEWEST=$(ls -t fc128diag/parthenon.out1.*.phdf 2>/dev/null|head -1)
if [ -n "$NEWEST" ]; then
  RM=$($PY - "$NEWEST" <<'PYEOF'
import sys,h5py,numpy as np
try:
 with h5py.File(sys.argv[1],"r") as h: print("%.6e"%float(np.array(h["prim"][:,0,...]).max()))
except: print("0.0")
PYEOF
)
  [ "$($PY -c "print(1 if $RM*$RHO0>=$STOP_CGS else 0)")" = "1" ] && { echo "FIRST CORE rho=$RM -> STOP"; echo "first core $RM $(date)" > STOP_128diag; exit 0; }
fi
# restart: first slot from the external t=1.100 diagnostic restart; later from its own
LATEST=$(ls -t fc128diag/parthenon.out2.*.rhdf 2>/dev/null|head -1)
if [ -n "$LATEST" ]; then RA="-r $LATEST"; else RA="-r fc128b/parthenon.out2.00002.rhdf"; fi
[ $N -lt 1 ] && sbatch --dependency=afterany:$SLURM_JOB_ID $0 && echo "successor queued"
echo "=== fc128diag slot $N $(date) binary=$(md5sum $BIN|cut -c1-8) restart=$RA ===" >> fc128diag/run.log
stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN $RA -t 00:50:00 \
  parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
  parthenon/mesh/numlevel=14 refinement/curr_max_level=4 parthenon/time/tlim=1.1037 parthenon/time/nlim=100000000 \
  parthenon/time/ncycle_out=20 parthenon/mesh/do_coalesced_comms=true parthenon/mesh/amr_check_interval=1000000000 \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 parthenon/output0/dt=-1 parthenon/output0/dn=200 \
  parthenon/output1/dn=50 hydro/ct_proj_diag=true 'parthenon/output1/variables=prim, grav.phi, ct.dEint' parthenon/output2/dn=100 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/ion_zeta=1.0e-16 diffusion/eta_ad_cap_code=3.0 -d fc128diag >> fc128diag/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> fc128diag/run.log
