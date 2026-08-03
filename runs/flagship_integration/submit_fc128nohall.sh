#!/bin/bash
#SBATCH --job-name=fc128nhl
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00   # diagnostic: the identical window took ~40 min in fc128noamr. Do NOT inherit the production chain 8h -- it blocks backfill scheduling.
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/%x_%j.out
set -o pipefail
# DEFECT-3 DISCRIMINATING TEST #2: HALL OFF (2026-07-29).
# Question: is the residual low-density hole caused by the Hall term under CT?
# Hall is the ONE non-ideal term fix B did not repair. CT_AddDiffusivePoynting rebuilt the
# diffusive Poynting flux from the CT edge EMF for Ohmic + ambipolar only; hall.cpp still
# deposits its own face-based cons.flux(IEN) Poynting term computed from a DIFFERENT stencil
# than the CT edge-EMF induction it pairs with -- exactly the bug class fix B addressed.
# Hall is also dispersive, unsplit-only under CT, and needs an Ohmic floor (0.05 here).
# CONTROLLED: identical to submit_fc128noamr.sh in every respect (same t=1.100 restart, same
# binary athenaPK_ctfull_faf89f87 with both fixes, FROZEN MESH via amr_check_interval=1e9,
# same tlim/cadence) EXCEPT `diffusion/hall=none`.
# PAIRED CONTROL ALREADY ON DISK: fc128noamr/parthenon.out1.00004.phdf (t=1.103633, Hall ON,
#   frozen mesh) -> rho_min(r<1)=6.8145e-02, 5 cells below 0.1, max ME/E=0.8349, 47 cells >0.5.
# READ-OUT at the matched epoch t~1.1036:
#   hole GONE (rho_min O(1), maxME/E ~0.1)  -> HALL UNDER CT IS DEFECT 3.
#   hole UNCHANGED                          -> Hall excluded; remaining suspects are the CT
#                                              interaction with self-gravity / radiation /
#                                              chemistry source terms.
# CAVEAT: turning Hall off changes the PHYSICS, not just the numerics, so a partial change is
# ambiguous; only a clean presence/absence is decisive. Frozen mesh under-resolves progressively,
# so this is valid only over the short window t=1.100 -> 1.1037. Diagnostic, not a physics run.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_ctfull_faf89f87
WD=/beegfs/u/bbg6470/athenapk/runs/flagship_integration; cd $WD
WRAP=$WD/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python; RHO0=5.467e-19; STOP_CGS=1.0e-12
N=$(cat fc128nohall_n 2>/dev/null||echo 0); N=$((N+1)); echo $N > fc128nohall_n
[ -f STOP_128nohall ] && { echo "STOP -> exit"; exit 0; }
mkdir -p fc128nohall
# first-core density check on the newest science dump
NEWEST=$(ls -t fc128nohall/parthenon.out1.*.phdf 2>/dev/null|head -1)
if [ -n "$NEWEST" ]; then
  RM=$($PY - "$NEWEST" <<'PYEOF'
import sys,h5py,numpy as np
try:
 with h5py.File(sys.argv[1],"r") as h: print("%.6e"%float(np.array(h["prim"][:,0,...]).max()))
except: print("0.0")
PYEOF
)
  [ "$($PY -c "print(1 if $RM*$RHO0>=$STOP_CGS else 0)")" = "1" ] && { echo "FIRST CORE rho=$RM -> STOP"; echo "first core $RM $(date)" > STOP_128nohall; exit 0; }
fi
# restart: first slot from the external t=1.100 diagnostic restart; later from its own
LATEST=$(ls -t fc128nohall/parthenon.out2.*.rhdf 2>/dev/null|head -1)
if [ -n "$LATEST" ]; then RA="-r $LATEST"; else RA="-r fc128b/parthenon.out2.00002.rhdf"; fi
[ $N -lt 1 ] && sbatch --dependency=afterany:$SLURM_JOB_ID $0 && echo "successor queued"
echo "=== fc128nohall slot $N $(date) binary=$(md5sum $BIN|cut -c1-8) restart=$RA ===" >> fc128nohall/run.log
stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN $RA -t 00:50:00 \
  parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
  parthenon/mesh/numlevel=14 refinement/curr_max_level=4 parthenon/time/tlim=1.1037 parthenon/time/nlim=100000000 \
  parthenon/time/ncycle_out=20 parthenon/mesh/do_coalesced_comms=true parthenon/mesh/amr_check_interval=1000000000 \
  parthenon/mesh/task_collection_timeout_in_seconds=1800 parthenon/output0/dt=-1 parthenon/output0/dn=200 \
  parthenon/output1/dn=50 parthenon/output2/dn=100 \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=400 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/eta_hall_cap_code=0.05 \
  diffusion/ion_zeta=1.0e-16 diffusion/hall=none diffusion/eta_ad_cap_code=3.0 -d fc128nohall >> fc128nohall/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> fc128nohall/run.log
