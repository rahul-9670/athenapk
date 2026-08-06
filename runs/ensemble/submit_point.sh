#!/bin/bash
#SBATCH --job-name=ens
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
set -o pipefail
#
# One ensemble member (flagship Phase 7 IC ensemble). Parameterized by env:
#   PDIR = this point's dir (holds fhc_ens.in from design.py); STOP_CGS = matched epoch (1e-13).
# Same tabulated-multigroup binary + physics as mg_prod_tab; only the sampled IC keys differ.
# Self-chaining + timeout-safe + density stop, mirroring mg_prod_tab/submit_mgtab.sh.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

: "${PDIR:?set PDIR}"; STOP_CGS="${STOP_CGS:-1.0e-13}"
# PINNED to the preserved hard link, NOT build_gpu/bin/athenaPK. That path is the scratch build
# slot: `make -C build_gpu` overwrites it, and with MAX_CHAIN=30 across 24 members this campaign
# can span up to 720 chained jobs over days. A rebuild mid-campaign would silently change the
# binary underneath half the ensemble, and the members already finished would not be comparable
# with the ones still running -- with nothing in the output to show it happened.
# 84a6d248 = 0d3a559 (multigroup + tabulated + rtsafe + the 2026-08-05/06 audit batch, incl. the
# N3 restart fix). Provenance frozen in docs/provenance/binary_84a6d248/.
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_PRESERVED_84a6d248
WRAP=$PDIR/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
MAX_CHAIN=30; PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python; RHO0=5.467e-19
cd $PDIR

N=$(cat chain_n 2>/dev/null || echo 0); N=$((N+1)); echo $N > chain_n
[ -f STOP_CHAIN ] && { echo "STOP_CHAIN -> exit"; exit 0; }
grep -q "Driver completed" run.log 2>/dev/null && { echo "completed -> exit"; exit 0; }
LATEST=$(ls -t $PDIR/parthenon.out2.*.rhdf 2>/dev/null | head -1)
if [ $N -ge 2 ] && [ -z "$LATEST" ]; then echo "no restart -> STOP"; echo "no-restart $(date)" > STOP_CHAIN; exit 0; fi
# matched-epoch density stop
NEWEST=$(ls -t $PDIR/parthenon.out1.*.phdf 2>/dev/null | head -1)
if [ -n "$NEWEST" ]; then
  RM=$($PY - "$NEWEST" <<'PYEOF'
import sys,h5py,numpy as np
try:
 with h5py.File(sys.argv[1],"r") as h: print("%.6e"%float(np.array(h["prim"][:,0,...]).max()))
except: print("0.0")
PYEOF
)
  [ "$($PY -c "print(1 if $RM*$RHO0>=$STOP_CGS else 0)")" = "1" ] && {
    echo "matched epoch reached -> STOP"; echo "epoch stop rho=$RM $(date)" > STOP_CHAIN; exit 0; }
fi
[ $N -lt $MAX_CHAIN ] && sbatch --dependency=afterany:$SLURM_JOB_ID --export=ALL,PDIR=$PDIR,STOP_CGS=$STOP_CGS $0 && echo "successor queued"

if [ -n "$LATEST" ]; then RA="-r $LATEST"; else RA="-i fhc_ens.in"; fi
echo "=== $(basename $PDIR) slot $N $(date) ===" >> $PDIR/run.log
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN $RA -t 11:30:00 \
  parthenon/mesh/do_coalesced_comms=true diffusion/integrator=rkl2 \
  diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=1000 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 parthenon/output2/dn=250 \
  >> $PDIR/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> $PDIR/run.log
