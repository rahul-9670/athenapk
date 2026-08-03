#!/bin/bash
#SBATCH --job-name=conv
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
# One resolution-ladder level for the flagship convergence study. Parameterized by env:
#   NJEANS  = Jeans resolution (4 / 8 / 16)
#   RUNDIR  = this level's run directory (e.g. .../convergence_ladder/nj8)
#   FIRST_CORE_STOP_CGS = matched stop epoch (default 1e-13 = first-core onset)
# Same IC/physics as mg_prod_tab (tabulated multigroup), only refinement/njeans differs.
# Self-chaining + timeout-safe + density stop, mirroring mg_prod_tab/submit_mgtab.sh.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

: "${NJEANS:?set NJEANS}"; : "${RUNDIR:?set RUNDIR}"
FIRST_CORE_STOP_CGS="${FIRST_CORE_STOP_CGS:-1.0e-13}"
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK   # multigroup+tabulated (NOT prod binary)
# 2026-07-30: the ladder must converge the FLAGSHIP configuration, not the older mg_prod_tab one.
# fhc_ladder.in is a byte-copy of runs/prod_flagship_test/fhc_flagship.in (md5 65451669), i.e.
# mg_prod_tab + WS-4 dust (nscalars=7, <physics> dust=true, <dust> block) + the audit items
# cap_diag / dust_coupling / mag_diag. Kept as a LOCAL copy so a later edit to the test deck
# cannot silently change what a running ladder is converging.
DECK=/beegfs/u/bbg6470/athenapk/runs/convergence_ladder/fhc_ladder.in
WRAP=$RUNDIR/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
MAX_CHAIN=40; PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python; RHO0=5.467e-19
cd $RUNDIR

N=$(cat chain_n 2>/dev/null || echo 0); N=$((N+1)); echo $N > chain_n
echo "conv nj=$NJEANS slot $N (job $SLURM_JOB_ID) $(date)"
[ -f STOP_CHAIN ] && { echo "STOP_CHAIN -> exit"; exit 0; }
grep -q "Driver completed" run.log 2>/dev/null && { echo "completed -> exit"; exit 0; }
LATEST=$(ls -t $RUNDIR/parthenon.out2.*.rhdf 2>/dev/null | head -1)
if [ $N -ge 2 ] && [ -z "$LATEST" ]; then echo "no restart -> STOP"; echo "no-restart $(date)" > STOP_CHAIN; exit 0; fi

# density stop (matched epoch)
NEWEST=$(ls -t $RUNDIR/parthenon.out1.*.phdf 2>/dev/null | head -1)
if [ -n "$NEWEST" ]; then
  RHOMAX=$($PY - "$NEWEST" <<'PYEOF'
import sys,h5py,numpy as np
try:
    with h5py.File(sys.argv[1],"r") as h: print("%.6e"%float(np.array(h["prim"][:,0,...]).max()))
except: print("0.0")
PYEOF
)
  if [ "$($PY -c "print(1 if $RHOMAX*$RHO0>=$FIRST_CORE_STOP_CGS else 0)")" = "1" ]; then
    echo "matched epoch reached ($RHOMAX*rho0>=$FIRST_CORE_STOP_CGS) -> STOP"
    echo "epoch stop nj=$NJEANS rho_max=$RHOMAX $(date)" > STOP_CHAIN; exit 0; fi
fi

[ $N -lt $MAX_CHAIN ] && sbatch --dependency=afterany:$SLURM_JOB_ID --job-name=conv_nj$NJEANS \
  --export=ALL,NJEANS=$NJEANS,RUNDIR=$RUNDIR,FIRST_CORE_STOP_CGS=$FIRST_CORE_STOP_CGS $0 && echo "successor queued"

echo "binary:"; md5sum $BIN
if [ -n "$LATEST" ]; then RA="-r $LATEST"; else RA="-i $DECK"; fi
echo "=== nj=$NJEANS slot $N start $(date) ===" >> $RUNDIR/run.log
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN $RA -t 11:30:00 \
  refinement/njeans=$NJEANS parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 parthenon/output2/dn=250 \
  diffusion/cap_diag=true hydro/mag_diag=true \
  >> $RUNDIR/run.log 2>&1
echo "RUN_EXIT $? $(date)"