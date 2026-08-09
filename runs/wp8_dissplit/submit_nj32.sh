#!/bin/bash
#SBATCH --job-name=wp8nj32
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp8_dissplit/%x_%j.out
set -o pipefail
#
# WP-8 FOURTH RUNG: njeans = 32. Self-chaining, density-stopped.
#
# WHY. With three rungs, dissA-smooth and Jsq-smooth are NON-MONOTONE (+77.6 %/-59.9 % and
# +139.2 %/-28.7 %), and nj4 has no business in a convergence fit at all -- 84.5 % of its current
# and 90.8 % of its ambipolar heating sit on grid-scale current, i.e. it is simply under-resolved.
# Dropping nj4 leaves TWO points, which cannot establish convergence. nj32 gives three clean ones
# (nj8, nj16, nj32) and is the only way to separate "not yet asymptotic" from "genuinely
# divergent" for the two quantities WP-8 still has open.
#
# WHY A FULL LEG AND NOT A RESTART. The other three rungs were restarted from the ladder's own
# checkpoints, which cost ~11 GPU-h total. There is no njeans=32 checkpoint, and one cannot be
# faked: restarting nj16 at higher njeans would carry nj16's refinement history through the whole
# earlier collapse, so it would not be an njeans=32 run. A ladder rung is defined by its
# refinement criterion applied from t=0.
#
# SIZING, measured rather than assumed:
#   * MEMORY. nj16 used 43.7 GiB/GPU with 1373 blocks on 5 ranks = 275 blocks/rank, i.e.
#     0.159 GiB/block/rank -- exactly the D1 law, so the law does apply to this deck. Capacity is
#     ~500 blocks/rank on an 80 GiB H100. Block count barely grows with njeans (nj8 -> nj16 at the
#     matched epoch: 1317 -> 1373 blocks, +4 %) because a finer Jeans criterion DEEPENS the
#     refinement rather than widening it. nj32 should land ~1400-1600 blocks = 45-51 GiB. Fits.
#   * TIME is the real cost: finer cells => smaller dt. nj16 reached rho = 1e-12 at cycle ~607
#     from t=0; nj32 should need roughly 2-3x that. Estimated 50-80 GPU-h, hence chaining.
#
# STOP AT 2e-12, not at the ladder's default. This rung exists only to be read at the matched
# epoch rho = 1e-12; 2e-12 clears it with margin and stops the run there instead of integrating
# into the second core, which is where the other legs spent most of their 12 h slots.
#
# 5 RANKS on purpose -- the same decomposition as the other three rungs.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

R=/beegfs/u/bbg6470/athenapk/runs/wp8_dissplit
L=/beegfs/u/bbg6470/athenapk/runs/convergence_ladder
D=$R/nj32
BIN=/beegfs/u/bbg6470/wp0_repro/athenapk/build_gpu_repro/bin/athenaPK   # 9a9ee292, round-3
DECK=$L/fhc_ladder.in
SPLIT=182915.67587342235
SHEET=0.1
STOP_CGS=2.0e-12
MAX_CHAIN=25
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
RHO0=5.467e-19

mkdir -p $D; cd $D
cp -f $R/wrap_slurm.sh $D/wrap_mod.sh          # SLURM-list pinning; see its header
N=$(cat chain_n 2>/dev/null || echo 0); N=$((N+1)); echo $N > chain_n
echo "=== WP-8 nj32 slot $N (job $SLURM_JOB_ID) $(date) host=$(hostname) ==="
echo "binary $(md5sum $BIN | cut -c1-12)   stop at rho >= $STOP_CGS g/cm3"

[ -f STOP_CHAIN ] && { echo "STOP_CHAIN present -> exit"; cat STOP_CHAIN; exit 0; }
grep -aq "Driver completed" run.log 2>/dev/null && { echo "already completed -> exit"; exit 0; }

# --- GPU preflight (see submit_rung.sh for the two failures that motivated it) ---
SEL="${CUDA_VISIBLE_DEVICES:-}"
if [ -n "$SEL" ]; then Q="nvidia-smi -i $SEL"; else Q="nvidia-smi"; fi
BUSY=$($Q --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | awk '$1>1000' | wc -l)
NVIS=$($Q --query-gpu=index --format=csv,noheader 2>/dev/null | wc -l)
echo "    preflight: $NVIS GPU(s) mine, $BUSY with >1 GiB in use"
if [ "$BUSY" -gt 0 ] || [ "$NVIS" -lt 5 ]; then
  echo "    ** PREFLIGHT FAIL -- contended allocation; requeueing this slot without consuming it."
  [ $N -lt $MAX_CHAIN ] && sbatch --exclude=g004 --job-name=wp8nj32 $0 && echo "    requeued"
  echo WP8_NJ32_PREFLIGHT_FAIL
  exit 0
fi

# --- density stop: has the PREVIOUS slot already passed the epoch? ---
NEWEST=$(ls -t $D/parthenon.out1.*.phdf 2>/dev/null | head -1)
if [ -n "$NEWEST" ]; then
  RHOMAX=$($PY - "$NEWEST" <<'PYEOF'
import sys, h5py, numpy as np
try:
    with h5py.File(sys.argv[1], "r") as h:
        print("%.6e" % float(np.array(h["prim"][:, 0, ...]).max()))
except Exception:
    print("0.0")
PYEOF
)
  echo "    newest snapshot rho_max = $RHOMAX code = $($PY -c "print('%.4e'%($RHOMAX*$RHO0))") g/cm3"
  if [ "$($PY -c "print(1 if $RHOMAX*$RHO0>=$STOP_CGS else 0)")" = "1" ]; then
    echo "STOP epoch reached -> done"; echo "epoch stop $(date) rho=$RHOMAX" > STOP_CHAIN
    echo WP8_NJ32_REACHED_EPOCH; exit 0; fi
fi

# queue the successor BEFORE running, so a timeout still chains (afterany)
[ $N -lt $MAX_CHAIN ] && sbatch --dependency=afterany:$SLURM_JOB_ID --exclude=g004 \
  --job-name=wp8nj32 $0 >/dev/null && echo "    successor queued"

LATEST=$(ls -t $D/parthenon.out2.*.rhdf 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then RA="-r $LATEST"; echo "    resuming from $(basename $LATEST)";
else RA="-i $DECK"; echo "    fresh start from t=0"; fi

MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
echo "=== slot $N start $(date) ===" >> $D/run.log
stdbuf -oL -eL mpirun -n 5 $MCA $D/wrap_mod.sh $BIN $RA -t 11:30:00 \
  refinement/njeans=32 parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
  diffusion/cap_diag=true hydro/mag_diag=true \
  hydro/mag_diag_rho_split=$SPLIT hydro/mag_diag_sheet_thresh=$SHEET \
  parthenon/output0/dt=1.0e-9 parthenon/output1/dn=25 parthenon/output2/dn=250 \
  >> $D/run.log 2>&1
echo "RUN_EXIT $? $(date)  cycles=$(grep -ac '^cycle=' $D/run.log)"
grep -a "^cycle=" $D/run.log | tail -1 | sed 's/ zone.*//;s/^/    /'
echo WP8_NJ32_SLOT_DONE
