#!/bin/bash
#SBATCH --job-name=ens
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:h100:4
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
set -o pipefail
#
# One ensemble member (flagship Phase 7 IC ensemble). Parameterized by env:
#   PDIR = this point's dir (holds fhc_ens.in from design.py); STOP_CGS = matched epoch (1e-13).
# Same tabulated-multigroup binary + physics as mg_prod_tab; only the sampled IC keys differ.
# Self-chaining + timeout-safe + density stop, mirroring mg_prod_tab/submit_mgtab.sh.
#
# SIZING, 2026-08-06 -- 5 GPUs/12 h was the WORST available choice; measured, not guessed.
#
# PACKING. The gpu partition is 4 nodes holding 8/8/8/7 H100 = 31 total. At 5 GPUs a job cannot
# be packed twice onto an 8-GPU node (8//5 = 1), so 5-GPU jobs give 4 concurrent members and
# leave 11 of 31 GPUs IDLE (35.5 % wasted). At 4 GPUs exactly two members fit per 8-GPU node:
# 7 concurrent, 3 idle (9.7 %). That is ~75 % more members in flight for the same hardware.
# (2 GPUs packs best of all -- 15 concurrent, 3.2 % waste -- but see MEMORY.)
#
# MEMORY, from runs/mg_prod_tab/gpumem.log (the same physics, 5 GPUs, 144 samples). Consumption
# is close to linear in blocks/rank: 28 GiB at 172 blk/rank, 42 GiB at 275 => ~0.16 GiB/block.
# This ensemble stops at rho=1e-13, which mg_prod_tab's own rho(t) table places near ~900 blocks
# (t=1.08 -> 862 blk), NOT its 1373-block endpoint at 1e-12. So per card:
#     4 GPUs -> ~225 blk/rank -> ~36 GiB of 79.2  (45 %, comfortable)
#     3 GPUs -> ~300 blk/rank -> ~48 GiB          (61 %, ok but packs worse: 8 concurrent)
#     2 GPUs -> ~450 blk/rank -> ~72 GiB          (91 % -- REJECTED, no headroom for AMR
#                                                  overshoot, and D1 shows growth is superlinear
#                                                  in depth)
#
# TIME. mg_prod_tab ran t=0 to rho=7.25e-12 -- well PAST this ensemble's 1e-13 stop -- in
# 2.39 h wall on 5 GPUs. A 12 h slot was ~5x the whole job. 4 h at 4 GPUs keeps ~1.3x margin on
# the measured time even after the 5->4 GPU slowdown, and short jobs backfill enormously better
# on a 4-node partition. Shortening is SAFE because this script self-chains: a member that needs
# longer simply takes another slot (MAX_CHAIN=30). The only cost of being wrong is one extra
# slot; the cost of 12 h was permanent queue starvation.
#
# 2026-08-06, SECOND SHRINK: 4 h -> 2 h. mg_prod_tab reached rho=7.25e-12 -- past this ensemble's
# 1e-13 stop -- in 2.39 h on 5 GPUs, so 2 h at 4 GPUs will usually need a second slot. That is the
# intended trade: chaining makes truncation free apart from one restart, and short jobs backfill
# far better on a 4-node partition. A shorter slot ALSO caps how far a member can overshoot the
# density stop (which is only checked BETWEEN slots), which lowers peak block count and therefore
# peak memory.
#
# 4 GPUs IS THE FLOOR -- do not shrink further. With the memory law now measured
# (D1: consumed = 0.159 GiB x blocks/rank, +/-0.8%, docs/validation/D1_gpu_memory_imbalance.md):
#   overshoot to mg_prod_tab's 1373-block endpoint  =>  4 GPUs: 343 blk/rank = 54.6 GiB (69%, safe)
#                                                       3 GPUs: 458 blk/rank = 72.8 GiB (92%, NO)
# and on throughput 3 GPUs is worse anyway: the partition (8/8/8/7 = 31 H100) fits 8 concurrent
# 3-GPU jobs = 24 GPUs busy, versus 7 concurrent 4-GPU jobs = 28 busy. Fewer GPUs per job here
# means LESS aggregate compute, not more.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

: "${PDIR:?set PDIR}"; STOP_CGS="${STOP_CGS:-1.0e-13}"
# Science-output cadence, in CYCLES. Default 50 = the value in every fhc_ens.in, so leaving OUT_DN
# unset reproduces the original campaign exactly. It MUST be passed on the command line rather than
# edited into the deck: a restarted run takes its parameters from the rhdf + CLI, not from the .in.
#
# Measured cadence behaviour (from the completed campaign, docs in README):
#   * BEFORE first core the collapse runs away and dn=50 spans 1.3-1.9 decades of rho per snapshot;
#   * AFTER first core dt collapses and dn=50 spans only 0.08-0.17 decades (point000, 00007->final).
# So coarse output only loses resolution across the first-core transition itself. That is what cost
# point021/point023 their 1e-12 measurement -- both had ample density margin (31.3x, 26.4x) but
# their cadence jumped straight over the acceptance window.
OUT_DN="${OUT_DN:-50}"
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
# matched-epoch density stop.
#
# FIXED 2026-08-06: this used to read ONLY the newest phdf and fall back to `print("0.0")` on
# any exception. A truncated newest snapshot -- exactly what a mid-write kill produces, and
# slots here end on a timer -- therefore read as rho=0, so the member RESTARTED instead of
# stopping, silently burning another whole slot past its own finish line. It now walks BACK to
# the newest READABLE snapshot and reports which file it used, so a corrupt tail costs one
# snapshot of resolution instead of the entire stop condition. It also returns -1 rather than 0
# when nothing is readable, so "no data yet" is distinguishable from "density is zero".
RM=$($PY - "$PDIR" <<'PYEOF'
import sys, glob, os, h5py, numpy as np
fs = sorted(glob.glob(os.path.join(sys.argv[1], "parthenon.out1.*.phdf")),
            key=os.path.getmtime, reverse=True)
for f in fs:
    try:
        with h5py.File(f, "r") as h:
            print("%.6e %s" % (float(np.array(h["prim"][:, 0, ...]).max()),
                               os.path.basename(f)))
        break
    except Exception:
        continue          # truncated / half-written -> fall back to the previous snapshot
else:
    print("-1.0 none")    # NO readable snapshot at all -- NOT the same as density zero
PYEOF
)
RMV=${RM%% *}; RMF=${RM##* }
if [ "$RMV" = "-1.0" ]; then
  echo "no readable snapshot yet (fresh start, or all truncated) -> continuing"
elif [ "$($PY -c "print(1 if $RMV*$RHO0>=$STOP_CGS else 0)")" = "1" ]; then
  echo "matched epoch reached (rho=$RMV code, from $RMF) -> STOP"
  echo "epoch stop rho=$RMV from $RMF $(date)" > STOP_CHAIN; exit 0
fi
[ $N -lt $MAX_CHAIN ] && sbatch --dependency=afterany:$SLURM_JOB_ID --export=ALL,PDIR=$PDIR,STOP_CGS=$STOP_CGS,OUT_DN=$OUT_DN $0 && echo "successor queued"

if [ -n "$LATEST" ]; then RA="-r $LATEST"; else RA="-i fhc_ens.in"; fi
echo "=== $(basename $PDIR) slot $N $(date) ===" >> $PDIR/run.log
stdbuf -oL -eL mpirun -n 4 $MCA $WRAP $BIN $RA -t 01:45:00 \
  parthenon/mesh/do_coalesced_comms=true diffusion/integrator=rkl2 \
  diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=1000 \
  diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 parthenon/output2/dn=250 \
  parthenon/output1/dn=$OUT_DN \
  >> $PDIR/run.log 2>&1
echo "RUN_EXIT $? $(date)" >> $PDIR/run.log
