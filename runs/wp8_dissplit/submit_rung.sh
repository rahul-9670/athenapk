#!/bin/bash
#SBATCH --job-name=wp8rung
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=02:40:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp8_dissplit/%x_%j.out
set -o pipefail
#
# WP-8 -- the nj8 and nj16 rungs of the dissO/dissA split. Submit with NJ=8 or NJ=16.
# DO NOT LAUNCH until the nj4 gate (submit_nj4_v2.sh) passes all three gates, in particular
# GATE A: the backported binary's OFF state must be BYTE-IDENTICAL to the pre-backport binary.
# Without that, everything measured here is contaminated by the binary change rather than by
# resolution, which is the one thing this study must not confuse.
#
# WHY A RESTART IS ENOUGH (no ladder re-run): `nonideal_eta` is {Cell, Derived, OneCopy} with no
# Restart flag, so eta lives in no output file -- but being Derived it is RECOMPUTED FROM STATE at
# init, so restarting regenerates it. Offline reconstruction is genuinely impossible here because
# diffusion/dust_coupling=true makes eta depend on the evolved grain state (f_dg, a_c), which no
# output block carries.
#
# BINARY = the WP-0 reproduction of the ladder binary 5ebddce0 (commit 29a7174 + archived patch),
# with ONLY the read-only mag_diag split columns backported (job 2495349). Ladder physics, plus
# the diagnostic. A newer binary would drag in the audit-fix chain, and A1 changes the refinement
# criterion -- block structure is precisely the variable under study.
#
# CYCLE TARGETS, derived from the ladder's own epoch_scan.txt by linear interpolation in cycle
# between the bracketing restarts, then rounded UP so each rung finishes PAST the matched epoch:
#   nj8 : restart out2.00001 (cyc 250, t=1.081307). rho=1e-12 at t~1.0898 => cycle ~463.
#         nlim=500 (t=1.091254, rho=4.66e-12) -- comfortably past. 250 cycles @ ~20 s = ~83 min.
#   nj16: restart out2.00002 (cyc 500, t=1.098885). rho=1e-12 at t~1.09949 => cycle ~607.
#         nlim=650 -- past. 150 cycles @ ~29 s = ~72 min.
# Internal -t leaves margin under the SLURM limit so the run ends cleanly and flushes its .hst.
#
# OUTPUT: out1 every 25 cycles. The .hst has NO density column, so the matched epoch has to be
# located from snapshots -- but these are only ~156 MB each here, so ~10 per leg is ~1.6 GB, which
# is nothing. out2 (restart) is suppressed: these legs must never chain, and a multi-GB rhdf is
# pure quota cost for a measurement that reads the .hst.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

: "${NJ:?set NJ (8 or 16)}"
R=/beegfs/u/bbg6470/athenapk/runs/wp8_dissplit
L=/beegfs/u/bbg6470/athenapk/runs/convergence_ladder
BIN=/beegfs/u/bbg6470/wp0_repro/athenapk/build_gpu_repro/bin/athenaPK
SPLIT=182915.67587342235      # 1e-13 g/cm3 = rho_crit, matching wp8_split_convergence.py
SHEET=0.1                     # measured: 0.3 selects 0% at nj8/nj16 (a null diagnostic)

case "$NJ" in
  # nj4's restart sits ON the matched epoch, so a handful of cycles suffices; it is included
  # here (rather than reusing the v2_on gate run) so all three rungs carry the round-3 columns
  # and come from one identical code path.
  4)  RST=$L/nj4/parthenon.out2.00001.rhdf;  NLIM=262; TLIM=00:30:00 ;;
  8)  RST=$L/nj8/parthenon.out2.00001.rhdf;  NLIM=500; TLIM=02:15:00 ;;
  16) RST=$L/nj16/parthenon.out2.00002.rhdf; NLIM=650; TLIM=02:15:00 ;;
  *)  echo "NJ must be 4, 8 or 16"; exit 2 ;;
esac

D=$R/nj$NJ; rm -rf $D; mkdir -p $D; cd $D
# wrap_slurm.sh, NOT the ladder's wrap_mod.sh. The latter derives the device from
# `nvidia-smi -L | wc -l`, which on g004 reported 6 devices for a 5-GPU allocation, so
# local_rank % NGPU put rank 2 on a card owned by another job (35,975 MiB) -> the OOM in
# job 2495385. wrap_slurm.sh indexes into the list SLURM actually exported. See its header.
cp $R/wrap_slurm.sh $D/wrap_mod.sh
echo "=== WP-8 rung nj$NJ, job $SLURM_JOB_ID $(date) host=$(hostname) ==="
echo "binary $(md5sum $BIN | cut -c1-12)   restart $RST   nlim=$NLIM"
nvidia-smi --query-gpu=index,memory.used --format=csv,noheader | sed 's/^/    pre-run /'

# ---- GPU PREFLIGHT -------------------------------------------------------------------------
# Attempt 1 (jobs 2495384/2495385) burned two slots on GPU contention, in two different ways:
#   nj8  on g002: cudaGetDeviceCount -> cudaErrorInitializationError on rank 4 (the 5th GPU).
#   nj16 on g004: Kokkos "failed to allocate 17.09 MiB" -- and its own pre-run dump showed a
#                 device holding 35,975 MiB of ANOTHER job's memory.
# `squeue -p gpu` lists no jobs while the nodes report 4/2/4/1 GPUs allocated, i.e. other users'
# work is invisible to me; contention cannot be predicted from the queue, only measured here.
# This also corrects an earlier overstatement of mine: GresEnforceBind=Yes did NOT mean no other
# job could be on a visible card -- the 36 GB reading is a direct measurement that it can.
#
# Both jobs still exited SLURM-COMPLETED 0:0 because the batch script itself succeeded, which is
# exactly why the failure has to be caught here rather than trusted to sacct.
# Scope the query to MY devices. A bare `nvidia-smi` is not reliably cgroup-filtered here: on
# g004 the batch step listed SIX devices for a five-GPU request, so a single foreign tenant on a
# card I was never given would make a whole-node check refuse a perfectly good allocation. -i
# restricts the query to the indices SLURM handed me.
NEED=5
SEL="${CUDA_VISIBLE_DEVICES:-}"
if [ -n "$SEL" ]; then Q="nvidia-smi -i $SEL"; else Q="nvidia-smi"; fi
BUSY=$($Q --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | awk '$1>1000' | wc -l)
NVIS=$($Q --query-gpu=index --format=csv,noheader 2>/dev/null | wc -l)
echo "    preflight: CUDA_VISIBLE_DEVICES='${SEL:-<unset>}' -> $NVIS GPU(s) mine, $BUSY with >1 GiB in use"
if [ "$BUSY" -gt 0 ] || [ "$NVIS" -lt "$NEED" ]; then
  echo "    ** PREFLIGHT FAIL -- refusing to start on a contended allocation."
  echo "    ** Resubmit; this is transient node contention, not a defect in the run."
  echo WP8_RUNG_PREFLIGHT_FAIL
  exit 0
fi

MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
stdbuf -oL -eL mpirun -n 5 $MCA $D/wrap_mod.sh $BIN -r $RST -t $TLIM \
  refinement/njeans=$NJ parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
  diffusion/cap_diag=true hydro/mag_diag=true \
  hydro/mag_diag_rho_split=$SPLIT hydro/mag_diag_sheet_thresh=$SHEET \
  parthenon/time/nlim=$NLIM parthenon/output0/dt=1.0e-9 \
  parthenon/output1/dn=25 parthenon/output2/dn=1000000 \
  > $D/run.log 2>&1
echo "RUN_EXIT $? $(date)  cycles=$(grep -ac '^cycle=' $D/run.log)"

echo "=== split engaged? ==="
grep -a "mag_diag .*split ON" $D/run.log | sed 's/^/    /' || echo "    ** SPLIT DID NOT ENGAGE"
echo "=== columns sane? ==="
/beegfs/u/bbg6470/venvs/analysis_env/bin/python $R/check_split.py $D/parthenon.out0.hst
echo "=== did this rung reach the matched epoch? ==="
/beegfs/u/bbg6470/venvs/analysis_env/bin/python - $D <<'PY'
import sys, glob, os, h5py, numpy as np
RHO0=5.467e-19
best=None
for f in sorted(glob.glob(os.path.join(sys.argv[1],"parthenon.out1.*.phdf"))):
    try:
        with h5py.File(f,"r") as h:
            r=float(np.array(h["prim"][:,0,...]).max()); t=float(h["Info"].attrs["Time"])
    except Exception: continue
    d=abs(np.log10(r*RHO0/1e-12))
    print(f"    {os.path.basename(f):26s} t={t:.6f} rho_max={r*RHO0:.4e} ({d:+.3f} dex)")
    if best is None or d<best[0]: best=(d,f,t,r*RHO0)
if best:
    print(f"    -> matched epoch: {os.path.basename(best[1])} t={best[2]:.6f} "
          f"rho={best[3]:.4e} ({best[0]:.3f} dex from 1e-12)")
    if best[0]>0.30: print("    ** WARNING >0.30 dex from target -- outside the ladder's acceptance window")
PY
echo WP8_RUNG_DONE
