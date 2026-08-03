#!/bin/bash
#SBATCH --job-name=prod_t4_full
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/%x_%j.out
set -o pipefail
# FHC production ladder t4_full — SECOND-CORE CAMPAIGN (launched 2026-07-05).
# eos=hydrogen (tabulated multi-Saha EOS) + belllin opacity + full non-ideal (AD+Ohmic+
# Hall) + chem + M1 RT, numlevel=20, binary athenaPK_eos.
#
# 2026-07-10 reconfiguration (RRZ utilisation flag + Ohmic dt-wall):
#  - 8 -> 4 GPUs: run is AMR/comm-bound, 4-GPU probe measured 1.20x better
#    GPU-seconds/zone-cycle and 36.5/79.7 GiB peak per card (wsec_AMR sink vanishes).
#  - mixed diffusion integrator (CLI overrides below + mirrored in fhc.in): parabolic
#    Ohm/AD via RKL2 STS, dispersive Hall EMF unsplit at whistler dt, Ohmic floor in the
#    RKL2 set. Lifts dt from the Ohmic-decoupling wall 4.5e-9 to the whistler/cap limit
#    (attribution: runs/dt_attrib, DEV_LOG 2026-07-10).
#  - "-t 11:30:00": clean walltime checkpoint-exit (COMPLETED, not TIMEOUT).
# 2026-07-12 speed changeover (nsys profile job 2321477 + A/B jobs, DEV_LOG 2026-07-11):
#  - binary athenaPK_eos_v2 (adds diffusion/rkl2_freeze_eta; old athenaPK_eos kept for
#    rollback; freeze-off dt regression vs old binary exact at cycle 71501).
#  - rkl2_freeze_eta=true: eta cache refreshed once per Strang half instead of every
#    RKL2 stage (PrecomputeNonidealEta was 65.7% of GPU kernel time) -> wsec_step
#    43.9 -> 20.1 s in the A/B at cycle 71500 (2.18x).
#  - rkl2_max_dt_ratio 400 -> 1000: cap was pinned every cycle at level 12; sweep at
#    cycle 71500 gave dt +59% for wsec +21% (net +31%); whistler binds near 1000 so
#    higher caps buy nothing. RKL2 eigenmode-validated at 1000x dt (0.68%).
#  - combined config probed on restart final.rhdf (cycle 71913) before this went live.
# 2026-07-12 night changeover (user-approved set, DEV_LOG 2026-07-12 late):
#  - binary athenaPK_eos_v4: fused ionization_chem evaluator (one Wardle-tensor solve
#    per cell instead of three; +7.4% in 3-way probe 2334442) + eta_ohm_cap_code param.
#    v2/v3 kept for rollback.
#  - diffusion/eta_ohm_cap_code=0.1 (PHYSICS, user signed off): caps ionization eta_O
#    in the diffusion-decoupled first-core interior. Priced by probe 2344652 at cycle
#    73000: dt 2.5e-8 -> 6.8e-8 (whistler/hydro-limited, STS ratio unpins 1000 -> ~50,
#    45 -> 11 stages), wsec_step 34.1 -> 15.1 s => 6.1x sim-time/wall vs no-cap.
#  - parthenon/output2/dn=250: restart cadence ~halved (dt is now 4-5x larger so wall
#    time between restarts shrank anyway; belt-and-braces for AMR-level jumps).
#  - 4 -> 5 GPUs (single node; g00x carry 8xH100): headroom for level 13+ block growth
#    (was 70-72 GiB/card = 87-88% on 4). NOTE a 5-GPU slot can queue longer than a
#    4-GPU one on a packed cluster; chain is afterany so worst case is wait, not loss.
#
# 2026-07-16 relaunch staging (chain stopped 2026-07-14 15:34 at ckpt 00163/cycle 77250,
# user-directed, for Hall retooling; DEV_LOG 2026-07-15):
#  - binary athenaPK_eos_v6 (md5 03a299cf): v5's diffusion/eta_hall_cap_code feature
#    (inactive without the input key) + the RKL2 STS named-pack fix — the {Independent}
#    packs in AddSTSTasks were clobbering rad.Er/Fr (outputs dead-zero since cycle
#    71000) and grav.phi (zero initial guess each cycle). Gas evolution bit-identical;
#    rad outputs restart-continuous post-fix. v4 = rollback.
#  - eta_hall_cap_code DELIBERATELY NOT SET: hallcap probe 2359338 (2026-07-17, pinned
#    to quarantined ckpt 00163) priced cap 0.3 = dt bit-identical (never binds) and
#    cap 0.1 = -2% (noise); binding-cell |eta_H| ~ 7e-3 is far below any safe cap.
#    Off-state check PASS (v5 no-cap bit-identical to v4). Feature stays in v6, inert.
#  - RELAUNCH 2026-07-17: resumes from 00142 (cycle 71000, last clean pre-bug state,
#    healthy Er) — post-bug outputs quarantined in quarantine_postbug_71000/. The redo
#    runs with correct thermodynamics (hotter core -> thermal K ionization -> lower
#    eta_H/eta_O expected; whistler wall may recede organically).
#
# SELF-CHAINING: each slot submits its successor (afterany dependency) BEFORE running,
# and every slot starts with guards so a dead run stops the chain within one slot:
#   - STOP_CHAIN file present            -> exit (touch STOP_CHAIN to stop manually)
#   - "Driver completed" in run.log      -> tlim reached, done, exit
#   - frozen time (t1-style corruption: <1e-10 code time over last 5000 cycles) -> STOP
#   - slot >=2 but no restart file       -> config broken, STOP
# Chain capped at MAX_CHAIN slots; extend with: echo 0 > chain_n
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_eos_v6
WDIR=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full
WRAP=$WDIR/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
MAX_CHAIN=30
cd $WDIR

# ---- chain bookkeeping ----
N=$(cat chain_n 2>/dev/null || echo 0); N=$((N+1)); echo $N > chain_n
echo "chain slot $N (job $SLURM_JOB_ID) $(date)"

# ---- guards: refuse to run/extend a finished or dead campaign ----
if [ -f STOP_CHAIN ]; then echo "STOP_CHAIN present -> exiting"; exit 0; fi
if grep -q "Driver completed" run.log 2>/dev/null; then
  echo "run already completed (tlim) -> exiting"; exit 0; fi
if [ -f run.log ]; then
  STUCK=$(grep "^cycle=" run.log | tail -5000 | awk -F'time=' \
    'NR==1{split($2,a," ");t0=a[1]} {split($2,a," ");t1=a[1];n++} END{if(n>=5000 && t1-t0<1e-10) print "yes"}')
  if [ "$STUCK" = "yes" ]; then
    echo "FROZEN TIME over last 5000 cycles (t1-style dead state) -> STOP_CHAIN"
    echo "frozen-time guard tripped, slot $N, $(date)" > STOP_CHAIN; exit 0; fi
fi
LATEST=$(ls -t $WDIR/parthenon.out2.*.rhdf 2>/dev/null | head -1)
if [ $N -ge 2 ] && [ -z "$LATEST" ]; then
  echo "slot $N but no restart exists -> earlier slot died before first dump -> STOP_CHAIN"
  echo "no-restart guard tripped, slot $N, $(date)" > STOP_CHAIN; exit 0; fi

# ---- submit successor BEFORE running (survives our own TIMEOUT kill) ----
if [ $N -lt $MAX_CHAIN ]; then
  sbatch --dependency=afterany:$SLURM_JOB_ID $WDIR/submit.sh && echo "successor queued"
else
  echo "chain cap MAX_CHAIN=$MAX_CHAIN reached; no successor (echo 0 > chain_n to extend)"
fi

# ---- run ----
echo "binary:"; md5sum $BIN
if [ -n "$LATEST" ]; then RA="-r $LATEST"; echo "RESUMING FROM $LATEST"
else RA="-i fhc.in"; echo "FRESH START FROM t=0"; fi
echo "=== slot $N start $(date) job $SLURM_JOB_ID ===" >> $WDIR/run.log
# GPU memory logger: t1 died of allocator growth under AMR churn; keep the evidence.
( while true; do echo "$(date +%s) $(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits | tr "\n" " ")" >> $WDIR/gpumem.log; sleep 60; done ) &
GPULOG_PID=$!
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN $RA -t 11:30:00 \
  parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 parthenon/output2/dn=250 \
  >> $WDIR/run.log 2>&1
RC=$?; kill $GPULOG_PID 2>/dev/null; echo "RUN_EXIT $RC $(date)"
tail -3 $WDIR/run.log
