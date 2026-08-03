#!/bin/bash
#SBATCH --job-name=mg_prod_tab
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/mg_prod_tab/%x_%j.out
set -o pipefail
#
# TABULATED MULTIGROUP RHD production campaign (n_group=3 + tabulated dust+gas opacity).
# Same FHC config as prod_v9 (glmMHD + self-grav + hydrogen EOS + chem + 3 non-ideal MHD +
# M1 RT) but with frequency-resolved multigroup radiation. SEPARATE dir + binary from prod_v9
# (prod_v9 uses athenaPK_eos_v9; this uses the multigroup+tabulated athenaPK). prod_v9 untouched.
#
# STOP CRITERION (user: "till first-core formation and some more time after"): the chain stops
# when the central density reaches FIRST_CORE_STOP_CGS (1e-12 g/cm^3 = one decade past the
# ~1e-13 first-core onset => first core formed AND evolved a decade). Checked each slot from the
# newest phdf before queuing the successor.
#
# SELF-CHAINING: each slot submits its successor (afterany) BEFORE running; guards stop a dead
# or finished chain within one slot. Extend the cap with: echo 0 > chain_n
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK   # multigroup + tabulated (NOT prod binary)
WDIR=/beegfs/u/bbg6470/athenapk/runs/mg_prod_tab
WRAP=$WDIR/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
MAX_CHAIN=40
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
RHO0_CGS=5.467e-19
FIRST_CORE_STOP_CGS=1.0e-12   # first-core onset ~1e-13; stop a decade past => "first core + some more"

cd $WDIR

# ---- chain bookkeeping ----
N=$(cat chain_n 2>/dev/null || echo 0); N=$((N+1)); echo $N > chain_n
echo "chain slot $N (job $SLURM_JOB_ID) $(date)"

# ---- guards: refuse to run/extend a finished or dead campaign ----
if [ -f STOP_CHAIN ]; then echo "STOP_CHAIN present -> exiting"; exit 0; fi
if grep -q "Driver completed" run.log 2>/dev/null; then
  echo "run already completed -> exiting"; exit 0; fi
if [ -f run.log ]; then
  STUCK=$(grep "^cycle=" run.log | tail -5000 | awk -F'time=' \
    'NR==1{split($2,a," ");t0=a[1]} {split($2,a," ");t1=a[1];n++} END{if(n>=5000 && t1-t0<1e-10) print "yes"}')
  if [ "$STUCK" = "yes" ]; then
    echo "FROZEN TIME over last 5000 cycles -> STOP_CHAIN"
    echo "frozen-time guard tripped, slot $N, $(date)" > STOP_CHAIN; exit 0; fi
fi
LATEST=$(ls -t $WDIR/parthenon.out2.*.rhdf 2>/dev/null | head -1)
if [ $N -ge 2 ] && [ -z "$LATEST" ]; then
  echo "slot $N but no restart exists -> earlier slot died pre-dump -> STOP_CHAIN"
  echo "no-restart guard tripped, slot $N, $(date)" > STOP_CHAIN; exit 0; fi

# ---- FIRST-CORE density stop: check newest phdf max density (code units) vs threshold ----
NEWEST_PHDF=$(ls -t $WDIR/parthenon.out1.*.phdf 2>/dev/null | head -1)
if [ -n "$NEWEST_PHDF" ]; then
  RHOMAX=$($PY - "$NEWEST_PHDF" <<'PYEOF'
import sys, h5py, numpy as np
try:
    with h5py.File(sys.argv[1],"r") as h:
        # prim component 0 = density; [block,comp,k,j,i]
        print("%.6e" % float(np.array(h["prim"][:,0,...]).max()))
except Exception as e:
    print("0.0")
PYEOF
)
  STOP=$($PY -c "rho0=$RHO0_CGS; print(1 if $RHOMAX*rho0 >= $FIRST_CORE_STOP_CGS else 0)")
  echo "first-core check: rho_max=$RHOMAX code = $($PY -c "print('%.3e'%($RHOMAX*$RHO0_CGS))") cgs (stop at $FIRST_CORE_STOP_CGS)"
  if [ "$STOP" = "1" ]; then
    echo "FIRST CORE + margin reached (rho_max*rho0 >= $FIRST_CORE_STOP_CGS) -> STOP_CHAIN"
    echo "first-core density stop, slot $N, rho_max=$RHOMAX code, $(date)" > STOP_CHAIN
    exit 0; fi
fi

# ---- submit successor BEFORE running (survives our own TIMEOUT kill) ----
if [ $N -lt $MAX_CHAIN ]; then
  sbatch --dependency=afterany:$SLURM_JOB_ID $WDIR/submit_mgtab.sh && echo "successor queued"
else
  echo "chain cap MAX_CHAIN=$MAX_CHAIN reached; extend with: echo 0 > chain_n"
fi

# ---- run ----
echo "binary:"; md5sum $BIN
if [ -n "$LATEST" ]; then RA="-r $LATEST"; echo "RESUMING FROM $LATEST"
else RA="-i fhc_mgtab.in"; echo "FRESH START FROM t=0"; fi
echo "=== slot $N start $(date) job $SLURM_JOB_ID ===" >> $WDIR/run.log
( while true; do echo "$(date +%s) $(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits | tr "\n" " ")" >> $WDIR/gpumem.log; sleep 60; done ) &
GPULOG_PID=$!
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN $RA -t 11:30:00 \
  parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 parthenon/output2/dn=250 \
  >> $WDIR/run.log 2>&1
RC=$?; kill $GPULOG_PID 2>/dev/null; echo "RUN_EXIT $RC $(date)"
tail -3 $WDIR/run.log
