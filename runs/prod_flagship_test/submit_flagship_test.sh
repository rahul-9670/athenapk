#!/bin/bash
#SBATCH --job-name=flagtest
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/prod_flagship_test/%x_%j.out
set -o pipefail
#
# FLAGSHIP PRODUCTION TEST (2026-07-30) -- one bounded slot, NOT self-chaining.
#
# Purpose: exercise the full production physics stack with the four flagship-audit items
# switched ON, before committing the real campaign:
#   item 1  diffusion/cap_diag=true     -- eta-cap activation accounting
#   item 2  hydro/mag_diag=true         -- J^2, current helicity, ME_tor/ME_pol, eta-weighted diss
#   item 5  diffusion/dust_coupling=true + <dust> evolve=true, on the pgen dust-IC fix
#   (item 4, the Omega.B polarity axis, is an ensemble-design change and is not exercised here)
#
# DIVERGENCE CONTROL: GLM/Dedner (fluid=glmmhd, no divergence_control key). CT is NOT used --
# measured cost of dropping CT is 1.6% on flux retention / 2.0% on mu_core (runs/inc7_gate).
#
# Deliberately NOT self-chaining: this is a diagnostic slot. It writes a restart, so the real
# campaign can either continue from it or start fresh.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WDIR=/beegfs/u/bbg6470/athenapk/runs/prod_flagship_test
WRAP=$WDIR/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

cd $WDIR
echo "=== flagship production test, job $SLURM_JOB_ID, $(date) ==="
echo "binary:"; md5sum $BIN

# Resume if a restart exists (so a re-submit after TIMEOUT continues), else fresh.
LATEST=$(ls -t $WDIR/parthenon.out2.*.rhdf 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then RA="-r $LATEST"; echo "RESUMING FROM $LATEST"
else RA="-i fhc_flagship.in"; echo "FRESH START FROM t=0"; fi

echo "=== start $(date) job $SLURM_JOB_ID ===" >> $WDIR/run.log
( while true; do echo "$(date +%s) $(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits | tr "\n" " ")" >> $WDIR/gpumem.log; sleep 60; done ) &
GPULOG_PID=$!
# Same CLI overrides as the mg_prod_tab production campaign, so throughput is comparable.
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN $RA -t 03:30:00 \
  parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 parthenon/output2/dn=250 \
  >> $WDIR/run.log 2>&1
RC=$?; kill $GPULOG_PID 2>/dev/null; echo "RUN_EXIT $RC $(date)"
tail -5 $WDIR/run.log
