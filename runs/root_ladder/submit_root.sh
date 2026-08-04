#!/bin/bash
#SBATCH --job-name=rootlad
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --output=%x_%j.out
# MUST be set explicitly: the gpu partition's DefaultTime is 01:00:00, while line ~70 passes
# `-t 05:30:00` to athenaPK. Without this the code plans for 5.5 h and SLURM kills at 1 h -- a
# 5.5x mismatch that stayed invisible while every leg happened to finish inside the hour.
# It surfaced when WP-2's cr300 leg (chat=999 km/s, ~3.3x the radiation substeps) was killed at
# exactly 01:00:03 with 80 of ~181 cycles done, and would have killed the 512^3 leg too (the
# original r512 took 3.06 h). Keep this >= the `-t` value below.
#SBATCH --time=06:00:00
set -o pipefail
#
# ROOT-GRID CONVERGENCE LADDER (2026-07-31, option "a").
#
# WHY THIS EXISTS: runs/convergence_ladder/ varies refinement/njeans (4/8/16) at a FIXED 256^3
# root grid. Measured on those runs: all three legs sit on the identical uniform 512-block root
# grid until t=0.91 (nj16) / 0.999 (nj8) / 1.049 (nj4) -- i.e. 83-96% of the evolved time, and
# the phase in which essentially all of the magnetic braking happens, is run at ONE resolution
# and is therefore untested. That phase costs only 5.8-9.9% of the wall time, so testing it is
# cheap. This ladder varies the BASE resolution (128/256/512) on a UNIFORM grid to t=1.0.
#
# The two knobs are close to orthogonal: refinement here is Jeans-driven, so the cell size
# ACHIEVED at a given density is set by njeans and is nearly independent of the root grid --
# njeans isolates the core, the root grid isolates the envelope.
#
# Parameterized by env:
#   NX     = root cells per dimension (128 / 256 / 512)
#   RUNDIR = this leg's directory
#   NRANK  = MPI ranks = GPUs (see the memory note below)
#   MBLK   = OPTIONAL meshblock edge (default: whatever the deck says, 32). Needed for the
#            512^3 rung: at 32^3 blocks it projects to 82.8 GiB/GPU against an 80 GiB card and
#            does NOT fit; 64^3 blocks carry proportionally fewer ghost cells and bring it to
#            69.8 GiB. Legitimate ONLY because this study is a UNIFORM grid with no AMR, so
#            meshblock size is pure domain decomposition, not discretization -- and WP-12
#            measured decomposition to have no effect (8 vs 4 ranks agreed to all 7 printed
#            figures on every history quantity across 179 rows).
#
# MEMORY SIZING (anchored on a measurement, not a guess): prod_flagship_test job 2431680 peaked
# at 47.0 GiB/GPU with 1408 blocks of 32^3 on 5 GPUs = 9.23e6 cells/GPU => 5.09 KiB/cell. That
# anchor is CONSERVATIVE for this study because it includes the deep-AMR multigrid hierarchy and
# AMR buffers, which a uniform grid does not carry. Projected peak per GPU:
#   128^3 =   64 blocks, 8 ranks  ->  8 blk/rank -> ~1.3 GiB
#   256^3 =  512 blocks, 8 ranks  -> 64 blk/rank -> ~10.7 GiB
#   512^3 = 4096 blocks, 16 ranks -> 256 blk/rank -> ~42.7 GiB   (on 8 ranks it would be
#                                    ~85 GiB > the 80 GiB H100 -- hence 2 nodes for this leg.)
#
# Single slot, NOT self-chaining (each leg is short). Resumes from the newest restart if one
# exists, else fresh from t=0.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

: "${NX:?set NX}"; : "${RUNDIR:?set RUNDIR}"; : "${NRANK:?set NRANK}"
# BIN is overridable so the ladder can be re-run on the WP-20 candidate. The deck now sets
# turb_ksample=k2, which binary 5ebddce0 does NOT know and silently ignores -- running the
# old binary against the new deck yields the OLD IC with no warning. The banner check below
# is the guard: it aborts rather than quietly producing a superseded initial condition.
# 2026-08-04: default swapped build_gpu_wp20 -> build_gpu_v4 (869c1d34 = v3 + B1 + B6).
# OFF-state gate-verified byte-identical to v3 on all 11 STATE columns; see submit_deep.sh for the
# full reasoning and the last-printed-digit diagnostic caveat. B1's clamp is DEFAULT OFF.
# NOTE this is only the FALLBACK default -- stage_switches.sh exports BIN explicitly and is
# deliberately still pinned to build_gpu_v2 (f181c0a1), because r128_sw and r256_sw already ran on
# that binary and the WP-7 root-grid ladder must stay a matched comparison. Do not "modernise" it.
BIN="${BIN:-/beegfs/u/bbg6470/athenapk/build_gpu_v4/bin/athenaPK}"
DECK=/beegfs/u/bbg6470/athenapk/runs/root_ladder/fhc_rootladder.in
WRAP=$RUNDIR/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
cd $RUNDIR

echo "root ladder NX=$NX NRANK=$NRANK (job $SLURM_JOB_ID) $(date)"
echo "nodes: $SLURM_JOB_NODELIST"
echo "binary:"; md5sum $BIN
echo "deck:"; md5sum $DECK
grep -q "Driver completed" run.log 2>/dev/null && { echo "already completed -> exit"; exit 0; }
LATEST=$(ls -t $RUNDIR/parthenon.out2.*.rhdf 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then RA="-r $LATEST"; echo "resuming from $LATEST"; else RA="-i $DECK"; fi

echo "=== NX=$NX start $(date) ===" >> $RUNDIR/run.log
stdbuf -oL -eL mpirun -n $NRANK $MCA $WRAP $BIN $RA -t 05:30:00 \
  parthenon/mesh/nx1=$NX parthenon/mesh/nx2=$NX parthenon/mesh/nx3=$NX \
  ${MBLK:+parthenon/meshblock/nx1=$MBLK parthenon/meshblock/nx2=$MBLK parthenon/meshblock/nx3=$MBLK} \
  parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 parthenon/output2/dn=250 \
  diffusion/cap_diag=true hydro/mag_diag=true \
  diffusion/hall_ohmic_floor_ratio=0.2 \
  hydro/boundary_flux_clamp=${CLAMP:-true} \
  ${OV:-} \
  >> $RUNDIR/run.log 2>&1
# OV = optional extra CLI overrides, appended LAST so they win. Empty when unset, so every
# existing caller (the root ladder, stage_gfix.sh) is unaffected. Used by WP-2's creduc sweep
# (runs/wp2_creduc/stage_wp2.sh) to reuse this exact configuration with one key varied.
#
# CLAMP (2026-08-04): B1's boundary-flux clamp is now ON BY DEFAULT for new production runs.
# Rationale -- production already uses DIODE boundaries, which block inflow VELOCITY, yet mass
# still climbs **+0.0975 %** over a full pre-collapse run (deep_amr, t = 0 -> 1.015, 204 rows)
# because the Riemann solver still yields an inward FLUX from the pressure/density gradient.
# Diode does not cover that; the clamp does, driving the drift to exactly 0.000000 % (GPU gate
# job 2452598). It is gated OFF-state byte-identical and falsified against a uniform wind, which
# confirmed OUTWARD flux passes untouched at exactly rho*v*A.
#
# NOT applied retroactively. Completed work (WP-7 ladder, WP-18 ensemble, WP-8) is internally
# consistent with the clamp OFF, and +0.0975 % mass is ~150x below WP-18's sigma = 16 % on the
# flux-retention observable, so re-running would cost thousands of GPU-hours to move nothing.
# `stage_switches.sh` therefore exports CLAMP=false to keep the WP-7 root-grid ladder matched --
# r128_sw and r256_sw already ran without it. Job 2449283 (root512sw) is unaffected either way:
# sbatch copies the batch script at submission, so it runs the pre-edit version.
echo "RUN_EXIT $? $(date)"

# GUARD: the corrected turbulence IC must actually have been used. If the binary predates
# turb_ksample it ignores the key silently and the run is on the SUPERSEDED IC.
if grep -q "k sampling *: k2" $RUNDIR/run.log; then
  echo "IC_CHECK ok: k2 sampler active ($(grep -m1 'E(k) ~' $RUNDIR/run.log | tr -s ' '))"
else
  echo "IC_CHECK **FAILED**: run.log does not report the k2 sampler."
  echo "  The binary likely predates turb_ksample and silently used the OLD IC. Results void."
  echo "ic-check-failed $(date)" > $RUNDIR/STOP_CHAIN
fi
