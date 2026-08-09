#!/bin/bash
#SBATCH --job-name=wp8split
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=00:40:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp8_dissplit/%x_%j.out
set -o pipefail
#
# WP-08 -- MEASURE THE dissO/dissA DENSITY SPLIT ON THE LADDER. Stage 1 of 2: nj4 only.
#
# THE STATE OF WP-08. It is the last work package still OPEN. mag-Jsq does not converge on the
# njeans ladder (-52.3 %, -39.3 % per rung) and mag-dissO is not even monotone (-76.5 %, +40.5 %),
# while mag-ME converges to 0.2 %. mag_diag.hpp's own analysis says why: these "volume integrals"
# are effectively POINT SAMPLES -- 90 % of int|J|^2 dV comes from a volume fraction of 7.6e-7 to
# 4.7e-6 -- taken in exactly the region whose resolution changes between rungs. Its prescribed
# remedy is the DENSITY SPLIT: integrate over regions defined by physics, not by the grid, so the
# core and envelope budgets converge (or fail) independently and visibly.
#
# Two splits have now been TRIED and both FAILED, which is why this one is worth the GPU time:
#   * density split applied to Jsq (2026-08-06): the low-density bin carries 97-99 % of the
#     integral at f_eff ~ 1e-7, i.e. the pathology under a new name (-51.6 %/-39.9 % vs the
#     global -52.3 %/-39.3 %). Its conclusion: Jsq needs a CURRENT-SHEET indicator, not density.
#   * current-sheet split applied to Jsq (2026-08-08, measured today): no bin converges;
#     smooth@0.1 is not monotone (+139.6 %, -28.8 %) and smooth@0.3 only reproduces the global
#     because its sheet bin is EMPTY (0 cells at nj8 and nj16).
# Both of those split **Jsq**, computed offline from the snapshots. NEITHER split has ever been
# applied to **dissO/dissA**, which is the quantity mag_diag.hpp actually prescribed the density
# split for -- because eta is not in any snapshot, so it cannot be done offline (see below).
#
# WHY THIS NEEDS A RUN AT ALL, AND WHY IT IS NOT A LADDER RE-RUN.
# The split columns ALREADY EXIST in the code (mag-dissO-hi/-lo, mag-dissA-hi/-lo, mag-Vhi,
# mag-dissOV/AV, hydro.cpp:2069+), gated on hydro/mag_diag_rho_split > 0. The ladder ran
# hydro/mag_diag=true and diffusion/cap_diag=true but NEVER set mag_diag_rho_split, so its .hst
# files carry only the global cols 27/28. Another "implemented and gated, but zero decks opted
# in" case.
# It cannot be recovered offline: `nonideal_eta` is registered {Cell, Derived, OneCopy} with NO
# Restart flag (hydro.cpp:1741), so eta is in neither the .phdf nor the .rhdf. Recomputing it in
# Python is blocked too -- the deck sets diffusion/dust_coupling=true, so eta depends on the
# EVOLVED grain state (f_dg, a_c), and no dust field appears in any output block. x_e alone
# (prim_scalar_4) is not enough.
# But `Derived` is exactly what makes this cheap: eta is RECOMPUTED FROM STATE at init, so simply
# RESTARTING regenerates it. No re-run needed -- ~14 GPU-h for all three rungs against ~180 to
# re-run the ladder.
#
# WHY nj4 FIRST (staged). nj4's restart out2.00001 sits ON the matched epoch: t=1.093652 vs the
# matched-epoch snapshot's t=1.093652, |dt| = 1.07e-07, rho_max = 9.797e-13 g/cm3 (+0.009 dex
# from the 1e-12 target). So it needs ~2 cycles, not ~205 (nj8) or ~113 (nj16). It therefore
# delivers a real data point AND verifies the premise the other two rungs depend on:
#   mag_diag_rho_split is a GetOrAddReal (hydro.cpp:2008), so the value 0.0 is written into the
#   restart file's stored deck. Whether a CLI override beats the stored deck on a RESTART is
#   exactly the mechanism that killed D1 job 2461421 in 9 s (existence guards are restart-
#   hostile: GetOrAdd defaults get baked into the restart). Verify it on the 0.8-GPU-h leg, not
#   on the 13-GPU-h pair.
#
# BINARY = the ladder's own, md5 5ebddce0, confirmed against "binary:" in conv_nj4_2431681.out.
# NOT today's default-flip binary: N13/N14 repointed the opacity and EOS table defaults, which
# would silently change the physics relative to the three rungs this must be compared against.
#
# WRITES TO A FRESH DIRECTORY. The ladder run dirs are NOT touched -- their newest restarts are
# the resume points and their out1 series is the published WP-08 data. cwd is $R/nj4, the restart
# is read by absolute path, and output1/output2 are pushed past the run length so this writes no
# 12 GB snapshot for a 2-cycle diagnostic.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

R=/beegfs/u/bbg6470/athenapk/runs/wp8_dissplit
L=/beegfs/u/bbg6470/athenapk/runs/convergence_ladder
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_PRESERVED_5ebddce0
RST=$L/nj4/parthenon.out2.00001.rhdf          # t=1.093652 cyc=250 -- ON the matched epoch
NJ=4
# rho_crit = 1e-13 g/cm3 in code units (rho0 = 5.467e-19), the same split
# wp8_split_convergence.py uses for the Jsq density split, so the two are directly comparable.
SPLIT=182915.67587342235

D=$R/nj4; mkdir -p $D; cd $D
cp -n $L/nj4/wrap_mod.sh $D/ 2>/dev/null || cp -n /beegfs/u/bbg6470/athenapk/runs/ensemble/design01/point000/wrap_mod.sh $D/
WRAP=$D/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

echo "=== WP-08 dissO/dissA split, leg nj$NJ, job $SLURM_JOB_ID $(date) host=$(hostname) ==="
echo "binary md5 $(md5sum $BIN | cut -c1-8)   (ladder binary = 5ebddce0)"
echo "restart    $RST"
echo "rho_split  $SPLIT code = 1e-13 g/cm3"
echo "--- allocated GPUs: SLURM_JOB_GPUS=${SLURM_JOB_GPUS:-<unset>} CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"
nvidia-smi --query-gpu=index,memory.used,memory.total --format=csv,noheader | sed 's/^/    pre-run /'

# nlim is a TOTAL cycle count, and the restart is at cycle 250 -> 256 runs 6 cycles. Six, not two,
# so the .hst has enough rows to see whether the split columns are steady rather than a single
# value that could be anything. output0/dt=1e-9 is far below the local dt (5.2e-5 at cycle 250)
# and so forces one history row PER CYCLE.
stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN -r $RST -t 00:30:00 \
  refinement/njeans=$NJ parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
  diffusion/cap_diag=true hydro/mag_diag=true \
  hydro/mag_diag_rho_split=$SPLIT \
  parthenon/time/nlim=256 \
  parthenon/output0/dt=1.0e-9 \
  parthenon/output1/dn=1000000 parthenon/output2/dn=1000000 \
  > $D/run.log 2>&1
echo "RUN_EXIT $? $(date)"

# ---- THE GATE: did the CLI override actually win over the restart's stored deck? ----
echo; echo "=== GATE 1: does the code report the split as ON? ==="
grep -a "mag_diag density split" $D/run.log | sed 's/^/    /' || echo "    (no line -- SPLIT DID NOT ENGAGE)"
echo "=== GATE 2: are the split columns in the .hst header? ==="
head -3 $D/parthenon.out0.hst 2>/dev/null | tr ' ' '\n' | grep -n "mag-" | sed 's/^/    /'
echo "=== GATE 3: are they non-zero, and do hi+lo reconstruct the global? ==="
/beegfs/u/bbg6470/venvs/analysis_env/bin/python $R/check_split.py $D/parthenon.out0.hst
echo "cycles_run=$(grep -ac '^cycle=' $D/run.log)"
echo WP8_SPLIT_DONE
