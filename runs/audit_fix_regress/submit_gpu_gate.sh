#!/bin/bash
#SBATCH --job-name=auditgpu
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=%x_%j.out
#
# GPU EVIDENCE for the parts of the 2026-08-05 audit batch that have none.
#
# WHY THIS EXISTS. docs/validation/gate_v5_audit_batch.md establishes that the batch is INERT on
# the production path (v5 vs v4: every reproducible history column bit-identical). But it also
# records, honestly, that two of the batch's SUBSTANTIVE fixes were never actually exercised on
# GPU, because the legs chosen for them could not express the fix:
#
#   * A1 (EOS-consistent refinement) was run at numlevel=3, njeans=8 -- and nbtotal sat at a
#     constant 64 in BOTH legs, i.e. THE MESH NEVER REFINED. "Identical" there is a null result.
#     That is the same trap already on record for A_multipole.in / B_swindle.in.
#   * N2/N3 (dust sublimation temperature from the EOS table; dust unit consolidation) was run in
#     a configuration with no dust column in the history file and gas at T ~ 10 K, while N2 only
#     bites at the ~1500 K sublimation threshold.
#
# This job runs the decks that DO discriminate -- the ones that already exist and were designed
# for exactly this -- so the two claims stop resting on CPU evidence alone.
#
#   F_a1discrim.in   njeans = 4.7 sits BETWEEN the ideal criterion (min n_J 4.537) and the
#                    tabulated one (4.925), so the two formulas provably disagree. The block count
#                    that separates them is the INITIAL one (cycle 0), not the final one:
#                        pre-fix  nbtotal = 120 120 120 120 120
#                        post-fix nbtotal =  64  64  64 120 120
#                    i.e. the fix stops the initial over-refinement, and by cycle 3 the collapse
#                    has legitimately driven BOTH to 120. v5 HAS the fix => expect 64 at cycle 0.
#                    (2026-08-05 CORRECTION: this check previously read the LAST block count out of
#                    the run log and reported "A1 FIX MISSING ON GPU" on job 2461424. That was a
#                    false alarm -- 120 is the correct final count in both legs. Verified against
#                    results/a1discrim_a1old vs a1discrim_a1new; the GPU history matched the
#                    post-fix reference exactly on the block sequence and to 1.2e-06 on the
#                    momentum sums. Read column 4 of the .hst, row 0.)
#   H_dust_ideal.in  dust ON with eos = adiabatic  -- exercises the dust package's own path.
#   I_dust_hydrogen  dust ON with eos = hydrogen   -- exercises the N2 table-temperature branch.
#   E_sinkmhd.in     sinks under MHD. This is A4's path: the pre-fix code scaled total energy by
#                    (1-phi) while leaving B alone, charging accreted gas for magnetic energy
#                    still in the cell, and ABORTED ON CYCLE 0 with "Got negative pressure" when
#                    ME > e_th. Nobody had ever run sinks under MHD, which is why it survived.
#
# WHAT THIS GATE IS NOT. It is not a bit-identity gate against the CPU references. GPU and CPU
# differ in reduction and atomic ordering -- gate_v5_audit_batch.md measured v4 failing to
# reproduce ITSELF on the cons-Pout* surface integrals -- so byte equality is the wrong test and
# would return a meaningless FAIL. These are EXECUTION + STRUCTURE gates: the run completes, no
# NaN, and the one structural number that distinguishes fixed from unfixed comes out right.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

H=/beegfs/u/bbg6470/athenapk/runs/audit_fix_regress
BIN=/beegfs/u/bbg6470/athenapk/build_gpu_v5/bin/athenaPK
WRAP=/beegfs/u/bbg6470/athenapk/runs/deep_amr/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

[ -x "$BIN" ] || { echo "MISSING $BIN -- aborting"; exit 1; }
echo "job $SLURM_JOB_ID $(date)"; md5sum $BIN
cd /beegfs/u/bbg6470/athenapk && git rev-parse HEAD

run_deck () {   # deck
  local D=$1
  local G=$H/gpu_$D; rm -rf $G; mkdir -p $G
  echo "=== $D ==="
  ( cd $G && stdbuf -oL -eL mpirun -n 1 $MCA $WRAP $BIN -i $H/$D.in > run.log 2>&1 )
  local RC=$?
  local NCYC=$(grep -c '^cycle=' $G/run.log)
  local NAN=$(grep -ci "nan\|inf" $G/run.log)
  echo "  exit=$RC cycles=$NCYC nan_or_inf_mentions=$NAN"
  # POSITIVE CHECK: exit=0 is not evidence on this cluster -- runs have exited 0 after dying in
  # seconds. A leg counts only if it actually cycled.
  if [ "$NCYC" -eq 0 ]; then
    echo "  **VOID: no cycles. First error --**"
    grep -A4 -iE "PARTHENON ERROR|Kokkos ERROR|what\(\):|negative pressure" $G/run.log |
      head -8 | sed 's/^/     /'
  fi
  grep -iE "Kokkos ERROR|failed to allocate|negative pressure" $G/run.log | head -3 |
    sed 's/^/  ERR: /'
}

for d in F_a1discrim H_dust_ideal I_dust_hydrogen E_sinkmhd; do run_deck $d; done

echo
echo "=== A1 STRUCTURAL CHECK (the discriminating number) ==="
# The fix changes the refinement criterion's sound speed from the ideal gamma to the EOS table.
# At njeans = 4.7 the two disagree at t=0, so the INITIAL block count separates them: 120 pre-fix,
# 64 post-fix. Both legs reach 120 by cycle 3 through genuine collapse, so the FINAL count is not a
# discriminator -- reading it is what produced the false "FIX MISSING" on job 2461424.
# .hst column 4 = nbtotal; first non-comment row = cycle 0.
NB=$(grep -v '^#' $H/gpu_F_a1discrim/parthenon.out0.hst | head -1 | awk '{printf "%d", $4}')
echo "  initial nbtotal = ${NB:-UNKNOWN}   (expect 64 = fix present; 120 = fix absent on GPU)"
if   [ "$NB" = "64"  ]; then echo "  => A1 CONFIRMED ON GPU"
elif [ "$NB" = "120" ]; then echo "  => A1 FIX MISSING ON GPU -- investigate the compile"
else echo "  => INCONCLUSIVE: could not parse the block count; read the log by hand"; fi

echo
echo "=== E_sinkmhd A4 CHECK ==="
if grep -qi "negative pressure" $H/gpu_E_sinkmhd/run.log; then
  echo "  => A4 REGRESSION: the MHD sink path still produces negative pressure"
else
  echo "  => no negative pressure under MHD sinks (A4 fix holds on GPU)"
fi
echo "done $(date)"
