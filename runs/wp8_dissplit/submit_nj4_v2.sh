#!/bin/bash
#SBATCH --job-name=wp8nj4v2
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=00:50:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp8_dissplit/%x_%j.out
set -o pipefail
#
# WP-8 nj4, ATTEMPT 2 -- on the BACKPORTED ladder binary.
#
# Attempt 1 (job 2495290) failed for a reason that had nothing to do with the restart: I ran
# athenaPK_PRESERVED_5ebddce0, the ladder's own binary, which was built 2026-07-30 while the split
# columns were written 2026-07-31+. `strings` confirms the key is simply not in it. ParameterInput
# accepts any CLI key and merely WARNS that it is adding it, so an unknown key is silently inert
# and looks exactly like one that was honoured. See submit_dissplit.sh for the full post-mortem.
#
# The fix was to backport the ~74 non-comment lines of read-only diagnostic onto the WP-0
# reproduction tree (commit 29a7174 + archived patch), which WP-0 had already proven rebuilds to a
# BYTE-IDENTICAL parthenon.out0.hst against 5ebddce0. That gives ladder physics AND the
# diagnostic, avoiding the audit-fix chain -- specifically A1, which changes the refinement
# criterion, i.e. exactly the variable a resolution-convergence study measures.
#
# TWO LEGS, because a backport claim needs a falsifier as well as a result.
#
#   off : backported binary, NO split keys. Its .hst must be BYTE-IDENTICAL to the pre-backport
#         binary's (athenaPK_prewp8, saved by build_wp8.sh before the rebuild). This is the whole
#         basis for restarting the ladder with a modified binary: the claim "read-only reductions,
#         bit-identical when off" is ASSERTED in mag_diag.hpp and has never been tested for this
#         backport. If it fails, the measurement is contaminated and must not be used -- and
#         nj8/nj16 must not be launched.
#   on  : backported binary + both split thresholds. The actual measurement.
#
# The off leg runs the ORIGINAL pre-backport binary too, so the comparison is a real A/B on this
# exact deck and restart rather than a comparison against a number from another day.
#
# rho_split = 1.83e5 code = 1e-13 g/cm3 = rho_crit, the same threshold
# docs/validation/scripts/wp8_split_convergence.py uses for the Jsq density split, so the in-code
# and offline splits are directly comparable.
# sheet_thresh = 0.1: MEASURED as the only value that discriminates at production resolution
# (0.3 selects 0.49 % at nj4 and exactly 0 % at nj8/nj16 -- a null diagnostic that looks like a
# working one). See WP08_dissipation_nonconvergence.md.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

R=/beegfs/u/bbg6470/athenapk/runs/wp8_dissplit
L=/beegfs/u/bbg6470/athenapk/runs/convergence_ladder
BLD=/beegfs/u/bbg6470/wp0_repro/athenapk/build_gpu_repro/bin
NEW=$BLD/athenaPK            # backported
OLD=$BLD/athenaPK_prewp8     # pre-backport (= WP-0's reproduction of the ladder binary)
RST=$L/nj4/parthenon.out2.00001.rhdf   # t=1.093652 cyc=250 -- ON the matched epoch
SPLIT=182915.67587342235
SHEET=0.1

echo "=== WP-8 nj4 attempt 2, job $SLURM_JOB_ID $(date) host=$(hostname) ==="
echo "NEW $(md5sum $NEW | cut -c1-12)"
echo "OLD $(md5sum $OLD | cut -c1-12)"
echo "--- allocated GPUs: SLURM_JOB_GPUS=${SLURM_JOB_GPUS:-<unset>}"
nvidia-smi --query-gpu=index,memory.used --format=csv,noheader | sed 's/^/    pre-run /'

MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
COMMON="refinement/njeans=4 parthenon/mesh/do_coalesced_comms=true
        diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2
        diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true
        diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16
        diffusion/cap_diag=true hydro/mag_diag=true
        parthenon/time/nlim=256 parthenon/output0/dt=1.0e-9
        parthenon/output1/dn=1000000 parthenon/output2/dn=1000000"

run() { # run <dir> <binary> [extra args...]
  local d=$R/$1 b=$2; shift 2
  rm -rf $d; mkdir -p $d; cd $d
  cp -n $L/nj4/wrap_mod.sh $d/ 2>/dev/null
  stdbuf -oL -eL mpirun -n 5 $MCA $d/wrap_mod.sh $b -r $RST -t 00:20:00 \
    $COMMON "$@" > $d/run.log 2>&1
  echo "  $1 exit=$? cycles=$(grep -ac '^cycle=' $d/run.log)"
}

run v2_off_new $NEW
run v2_off_old $OLD
run v2_on      $NEW hydro/mag_diag_rho_split=$SPLIT hydro/mag_diag_sheet_thresh=$SHEET

cd $R
echo
echo "=== GATE A (falsifier): OFF state must be BYTE-IDENTICAL across the backport ==="
if cmp -s $R/v2_off_new/parthenon.out0.hst $R/v2_off_old/parthenon.out0.hst; then
  echo "    PASS -- backported binary reproduces the pre-backport .hst byte for byte"
else
  echo "    ** FAIL -- the backport CHANGED the off-state result. Do NOT use this measurement,"
  echo "    ** and do NOT launch nj8/nj16. First differing bytes:"
  cmp $R/v2_off_new/parthenon.out0.hst $R/v2_off_old/parthenon.out0.hst | sed 's/^/       /'
fi
echo "=== GATE B: did the split engage? ==="
grep -a "mag_diag .*split ON" $R/v2_on/run.log | sed 's/^/    /' || echo "    (no line -- SPLIT DID NOT ENGAGE)"
echo "=== GATE C: columns present, both bins populated, hi+lo == global ==="
/beegfs/u/bbg6470/venvs/analysis_env/bin/python $R/check_split.py $R/v2_on/parthenon.out0.hst
echo "=== bonus: in-code sheet split vs the offline Python (83.22 % grid-scale share at nj4) ==="
/beegfs/u/bbg6470/venvs/analysis_env/bin/python - $R/v2_on/parthenon.out0.hst <<'PY'
import sys, re, numpy as np
p=sys.argv[1]; cols=[]
for line in open(p):
    if line.startswith("#") and "[" in line: cols=re.findall(r"\[\d+\]=(\S+)",line); break
d=np.loadtxt(p,comments="#")
d=d[None,:] if d.ndim==1 else d
try:
    g=lambda n: d[-1, cols.index(n)]
    sh, sm = g("mag-Jsq-sheet"), g("mag-Jsq-smooth")
    tot=sh+sm
    print(f"    sheet={sh:.4e}  smooth={sm:.4e}  grid-scale share={sh/tot*100:.2f} %"
          f"   (offline Python got 83.22 % at nj4)")
except ValueError:
    print("    sheet columns absent")
PY
echo WP8_NJ4V2_DONE
