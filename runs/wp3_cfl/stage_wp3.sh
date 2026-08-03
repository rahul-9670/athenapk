#!/bin/bash
# WP-3 — temporal (CFL) convergence of the PRODUCTION configuration.
#
# WHAT WAS MISSING. `parthenon/time/cfl = 0.3` has never been tested against a smaller value on
# the production deck. The one prior attempt used a cheap CPU deck that was out of equilibrium and
# reported KE differing by a FACTOR OF 46 between cfl = 0.3 and 0.0375 -- a number that says
# nothing about production and was never transferable.
#
# DESIGN. Three legs at 256^3 on the new production physics (binary f181c0a1, diode fluid BCs,
# relative gravity residual), differing ONLY in `parthenon/time/cfl`:
#
#     cfl = 0.3    -- production; this leg ALREADY EXISTS as root_ladder/r256_sw, so it is not
#                     re-run. Reusing it is what makes this study cheap.
#     cfl = 0.15   -- half
#     cfl = 0.075  -- quarter
#
# Only `parthenon/time/cfl` is varied. `diffusion/cfl` and `radiation/cfl` are deliberately left
# at their production values (0.3, 0.4): those govern the explicit non-ideal and M1 substeps, and
# changing them too would confound "is the time integration converged" with "are the operator-split
# substep counts converged", which are separate questions. Note that under RKL2 the diffusion
# substep structure still follows dt through `rkl2_max_dt_ratio`, so halving the hydro CFL does
# change the STS pattern -- that is a property of the production configuration and is in scope.
#
# JUDGE AT t = 0.90, NOT AT THE t = 1.0 ENDPOINT. WP-7 established that the last 0.01 t0 of these
# uniform-grid runs sits on the collapse singularity (dt falls 1-2 decades) and that comparisons
# taken there measure where each leg happened to stall. See docs/validation/WP07_root_grid_ladder.md.
#
# COST. r256_sw runs ~40 min at cfl = 0.3 on 6 GPUs, so budget ~2x and ~4x that for the two legs.
# Both use 3 GPUs, not 6 -- WP-12 makes rank count a free parameter, and `MaxNodes=1` means a job
# cannot span nodes, so a 6-GPU request sits unschedulable whenever no single node has 6 free.
# Memory is not the binding constraint: r256_gfix measured 13.5 GiB/GPU on 6 ranks, i.e. ~27
# GiB/GPU on 3, against an 80 GiB card.
#
# Usage:  ./stage_wp3.sh        -> print
#         ./stage_wp3.sh go     -> submit
set -o pipefail
LAD=/beegfs/u/bbg6470/athenapk/runs/root_ladder
HERE=/beegfs/u/bbg6470/athenapk/runs/wp3_cfl
WRAPSRC=/beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh
BIN=/beegfs/u/bbg6470/athenapk/build_gpu_v2/bin/athenaPK
[ -r "$WRAPSRC" ] || { echo "FATAL: no GPU-pinning wrapper at $WRAPSRC"; exit 1; }
[ -x "$BIN" ]     || { echo "FATAL: no binary at $BIN"; exit 1; }
echo "binary: $(md5sum $BIN)"
echo "deck:   $(md5sum $LAD/fhc_rootladder.in)"
echo "reference leg (cfl=0.3) = $LAD/r256_sw -- not re-run"

for CFL in 0.15 0.075; do
  DIR=$HERE/cfl${CFL}
  mkdir -p $DIR
  install -m 755 $WRAPSRC $DIR/wrap_mod.sh || { echo "FATAL: cannot install wrapper"; exit 1; }
  CMD="sbatch --job-name=wp3c${CFL} --nodes=1 --ntasks=3 --gres=gpu:h100:3 \
--export=ALL,NX=256,RUNDIR=$DIR,NRANK=3,BIN=$BIN,OV=parthenon/time/cfl=$CFL \
$LAD/submit_root.sh"
  if [ "${1:-}" = "go" ]; then echo "+ $CMD"; eval $CMD; else echo "would run: $CMD"; fi
done
