#!/bin/bash
# Launch every ensemble member of a design (one self-chaining GPU job per point).
# Usage: ./launch_ensemble.sh <design_dir> [STOP_CGS]   (default STOP_CGS=1e-13 = first-core epoch)
# Throttle-friendly: submits all points; SLURM queues them. prod_v9 untouched (separate dirs+binary).
set -euo pipefail
DDIR="${1:?usage: launch_ensemble.sh <design_dir> [STOP_CGS]}"
STOP_CGS="${2:-1.0e-13}"
HERE=/beegfs/u/bbg6470/athenapk/runs/ensemble
for pdir in "$DDIR"/point*/; do
  [ -f "$pdir/fhc_ens.in" ] || continue
  jid=$(sbatch --job-name="ens_$(basename "$pdir")" \
        --export=ALL,PDIR="$(readlink -f "$pdir")",STOP_CGS="$STOP_CGS" \
        "$HERE/submit_point.sh" | awk '{print $NF}')
  echo "  $(basename "$pdir") -> job $jid"
done
echo "launched $(ls -d "$DDIR"/point*/ 2>/dev/null | wc -l) ensemble members (STOP at rho=$STOP_CGS cgs)"
