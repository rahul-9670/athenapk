#!/bin/bash
# WP-2 — reduced-speed-of-light insensitivity.  Staged 2026-08-02; launch with ./stage_wp2.sh go
#
# Semantics settled first, in docs/validation/WP02_creduc_semantics.md:
#   creduc is a DIVISOR.  chat = c_code/creduc, c_code = c_light/v_unit = 1.5779e6 code v.
#     creduc=300  -> chat =  999 km/s   (MORE accurate, MORE expensive)
#     creduc=1000 -> chat =  300 km/s   (production)
#     creduc=3000 -> chat =  100 km/s   (MORE approximate, cheapest)
#   Convergence is toward SMALLER creduc.  Acceptance = |300 leg - 1000 leg| < WP-18 sigma (16%).
#
# Everything else is identical to the WP-7 root-ladder configuration at 256^3 uniform, and all
# legs share one seed, so this is a PAIRED comparison.
#
# COST IS NOT EQUAL ACROSS LEGS: chat sets the radiation CFL, so the 300 leg takes ~3.3x more
# radiation substeps per hydro step than production.  Launch 3000 and 1000 first; 300 last.
set -euo pipefail
HERE=/beegfs/u/bbg6470/athenapk/runs/wp2_creduc
BIN=${BIN:-/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK}

# `grep -qa` on the file, not `strings | grep -q`: under pipefail grep -q's early exit gives
# strings a SIGPIPE and the check would fail on a GOOD binary.
grep -qa grav_rho_floor $BIN \
  || { echo "ABORT: $BIN predates the WP-13 gravity fix -- not comparable with the re-baseline"; exit 1; }
echo "binary: $(md5sum $BIN)"

for CR in 3000 1000 300; do
  D=$HERE/cr$CR
  mkdir -p $D
  cp -n $HERE/fhc.in $D/ 2>/dev/null || true
  # Wrapper source is runs/wrap_mod.sh (mode 0644 there -> must be installed 755). The old path
  # root_ladder/wrap_mod.sh does not exist and the `|| true` hid that: job 2446543 died in 1 s
  # with "prterun ... lacked permissions to execute". Hard-fail instead.
  install -m 755 /beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh $D/wrap_mod.sh
  printf 'cr%-6s creduc=%-6s chat=%.0f km/s\n' "$CR" "$CR" \
    "$(/beegfs/u/bbg6470/venvs/analysis_env/bin/python -c "print(2.99792458e10/1.9e4/$CR*1.9e4/1e5)")"
  if [ "${1:-}" = "go" ]; then
    # 4 ranks, not 8. Decomposition does not affect results (WP-12: 8 vs 4 agreed to all 7
    # printed figures on 179 rows), 256^3 at 4 ranks is 20.4 GiB/GPU against an 80 GiB card, and
    # 4 ranks uses 0.74x the GPU-hours of 8 (measured: 15.3 s/cyc x 4 vs 10.3 s/cyc x 8).
    # Smaller jobs also pack into partial-node gaps instead of waiting for a whole free node.
    # --ntasks/--gres must be on the CLI: submit_root.sh's header sets neither.
    sbatch --job-name=wp2cr$CR --nodes=1 --ntasks=4 --gres=gpu:h100:4 \
      --export=ALL,NX=256,RUNDIR=$D,NRANK=4,BIN=$BIN,OV="radiation/creduc=$CR" \
      /beegfs/u/bbg6470/athenapk/runs/root_ladder/submit_root.sh
  fi
done

[ "${1:-}" = "go" ] || echo $'\nSTAGED ONLY. Re-run with:  ./stage_wp2.sh go'
