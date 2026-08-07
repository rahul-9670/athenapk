#!/bin/bash
# Extend the ensemble members that stopped short of the 1e-12 measurement window.
#
# WHY. The 2026-08-07 two-epoch analysis (job 2478818) kept 10 of 24 members at 1e-13 and only
# 5 of 24 at 1e-12. At 1e-12 the exclusions are NOT a measurement problem: 17 members simply never
# reached the acceptance window (rho in [5.01e-13, 2.00e-12], i.e. --max-log-dist 0.3 around
# 1e-12), because the first campaign stopped every member at STOP_CGS=1e-13. Those 17 are listed
# below and are extended here to STOP_CGS=2.0e-12, which is where point000 -- the member that made
# the trip unaided -- actually ended.
#
# CADENCE. Passing OUT_DN=25 (deck default is 50 cycles). Measured on the completed campaign:
#   * BEFORE first core, dn=50 spans 1.3-1.9 decades of rho per snapshot -- far too coarse;
#   * AFTER first core, dn=50 already spans only 0.08-0.17 decades (point000, 00007 -> final).
# Every member extended here restarts at >= 1.1x rho_crit, i.e. ALREADY PAST the transition, so
# its remaining leg sits entirely in the well-resolved regime and dn=50 would have sufficed.
# dn=25 (~0.04-0.09 dec) is bought as margin in case a member's dt recovers faster than point000's.
#
# STORAGE, checked before launch rather than after. Snapshots near this epoch are ~4.07 GiB
# (point000, out1.00011+). point000's 1.43e-13 -> 2.00e-12 leg was 9 snapshots at dn=50 = ~450
# cycles, so at dn=25 expect ~18 snapshots = ~73 GiB per member, ~1.25 TiB over 17 members.
# Quota headroom at launch was 4536 GiB of 10998 GiB (59% used), so this is ~27% of what is free.
# dn=5 was REJECTED for exactly this reason: ~244 GiB/member = ~4.15 TiB would have consumed
# essentially all remaining quota.
#
# WATCHDOG. watchdog_epoch_stop.py MUST be run with the SAME STOP_CGS as this launch. It stops a
# member once its newest snapshot is at or past STOP_CGS; left at the old 1e-13 it would terminate
# every one of these members immediately, on their very first poll, before they advanced at all.
#
# Usage: ./extend_short_members.sh [STOP_CGS] [OUT_DN]
set -euo pipefail
STOP_CGS="${1:-2.0e-12}"
OUT_DN="${2:-25}"
HERE=/beegfs/u/bbg6470/athenapk/runs/ensemble
DDIR=$HERE/design01

# The 17 members excluded at 1e-12 for never reaching the window (job 2478818 log).
# point021/point023 are deliberately NOT here: they have ample margin (31.3x, 26.4x) and are
# already PAST 1e-12, so extending them achieves nothing -- they were excluded for overshooting
# the window by 0.50 and 0.42 decades. Recovering those two needs a re-run of their final leg at
# fine cadence from an earlier restart, which is a different operation.
MEMBERS="003 004 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020"

echo "extending to STOP_CGS=$STOP_CGS with output cadence dn=$OUT_DN"
for p in $MEMBERS; do
  pdir=$DDIR/point$p
  [ -f "$pdir/fhc_ens.in" ] || { echo "  point$p: NO DECK -- skipped"; continue; }
  newest=$(ls -t "$pdir"/parthenon.out2.*.rhdf 2>/dev/null | head -1)
  [ -n "$newest" ] || { echo "  point$p: NO RESTART -- skipped"; continue; }

  # Archive rather than delete: the marker records where epoch 1 stopped, which is evidence.
  # submit_point.sh exits immediately if STOP_CHAIN exists, so it has to be moved aside.
  [ -f "$pdir/STOP_CHAIN" ] && mv "$pdir/STOP_CHAIN" "$pdir/STOP_CHAIN.epoch1e13"

  jid=$(sbatch --job-name="ext_point$p" \
        --export=ALL,PDIR="$pdir",STOP_CGS="$STOP_CGS",OUT_DN="$OUT_DN" \
        "$HERE/submit_point.sh" | awk '{print $NF}')
  echo "  point$p -> job $jid  (restart $(basename "$newest"))"
done
echo
echo "NOW START THE WATCHDOG WITH THE MATCHING THRESHOLD:"
echo "  STOP_CGS=$STOP_CGS nohup setsid \$PY $HERE/watchdog_epoch_stop.py &"
