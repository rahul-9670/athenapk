#!/bin/bash
# THE SINGLE AUTHORISED RE-BASELINE: root ladder on the new physics.
#
#   binary  49d9c257 -> f181c0a1  (build_gpu_v2; OFF-state gate PASSED byte-identical on this
#                                  exact deck, job 2448500 -- so the binary alone changes nothing)
#   fluid BCs   outflow -> diode          (B1)
#   gravity     absolute -> relative      (B3)
#
# Both switches are result-changing and both are now IN THE DECK (fhc_rootladder.in), so this is
# the reference every later WP is measured against. Legs are named *_sw to keep them distinct from
# r*_gfix (the 49d9c257 + outflow + absolute baseline that Gate A validated).
#
# Rank counts follow WP-12 (decomposition-invariant) and the 74 %-efficiency measurement: 8 ranks
# BURNS GPU-hours, so use the smallest count that fits. 512^3 needs MBLK=64 to fit 80 GiB cards.
#
# Usage:  ./stage_switches.sh          -> print what would be submitted
#         ./stage_switches.sh go       -> submit
#         ./stage_switches.sh go 512   -> submit only the 512 rung
set -o pipefail
HERE=/beegfs/u/bbg6470/athenapk/runs/root_ladder
WRAPSRC=/beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh
BIN=/beegfs/u/bbg6470/athenapk/build_gpu_v2/bin/athenaPK
[ -r "$WRAPSRC" ] || { echo "FATAL: no GPU-pinning wrapper at $WRAPSRC"; exit 1; }
[ -x "$BIN" ]     || { echo "FATAL: no binary at $BIN"; exit 1; }
echo "binary: $(md5sum $BIN)"
echo "deck:   $(md5sum $HERE/fhc_rootladder.in)"

#      NX   DIR        NRANK  MBLK
# 2026-08-03: 256 rung 6 -> 3 ranks. Not a physics change (WP-12: decomposition-invariant) and
# not a cost increase (4->8 scales at 74 %, so FEWER ranks is cheaper per GPU-hour). It is a
# SCHEDULING fix: `MaxNodes=1` means a job cannot span nodes, and with g001 full and g002/g003/g004
# holding 3/3/1 free GPUs, a 6-GPU request was unschedulable indefinitely. Memory is not the
# constraint -- r256_gfix measured 13.5 GiB/GPU on 6 ranks (81.0 GiB total), so 3 ranks is ~27
# GiB/GPU against an 80 GiB card. 512 stays at 8: with MBLK=64 it needs ~69.8 GiB/GPU there.
LEGS="128  r128_sw     3      -
      256  r256_sw     3      -
      512  r512_sw     8      64"

ONLY=${2:-}
while read -r NX DIR NRANK MBLK; do
  [ -z "$NX" ] && continue
  [ -n "$ONLY" ] && [ "$ONLY" != "$NX" ] && continue
  mkdir -p $HERE/$DIR
  install -m 755 $WRAPSRC $HERE/$DIR/wrap_mod.sh || { echo "FATAL: cannot install wrapper"; exit 1; }
  CMD="sbatch --job-name=root${NX}sw --nodes=1 --ntasks=$NRANK --gres=gpu:h100:$NRANK \
--export=ALL,NX=$NX,RUNDIR=$HERE/$DIR,NRANK=$NRANK,BIN=$BIN"
  [ "$MBLK" != "-" ] && CMD="$CMD,MBLK=$MBLK"
  CMD="$CMD $HERE/submit_root.sh"
  if [ "${1:-}" = "go" ]; then echo "+ $CMD"; eval $CMD; else echo "would run: $CMD"; fi
done <<< "$LEGS"
