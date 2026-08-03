#!/bin/bash
cd /beegfs/u/bbg6470/athenapk/runs/flagship_integration
BUILD=2414619; BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
before=a040a8d71f0b5a7dee7dac42d8005bcb
# 1) wait for build
for i in $(seq 1 240); do
  [ -z "$(squeue -j $BUILD -h -o %T 2>/dev/null)" ] && break; sleep 15
done
md=$(md5sum "$BIN" 2>/dev/null | awk '{print $1}')
if [ "$md" = "$before" ] || [ -z "$md" ]; then echo "BUILD FAILED/unchanged (md5=$md)"; exit 0; fi
echo "BUILD OK new binary md5=$md; launching flagship Boris test"
# 2) submit flagship boris test
rm -rf probe_boris; JID=$(sbatch --parsable probe_boris.sh)
echo "probe_boris job=$JID"
HST=probe_boris/parthenon.out0.hst; EVENT=1.10354; PAST=1.10362
# 3) watch verdict
for i in $(seq 1 400); do
  st=$(squeue -j $JID -h -o %T 2>/dev/null)
  crash=$(grep -ciE "illegal|cudaError|Aborted|signal" probe_boris/run.log 2>/dev/null)
  if [ "${crash:-0}" -gt 0 ]; then echo "VERDICT: CRASH - $(grep -iE 'illegal|cudaError|Aborted' probe_boris/run.log | head -1)"; exit 0; fi
  if [ -f "$HST" ]; then
    nan=$(awk -v e=$EVENT 'NR>2 && $1>=e && ($6=="-nan"){print $1; exit}' "$HST")
    [ -n "$nan" ] && { echo "VERDICT: BORIS FAILED - NaN at t=$nan (still blows at event)"; exit 0; }
    tmax=$(awk 'NR>2 && $6!="-nan"{t=$1} END{print t}' "$HST")
    ok=$(awk -v t="$tmax" -v p=$PAST 'BEGIN{print (t>=p)?1:0}')
    [ "$ok" = "1" ] && { echo "VERDICT: BORIS WORKS - clean past the event to t=$tmax (no NaN); cycles/dt in run.log"; exit 0; }
  fi
  [ -z "$st" ] && { echo "VERDICT: job ended; last hst=[$(awk 'NR>2{l=$0}END{print l}' $HST 2>/dev/null)] $(grep RUN_EXIT probe_boris/run.log|tail -1)"; exit 0; }
  sleep 20
done
echo "watcher timeout"
