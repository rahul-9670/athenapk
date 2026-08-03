#!/bin/bash
cd /beegfs/u/bbg6470/athenapk/runs/flagship_integration
JOB=2414468; HST=probe_fofc/parthenon.out0.hst; EVENT=1.10354; PAST=1.10362
for i in $(seq 1 480); do   # up to ~2h at 15s
  # job gone?
  st=$(squeue -j $JOB -h -o "%T" 2>/dev/null)
  if [ -f "$HST" ]; then
    # NaN past the event?
    nanrow=$(awk -v e=$EVENT 'NR>2 && $1>=e && ($6=="-nan"){print $1; exit}' "$HST")
    if [ -n "$nanrow" ]; then echo "VERDICT: FOFC FAILED - NaN at t=$nanrow (still blows at the event)"; exit 0; fi
    # crossed past the event clean?
    tmax=$(awk 'NR>2 && $6!="-nan"{t=$1} END{print t}' "$HST")
    ok=$(awk -v t="$tmax" -v p=$PAST 'BEGIN{print (t>=p)?1:0}')
    if [ "$ok" = "1" ]; then echo "VERDICT: FOFC WORKS - clean past the event, reached t=$tmax (no NaN through t>=$PAST)"; exit 0; fi
  fi
  if [ -z "$st" ]; then
    # job ended; report final state
    fn=$(awk 'NR>2{last=$0} END{print last}' "$HST" 2>/dev/null)
    exitline=$(grep RUN_EXIT probe_fofc/run.log 2>/dev/null | tail -1)
    echo "VERDICT: job ended. last hst=[$fn] $exitline"; exit 0
  fi
  sleep 15
done
echo "VERDICT: watcher timed out after ~2h; job state=$st"
