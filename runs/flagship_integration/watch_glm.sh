#!/bin/bash
cd /beegfs/u/bbg6470/athenapk/runs/flagship_integration
JOB=2414820; HST=probe_glm/parthenon.out0.hst; EVENT=1.10354; PAST=1.10365
for i in $(seq 1 400); do
  st=$(squeue -j $JOB -h -o %T 2>/dev/null)
  crash=$(grep -ciE "illegal|cudaError|Aborted|FATAL|signal 7|Bus" probe_glm/run.log 2>/dev/null)
  [ "${crash:-0}" -gt 0 ] && { echo "VERDICT: GLM CRASH/incompat - $(grep -iE 'illegal|cudaError|Aborted|FATAL|Message:' probe_glm/run.log | head -2)"; exit 0; }
  if [ -f "$HST" ]; then
    nan=$(awk -v e=$EVENT 'NR>2 && $1>=e && ($6=="-nan"){print $1; exit}' "$HST")
    [ -n "$nan" ] && { echo "VERDICT: GLM ALSO NaN at t=$nan (wall baked into cyc1000 CT state; need GLM-from-earlier)"; exit 0; }
    tmax=$(awk 'NR>2 && $6!="-nan"{t=$1;dt=$2} END{print t" "dt}' "$HST")
    tt=$(echo $tmax|awk '{print $1}')
    ok=$(awk -v t="$tt" -v p=$PAST 'BEGIN{print (t>=p)?1:0}')
    [ "$ok" = "1" ] && { echo "VERDICT: GLM WORKS - clean past the event to t=$tmax (dt shown); CT was the cause"; exit 0; }
  fi
  [ -z "$st" ] && { echo "VERDICT: job ended; last hst=[$(awk 'NR>2{l=$0}END{print l}' $HST 2>/dev/null)]; $(grep RUN_EXIT probe_glm/run.log|tail -1)"; exit 0; }
  sleep 20
done
echo "watcher timeout"
