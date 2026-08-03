#!/bin/bash
JOB=2414464; BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
before=750dacf7c17844de9a4f3b175345271b
for i in $(seq 1 240); do   # up to ~1h
  st=$(squeue -j $JOB -h -o "%T" 2>/dev/null)
  if [ -z "$st" ]; then
    md=$(md5sum "$BIN" 2>/dev/null | awk '{print $1}')
    out=$(ls -t /beegfs/u/bbg6470/athenapk/runs/*build_gpu*.out /beegfs/u/bbg6470/athenapk/runs/*.out 2>/dev/null | head -1)
    err=$(grep -ciE "error:|Error 1|undefined reference" $out 2>/dev/null)
    if [ "$md" != "$before" ] && [ -n "$md" ]; then echo "BUILD DONE - new binary md5=$md (changed from $before), build errors=$err"; 
    else echo "BUILD JOB ENDED but binary md5=$md unchanged/missing (build may have FAILED); errors=$err; log=$out"; fi
    exit 0
  fi
  sleep 15
done
echo "build watcher timeout; job state=$st"
