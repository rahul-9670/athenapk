#!/bin/bash
# Watch the two validation jobs launched 2026-08-08 and report their VERDICTS, not just their exit
# states. Emits one line per job state change, then the analysis verdict for each, then exits.
#
# Coverage: the terminal-state grep covers COMPLETED, FAILED, TIMEOUT, CANCELLED and NODE_FAIL, so
# a crash is as loud as a success. Silence here means "still queued or running", nothing else.
J13=${J13:-2490880}   # WP-13b GPU+AMR restart reproducibility
J17=${J17:-2490875}   # WS-1 inc6 GPU sink accretion
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
R=/beegfs/u/bbg6470/athenapk/runs
declare -A done_j
prev=""
while true; do
  line=""
  for pair in "WP13b:$J13" "WP17gpu:$J17"; do
    name=${pair%%:*}; jid=${pair##*:}
    st=$(squeue -j $jid -h -o %T 2>/dev/null)
    [ -z "$st" ] && st=$(sacct -j $jid -n -X -o State 2>/dev/null | head -1 | tr -d ' ')
    [ -z "$st" ] && st=UNKNOWN
    line="$line $name:$st"
    case "$st" in
      COMPLETED|FAILED|TIMEOUT|CANCELLED*|NODE_FAIL|OUT_OF_MEMORY)
        if [ -z "${done_j[$name]}" ]; then
          done_j[$name]=1
          echo "$(date +%H:%M) $name job $jid finished: $st"
          if [ "$name" = "WP13b" ]; then
            echo "--- WP-13b status file ---"; cat $R/wp13b_gpu_restart/status 2>/dev/null
            echo "--- WP-13b verdict ---"
            $PY $R/wp13b_gpu_restart/compare_wp13b.py 2>/dev/null | grep -E "VERDICT|floor|divergence|PASS|FAIL|BLOCK COUNT|cycle  :"
          else
            echo "--- WS-1 inc6 verdict ---"
            grep -E "sink mass|NO ACCRETION|ABSENT|plateau|VERDICT|PASS|FAIL|exit=" $R/wp17_sinks/wp17gpu_${jid}.out 2>/dev/null | tail -20
          fi
        fi ;;
    esac
  done
  [ "$line" != "$prev" ] && { echo "$(date +%H:%M)$line"; prev="$line"; }
  [ -n "${done_j[WP13b]}" ] && [ -n "${done_j[WP17gpu]}" ] && { echo "$(date +%H:%M) both validation jobs done"; exit 0; }
  sleep 120
done
