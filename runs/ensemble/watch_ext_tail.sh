#!/bin/bash
# Watch the tail of the 2e-12 extension campaign: the members that have not yet
# written STOP_CHAIN. Emits one line per state CHANGE, plus a heartbeat, and exits
# when every watched member has stopped.
#
# Chain-proof: a chained successor inherits `#SBATCH --job-name=ens`, so members are
# resolved by WorkDir, never by job name. A member counts as in flight if ANY of its
# squeue rows exists; RUNNING if any row is RUNNING.
D=/beegfs/u/bbg6470/athenapk/runs/ensemble/design01
WATCH="${WATCH:-004 015 016}"
# NOISE. The change key deliberately EXCLUDES the snapshot count: a member writes a snapshot
# every ~8-10 min, so keying on it turned every ordinary write into a notification. Only the
# per-member STATE (STOP/RUN/PEND/GONE) is a change worth waking for; progress rides along on a
# half-hourly heartbeat instead.
prev=""; beat=0
while true; do
  q=$(squeue -u bbg6470 -o "%T|%Z" -h 2>/dev/null)
  cur=""; prog=""; done_n=0; total=0; stalled=""
  for p in $WATCH; do
    total=$((total+1)); d=$D/point$p
    if [ -f "$d/STOP_CHAIN" ]; then
      st=STOP; done_n=$((done_n+1))
    else
      rows=$(echo "$q" | grep -F "point$p" || true)
      if [ -z "$rows" ]; then st=GONE; stalled="$stalled point$p";
      elif echo "$rows" | grep -q RUNNING; then st=RUN; else st=PEND; fi
    fi
    cur="$cur point$p:$st"
    prog="$prog point$p:$st:$(ls "$d"/parthenon.out1.*.phdf 2>/dev/null | wc -l)"
  done
  if [ "$cur" != "$prev" ] || [ $((beat % 10)) -eq 0 ]; then
    echo "$(date +%H:%M) $prog  [done $done_n/$total quota $(rrz-quota 2>/dev/null | awk '/beegfs/&&/GiB/{print $2$3}')]"
    prev="$cur"
  fi
  [ -n "$stalled" ] && echo "$(date +%H:%M) ALERT no job in flight and no STOP_CHAIN:$stalled"
  [ "$done_n" -eq "$total" ] && { echo "$(date +%H:%M) ALL $total EXTENSION MEMBERS STOPPED"; exit 0; }
  beat=$((beat+1)); sleep 180
done
