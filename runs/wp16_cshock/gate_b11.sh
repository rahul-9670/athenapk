#!/bin/bash
#SBATCH --job-name=gateb11
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
# OFF-STATE GATE for the B11 Hall-floor startup notice (build_cpu d7d28f11 -> 2ddea223).
# The change is a GetOrAddReal + a std::cout, and the key it reads was already registered by the
# branch above it -- so it should be provably inert. "Should be" is not a gate, so: rerun the
# exact deck that d7d28f11 ran in stock_hall_whistler_glm/ and byte-compare the history.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads
H=/beegfs/u/bbg6470/athenapk/runs/wp16_cshock
G=$H/gate_b11_new; rm -rf $G; mkdir -p $G
echo "new binary: $(md5sum /beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK)"
( cd $G && /beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK \
    -i /beegfs/u/bbg6470/athenapk/inputs/hall_whistler_glm.in > run.log 2>&1 )
echo "run exit=$?"
OLD=$(ls $H/stock_hall_whistler_glm/*.hst 2>/dev/null | head -1)
NEW=$(ls $G/*.hst 2>/dev/null | head -1)
echo "old(d7d28f11)=$OLD"; echo "new(2ddea223)=$NEW"
if [ -n "$OLD" ] && [ -n "$NEW" ] && cmp -s "$OLD" "$NEW"; then
  echo "PASS: history BYTE-IDENTICAL -> the B11 notice is an OFF-state no-op."
else
  echo "FAIL: history differs (or a file is missing)."; diff "$OLD" "$NEW" | head -10
fi
echo "notice present in the new run: $(grep -c 'NOTE \[Hall\]' $G/run.log)"
