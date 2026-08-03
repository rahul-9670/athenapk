#!/bin/bash
#SBATCH --job-name=ct_coal_repro
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:40:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/ct_coalesced_fix/%x_%j.out
set -o pipefail
# Coalesced-comms + CT reproducer + correctness-gate baseline, on a COMPUTE NODE (fast, no
# front-end contention). Runs the minimal CT+self-gravity config 3 ways:
#   (1) coalesced OFF  -> gate REFERENCE (must complete; save phdf)
#   (2) coalesced ON   -> REPRODUCE the crash (pre-fix: aborts; post-fix: completes + matches ref)
# 8 MPI ranks so there is real inter-rank comm (the "receiving buffer" path).
source ~/athenapk_env.sh
export OMP_NUM_THREADS=1 PMIX_MCA_gds=hash OMPI_MCA_io=romio341
BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/ct_coalesced_fix
cd $WD
echo "binary:"; md5sum $BIN

echo "=== (1) coalesced OFF (gate reference) $(date) ==="
rm -rf off; mkdir -p off
srun --mpi=pmix -n 8 $BIN -i repro.in parthenon/mesh/do_coalesced_comms=false \
  -d off > off/log 2>&1
echo "OFF exit=$? cycles=$(grep -cE '^cycle=' off/log) driver=$(grep -c 'Driver completed' off/log)"

echo "=== (2) coalesced ON (reproduce / post-fix check) $(date) ==="
rm -rf on; mkdir -p on
srun --mpi=pmix -n 8 $BIN -i repro.in parthenon/mesh/do_coalesced_comms=true \
  -d on > on/log 2>&1
echo "ON exit=$? cycles=$(grep -cE '^cycle=' on/log) recv-buffer-err=$(grep -c 'receiving buffer' on/log) driver=$(grep -c 'Driver completed' on/log)"
echo "=== DONE $(date) ==="
