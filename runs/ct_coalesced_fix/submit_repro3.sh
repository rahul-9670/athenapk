#!/bin/bash
#SBATCH --job-name=ct_coal_r3
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:50:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/ct_coalesced_fix/%x_%j.out
set -o pipefail
# CT + self-gravity + AMR (multi-level => coarse-fine boundaries => flux-correction path with CT
# Edge flux field). Tests whether the flxcor coalesced path is the coalesced+CT trigger.
# Distinct per-run dirs (off3/on3). 8 ranks. nlim raised so AMR refines + ~40 cycles reached.
source ~/athenapk_env.sh
export OMP_NUM_THREADS=1 PMIX_MCA_gds=hash OMPI_MCA_io=romio341
BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/ct_coalesced_fix
cd $WD; echo "binary:"; md5sum $BIN

echo "=== (1) coalesced OFF (ref) $(date) ==="
rm -rf off3; mkdir -p off3
srun --mpi=pmix -n 8 $BIN -i repro3.in parthenon/mesh/do_coalesced_comms=false -d off3 > off3/log 2>&1
echo "OFF exit=$? cycles=$(grep -cE '^cycle=' off3/log) driver=$(grep -c 'Driver completed' off3/log) nblk_max=$(grep -oE 'nbtotal=[0-9]+' off3/log | tail -1)"

echo "=== (2) coalesced ON (reproduce) $(date) ==="
rm -rf on3; mkdir -p on3
srun --mpi=pmix -n 8 $BIN -i repro3.in parthenon/mesh/do_coalesced_comms=true -d on3 > on3/log 2>&1
echo "ON exit=$? cycles=$(grep -cE '^cycle=' on3/log) recv-buffer-err=$(grep -c 'receiving buffer' on3/log) driver=$(grep -c 'Driver completed' on3/log)"
echo "=== DONE $(date) ==="
