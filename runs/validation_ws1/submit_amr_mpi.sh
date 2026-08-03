#!/bin/bash
#SBATCH --job-name=ws1_sinkamr
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:40:00
#SBATCH --output=%x_%j.out
# WS-1 increment 1: inert sink through dynamic AMR across 32 MPI ranks; then a restart
# leg for the bit-identical check. CPU/std.
set -o pipefail
source ~/athenapk_env.sh
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
WDIR=/beegfs/u/bbg6470/athenapk/runs/validation_ws1/amr_mpi
rm -rf "$WDIR" && mkdir -p "$WDIR" && cd "$WDIR"

echo "=== FRESH run to t=5 (32 ranks, adaptive) ==="
stdbuf -oL -eL srun --mpi=pmix -n 32 "$BIN" -i ../sink_amr.in >run.log 2>&1
echo "fresh exit=$?"; grep -iE "Driver completed|error|abort|nan" run.log | tail -3

# Restart leg: resume from the ~t=2.5 restart into a subdir, run to t=5, for bit-identical check.
RST=$(ls -t "$WDIR"/parthenon.out1.*.rhdf 2>/dev/null | grep -v final | sort | sed -n '2p')
[ -z "$RST" ] && RST=$(ls "$WDIR"/parthenon.out1.00001.rhdf 2>/dev/null)
echo "=== RESTART leg from $RST ==="
mkdir -p "$WDIR/restart" && cd "$WDIR/restart"
stdbuf -oL -eL srun --mpi=pmix -n 32 "$BIN" -r "$RST" >rr.log 2>&1
echo "restart exit=$?"; grep -iE "Driver completed|error|abort|nan" rr.log | tail -3
echo "ALL_DONE"
