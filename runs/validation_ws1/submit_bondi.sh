#!/bin/bash
#SBATCH --job-name=ws1_bondi
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
set -o pipefail
source ~/athenapk_env.sh
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
WDIR=/beegfs/u/bbg6470/athenapk/runs/validation_ws1/bondi_run
rm -rf "$WDIR" && mkdir -p "$WDIR" && cd "$WDIR"
stdbuf -oL -eL srun --mpi=pmix -n 32 "$BIN" -i ../bondi.in >run.log 2>&1
echo "EXIT $? $(date)"; grep -iE "Driver completed|error|abort|nan" run.log | tail -3
