#!/bin/bash
#SBATCH --job-name=ws1_2body
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:40:00
#SBATCH --output=%x_%j.out
set -o pipefail
source ~/athenapk_env.sh
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
WDIR=/beegfs/u/bbg6470/athenapk/runs/validation_ws1/twobody_run
rm -rf "$WDIR" && mkdir -p "$WDIR" && cd "$WDIR"
stdbuf -oL -eL srun --mpi=pmix -n 8 "$BIN" -i ../twobody.in >run.log 2>&1
echo "EXIT $? $(date)"; grep -iE "Driver completed|error|abort|nan" run.log | tail -2
