#!/bin/bash
#SBATCH --job-name=ws1_gate4
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --output=%x_%j.out
set -o pipefail
source ~/athenapk_env.sh
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
SUB=/beegfs/u/bbg6470/athenapk/runs/validation_ws1/submit_gate4.sh
WDIR=/beegfs/u/bbg6470/athenapk/runs/validation_ws1/gate4_run
mkdir -p "$WDIR"; cd "$WDIR"
if grep -q "Driver completed" run.log 2>/dev/null; then echo "done (tlim)"; exit 0; fi
N=$(cat chain_n 2>/dev/null || echo 0); N=$((N+1)); echo $N > chain_n
if [ "$N" -lt 6 ]; then sbatch --dependency=afterany:$SLURM_JOB_ID "$SUB" && echo "successor queued"; fi
LATEST=$(ls -t "$WDIR"/parthenon.out2.*.rhdf 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then RA="-r $LATEST"; echo "RESUME $LATEST"; else RA="-i ../dt_recovery_fixed.in"; echo "FRESH"; fi
stdbuf -oL -eL srun --mpi=pmix -n 32 "$BIN" $RA >> run.log 2>&1
echo "EXIT $? $(date)"; grep -iE "created sink" run.log | tail -3
