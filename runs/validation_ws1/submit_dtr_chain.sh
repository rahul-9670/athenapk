#!/bin/bash
#SBATCH --job-name=ws1_dtr
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
WDIR=/beegfs/u/bbg6470/athenapk/runs/validation_ws1/dtr_run
SUB=/beegfs/u/bbg6470/athenapk/runs/validation_ws1/submit_dtr_chain.sh
mkdir -p "$WDIR"; cd "$WDIR"
# stop conditions: sink already formed (goal reached) or tlim done
if grep -q "Driver completed" run.log 2>/dev/null; then echo "done (tlim)"; exit 0; fi
if grep -q "created sink" run.log 2>/dev/null; then echo "sink formed -> goal reached, stop chain"; exit 0; fi
N=$(cat chain_n 2>/dev/null || echo 0); N=$((N+1)); echo $N > chain_n
if [ "$N" -lt 8 ]; then sbatch --dependency=afterany:$SLURM_JOB_ID "$SUB" && echo "successor queued"; fi
LATEST=$(ls -t "$WDIR"/parthenon.out2.*.rhdf 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then RA="-r $LATEST"; echo "RESUME $LATEST"; else RA="-i ../dt_recovery.in"; echo "FRESH"; fi
stdbuf -oL -eL srun --mpi=pmix -n 32 "$BIN" $RA >> run.log 2>&1
echo "EXIT $? $(date)"; grep -iE "created sink" run.log | tail -2
