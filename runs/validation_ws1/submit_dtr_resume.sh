#!/bin/bash
#SBATCH --job-name=ws1_dtrR
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:25:00
#SBATCH --output=%x_%j.out
set -o pipefail
source ~/athenapk_env.sh
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
WDIR=/beegfs/u/bbg6470/athenapk/runs/validation_ws1/dtr_resume
rm -rf "$WDIR" && mkdir -p "$WDIR" && cd "$WDIR"
# resume the deep-collapse state with a FIXED modest rho_sink so a sink forms NOW and
# accretion caps the core -> dt should recover. Also cap tlim so it stops after recovery.
stdbuf -oL -eL srun --mpi=pmix -n 32 "$BIN" -r /beegfs/u/bbg6470/athenapk/runs/validation_ws1/dtr_run/parthenon.out2.00016.rhdf \
  sinks/rho_sink_code=4000.0 parthenon/time/tlim=1.5 parthenon/output2/dn=100000 >run.log 2>&1
echo "EXIT $? $(date)"; grep -iE "created sink|Driver completed" run.log | head -3
