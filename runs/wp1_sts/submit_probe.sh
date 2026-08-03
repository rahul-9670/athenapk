#!/bin/bash
#SBATCH --job-name=wp1sts
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:h100:4
#SBATCH --cpus-per-task=8
#SBATCH --time=00:40:00
#SBATCH --exclusive
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp1_sts/slurm_%x_%j.out
set -o pipefail
source ~/athenapk_env.sh
module load cuda/12.5.1
R=/beegfs/u/bbg6470/athenapk/runs/wp1_sts/$LEG
B=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
W=/beegfs/u/bbg6470/athenapk/runs/wp1_sts/wrap_alloc.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
export OMP_NUM_THREADS=1 PMIX_MCA_gds=hash OMPI_MCA_io=romio341
echo "LEG=$LEG OVERRIDES=$OV  BINARY: $(md5sum $B)"
# self-resuming: pick the newest restart if one exists
LATEST=$(ls -1t $R/*.rhdf 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then RA="-r $LATEST"; echo "resuming from $LATEST"; else RA="-i $R/fhc.in"; fi
cd $R && stdbuf -oL -eL mpirun -n 4 $MCA $W $B $RA $OV > $R/run_$SLURM_JOB_ID.log 2>&1
echo "exit=$?" >> $R/status
