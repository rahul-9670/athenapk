#!/bin/bash
#SBATCH --job-name=wp13gpu
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp13_gpu/slurm_%j.out
set -o pipefail
source ~/athenapk_env.sh
module load cuda/12.5.1
R=/beegfs/u/bbg6470/athenapk/runs/wp13_gpu
B=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
W=$R/wrap_alloc.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
export OMP_NUM_THREADS=1 PMIX_MCA_gds=hash OMPI_MCA_io=romio341
echo "BINARY md5: $(md5sum $B)"

# Leg A: 18 cycles straight through.
cd $R/straight && stdbuf -oL -eL mpirun -n 2 $MCA $W $B -i $R/straight/fhc.in \
  parthenon/time/nlim=18 parthenon/output2/dn=6 > $R/straight/run.log 2>&1
echo "straight exit=$?" >> $R/status

# Leg B: 6 cycles, stop, restart from the cycle-6 dump, finish to 18.
cd $R/split && stdbuf -oL -eL mpirun -n 2 $MCA $W $B -i $R/split/fhc.in \
  parthenon/time/nlim=6 parthenon/output2/dn=6 > $R/split/run1.log 2>&1
echo "split1 exit=$?" >> $R/status
RST=$(ls -1t $R/split/*.rhdf 2>/dev/null | head -1)
cp "$RST" $R/rst.rhdf && echo "restart_from=$RST" >> $R/status
cd $R/split && stdbuf -oL -eL mpirun -n 2 $MCA $W $B -r $R/rst.rhdf \
  parthenon/time/nlim=18 > $R/split/run2.log 2>&1
echo "split2 exit=$?" >> $R/status
echo DONE >> $R/status
