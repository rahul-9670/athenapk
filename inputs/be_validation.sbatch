#!/bin/bash
#SBATCH --job-name=be_validation
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=be_validation_%j.out
#SBATCH --error=be_validation_%j.err

# Load AthenaPK environment
source ~/athenapk_env.sh

# Change to the run directory
cd /beegfs/u/bbg6470/athenapk_runs/be_validation

# Record start info
echo "===== Job info ====="
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $SLURM_JOB_NODELIST"
echo "Start:        $(date)"
echo "Working dir:  $(pwd)"
echo "===================="

# Launch
mpirun -n ${SLURM_NTASKS} /beegfs/u/bbg6470/athenapk/build/bin/athenaPK \
    -i collapse_be_laptop_test.in

# Record end
echo "===== End ====="
echo "End:          $(date)"
echo "==============="
