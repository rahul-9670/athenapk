#!/bin/bash
#SBATCH --job-name=wp22eta
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
# WP-22 part 2 -- reconstruct the PHYSICAL eta_A from a deep production snapshot.
# Runs on SLURM, not the front-end: the front-end is a 1-CPU cgroup and each prim component of
# prod_v9's 2402-block snapshot is ~630 MB off BeeGFS.
set -o pipefail
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
S=/beegfs/u/bbg6470/athenapk/docs/validation/scripts/wp22_eta_physical.py
echo "job $SLURM_JOB_ID $(date)"
for F in /beegfs/u/bbg6470/athenapk/runs/prod_v9/parthenon.out1.00096.phdf \
         /beegfs/u/bbg6470/athenapk/runs/prod_t4_full/parthenon.out1.01713.phdf; do
  echo "================================================================"
  $PY $S $F
done
echo "done $(date)"
