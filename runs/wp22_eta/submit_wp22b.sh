#!/bin/bash
#SBATCH --job-name=wp22b
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
set -o pipefail
/beegfs/u/bbg6470/venvs/analysis_env/bin/python \
  /beegfs/u/bbg6470/athenapk/docs/validation/scripts/wp22_eta_applied.py \
  /beegfs/u/bbg6470/athenapk/runs/wp22_eta/dump/parthenon.out1.00096.phdf
echo "done $(date)"
