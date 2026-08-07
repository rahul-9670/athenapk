#!/bin/bash
#SBATCH -J ens_uq
#SBATCH -p std
#SBATCH -A banerjee_std
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -o /beegfs/u/bbg6470/athenapk/runs/ensemble/uq_both_epochs.slurm.%j.log
#SBATCH -e /beegfs/u/bbg6470/athenapk/runs/ensemble/uq_both_epochs.slurm.%j.log

# Ensemble UQ at both candidate measurement epochs.
#
# Serial and almost entirely IO-bound: ~197 snapshot reads at ~24 s each (3 s of CPU per 8 min of
# wall), so this asks for 1 task and a short walltime rather than a node. It lives in the batch
# queue rather than on the front-end because the run is ~2 h and a session-parented front-end
# process is reapable -- two earlier attempts were killed that way before writing a single line.
#
# Restart-safe by construction: uq_both_epochs.py caches rho_max per snapshot keyed by
# (path, size, mtime_ns) and flushes every 10 entries, so a requeue resumes from the cache and
# re-reads at most the last few snapshots. Re-running this script is therefore always safe.

set -euo pipefail
source /home/bbg6470/athenapk_env.sh

PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
cd /beegfs/u/bbg6470/athenapk/runs/ensemble

echo "host=$(hostname)  job=${SLURM_JOB_ID:-none}  start=$(date -Is)"
$PY -u uq_both_epochs.py design01
echo "end=$(date -Is)  rc=$?"
