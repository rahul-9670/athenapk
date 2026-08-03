#!/bin/bash
#SBATCH --job-name=rt_stream
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --output=%x_%j.out
# WS-3b gate: free-streaming M1 pulse, donor-cell (dc) vs PLM. Same nlim => same final time;
# compare the FWHM of the advected Er pulse. PLM must cut the FWHM growth by >=2x vs dc.
set -o pipefail
source ~/athenapk_env.sh
export OMP_NUM_THREADS=8 OMPI_MCA_io=romio341
BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
cd /beegfs/u/bbg6470/athenapk/runs/validation_rt
for R in dc plm; do
  d=s_$R; rm -rf "$d"; mkdir -p "$d"; cd "$d"
  echo "=== reconstruction=$R ==="
  $BIN -i ../stream.in radiation/reconstruction=$R >run.log 2>&1
  echo "recon=$R exit $? cycles=$(grep -cE 'cycle=' run.log)"
  cd ..
done
echo "DONE $(date)"
