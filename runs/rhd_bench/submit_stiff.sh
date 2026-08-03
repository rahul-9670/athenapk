#!/bin/bash
#SBATCH --job-name=rhd_stiff
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/rhd_bench/%x_%j.out
set -o pipefail
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341 RAD_PRINT_NSUB=1
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/rhd_bench; cd $WD
WRAP=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR -x RAD_PRINT_NSUB"
# stiffer RSLA: creduc=10 (chat = c/10) -> radiation CFL tighter -> should subcycle -> fair stress
for CR in 10 100; do
  echo "=== creduc=$CR 256^3 1 H100 $(date) ==="; rm -rf outc$CR; mkdir outc$CR
  mpirun -n 1 $MCA $WRAP $BIN -i rhd_bench.in radiation/creduc=$CR parthenon/time/nlim=20 -d outc$CR > outc$CR/log 2>&1
  echo "creduc=$CR: steady wsec_step throughput + nsub:"
  grep -E "^cycle=" outc$CR/log | grep -oE "zone-cycles/wsec_step=[0-9.e+-]+" | tail -4
  grep "RAD_NSUB" outc$CR/log | tail -2
done
echo "DONE"
