#!/bin/bash
#SBATCH --job-name=mg_bench
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
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader|head -1)"
# multigroup M1: n_group=3, constant opacity (isolate multigroup cost from tabulated-opacity cost)
echo "=== n_group=3 (multigroup) 256^3 1 H100 $(date) ==="; rm -rf outmg; mkdir outmg
mpirun -n 1 $MCA $WRAP $BIN -i rhd_bench.in radiation/n_group=3 parthenon/time/nlim=30 -d outmg > outmg/log 2>&1
echo "--- n_group=3 throughput + nsub ---"
grep -E "^cycle=" outmg/log | grep -oE "zone-cycles/wsec_step=[0-9.e+-]+" | tail -5
grep RAD_NSUB outmg/log | tail -2
echo "DONE"
