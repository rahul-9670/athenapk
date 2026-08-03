#!/bin/bash
#SBATCH --job-name=hydro_bench
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:12:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/rhd_bench/%x_%j.out
set -o pipefail
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/rhd_bench; cd $WD
WRAP=/beegfs/u/bbg6470/athenapk/runs/flagship_integration/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
echo "=== 256^3 PURE HYDRO (radiation off), 1 H100 $(date) ==="; rm -rf outh; mkdir outh
mpirun -n 1 $MCA $WRAP $BIN -i rhd_bench.in physics/radiation=false parthenon/time/nlim=60 -d outh > outh/log 2>&1
grep -E "^cycle=" outh/log | grep -oE "cycle=[0-9]+ .*wsec_step=[0-9.e+-]+ zone-cycles/wsec_step=[0-9.e+-]+" | tail -6
grep -E "zone-cycles/wallsecond|walltime" outh/log
echo "DONE"
