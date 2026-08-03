#!/bin/bash
#SBATCH --job-name=rt_equil
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:20:00
#SBATCH --output=%x_%j.out
# WS-3a gate (i): Planck/Rosseland split equilibration + rate-scaling test.
# Uniform gas + uniform radiation (0-D matter coupling). Weak per-step coupling (creduc=30)
# so the relaxation is a resolved exponential; ratio=2 must relax ~2x faster than ratio=1
# (tau ~ 1/kappa_P) and both reach the SAME equilibrium (opacity-independent).
set -o pipefail
source ~/athenapk_env.sh
export OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
cd /beegfs/u/bbg6470/athenapk/runs/validation_rt
OV="parthenon/mesh/nx1=8 parthenon/mesh/nx2=8 parthenon/mesh/nx3=8 \
parthenon/meshblock/nx1=8 parthenon/meshblock/nx2=8 parthenon/meshblock/nx3=8 \
radiation/kappa_a_code=0.02 radiation/creduc=30.0 \
parthenon/time/nlim=50 parthenon/time/tlim=100000 \
parthenon/output0/dt=1e9 parthenon/output1/dt=1e-9"
for R in 1.0 2.0; do
  d=x${R%.*}; rm -rf "$d"; mkdir -p "$d"; cd "$d"
  echo "=== ratio=$R ==="
  $BIN -i ../equil.in radiation/planck_ross_ratio=$R $OV >run.log 2>&1
  echo "ratio=$R exit $? cycles=$(grep -cE 'cycle=' run.log)"
  cd ..
done
echo "DONE $(date)"
