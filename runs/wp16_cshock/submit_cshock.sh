#!/bin/bash
#SBATCH --job-name=wp16cs
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#
# WP-16 — close the ambipolar C-shock convergence gate.
#
# WHAT WAS ALREADY ESTABLISHED (DEV_LOG WS-5b): the harness works, and with gamma=1.4 +
# weak shock (vxl=2, ambipolar_coeff_code=0.15) + thin-2D (nx2=4; the pgen SEGFAULTS at nx2=1,
# a real bug) the AD C-shock runs clean at 128/256/512 and forms a continuous C-type structure
# whose thickness scales correctly with the ambipolar coefficient.
#
# WHY IT WAS NOT CLOSED: at tlim=40 the measured self-convergence order was only ~0.35 (v_x) /
# 0.65 (B_y). The diagnosis was that a thickness-~12 C-shock has a long relaxation time and the
# three resolutions relax at DIFFERENT rates, so the L1 difference between them is a transient,
# not truncation error. That is a hypothesis about the measurement, and it makes a falsifiable
# prediction: run far longer, confirm d/dt -> 0, and the order should rise toward 2. If it does
# NOT rise once the structure is genuinely steady, the low order is real and is a finding.
#
# THIS RUN: tlim = 400 (10x), with frequent history output so d/dt can be checked directly
# rather than assumed. Three resolutions for the self-convergence triple.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
DECK=/beegfs/u/bbg6470/athenapk/runs/validation/cshock_ad.in
H=/beegfs/u/bbg6470/athenapk/runs/wp16_cshock
TLIM=${TLIM:-400.0}
echo "WP-16 C-shock relaxation, tlim=$TLIM, job $SLURM_JOB_ID $(date)"; md5sum $B

for N in 128 256 512; do
  G=$H/n$N; rm -rf $G; mkdir -p $G
  env -C $G $B -i $DECK \
    parthenon/mesh/nx1=$N parthenon/meshblock/nx1=$((N/4)) \
    parthenon/time/tlim=$TLIM \
    > $G/run.log 2>&1
  echo "n=$N exit=$? $(date +%H:%M:%S)"
done
echo "done $(date)"
