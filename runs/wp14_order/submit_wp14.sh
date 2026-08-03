#!/bin/bash
#SBATCH --job-name=wp14
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#
# WP-14 order-of-accuracy ladder. Submitted to SLURM rather than run on the front-end: the
# front-end shell sits in a 1-CPU cgroup (`nproc` = 1) with a load average already >2, so the
# 256x128x128 rung is both slow there and antisocial. build_cpu is a Kokkos OpenMP build, so
# threads help; L1 error is a volume-weighted sum whose reduction order can shift the last
# ulp, which is irrelevant at the 1e-8 error level we are fitting an ORDER to.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/wp14_order
echo "WP-14 ladder job $SLURM_JOB_ID  waves='${WAVES}'  threads=$OMP_NUM_THREADS  $(date)"
md5sum $B

for WF in ${WAVES:?set WAVES}; do
  for N in 16 32 64 128; do
    G=$H/wf${WF}_n${N}
    # Do not clobber a rung that already completed (the pilot ran wave 0 on the front-end).
    if [ -s $G/linearwave-errors.dat ]; then echo "wf=$WF n=$N already done -> skip"; continue; fi
    rm -rf $G; mkdir -p $G
    M=16; [ $N -le 16 ] && M=8
    env -C $G $B -i $H/lw_mhd.in \
      problem/linear_wave/wave_flag=$WF \
      parthenon/mesh/nx1=$((N*2)) parthenon/mesh/nx2=$N parthenon/mesh/nx3=$N \
      parthenon/meshblock/nx1=$M parthenon/meshblock/nx2=$M parthenon/meshblock/nx3=$M \
      > $G/run.log 2>&1
    echo "wf=$WF n=$N exit=$? $(date +%H:%M:%S)"
  done
done
echo "WP-14 ladder done $(date)"
