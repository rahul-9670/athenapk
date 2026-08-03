#!/bin/bash
#SBATCH --job-name=wp14ent
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# WP-14, entropy family (wave_flag = 3) — RERUN. The first attempt (job 2448312) was a NULL TEST.
#
# ROOT CAUSE. src/pgen/linear_wave_mhd.cpp:164-168 reinterprets `tlim` as a number of wave
# periods:  tlim <- lambda / |ev[wave_flag]| * tlim.  The entropy mode's eigenvalue IS the
# background flow speed, and lw_mhd.in sets `vflow = 0.0`, so ev[3] = 0 and
#   tlim = lambda / 0 = inf.
# The run.logs confirm it verbatim: `tlim=inf nlim=100000`. The n=16 and n=32 rungs therefore
# terminated on nlim after t = 1406 and t = 703 (not one period), and n=64 was still running at
# t = 71.7 / 20396 cycles when the 4 h job wall-clocked. The 1e-12 "errors" those rungs reported
# are how far a STATIONARY profile drifts over 100k cycles -- a real (and reassuring) number, but
# not a truncation error, which is why the observed order came out -0.16.
#
# FIX: give the entropy mode something to advect on. vflow = 1.0 makes ev[3] = 1, lambda = 1
# (the box is built so the oblique wavelength is exactly 1.0), hence tlim = 1.0 = one period and
# the exact solution at t=tlim is again the initial condition.
#
# NOTE ON COMPARABILITY: vflow != 0 shifts the background state relative to the fast/Alfven/slow
# ladders, which ran at vflow = 0. That is unavoidable -- the entropy mode is identically static
# without a background flow -- and it is the standard way this family is tested. The measured
# order is still the order of the same scheme.
#
# Families 4/5/6 (slow+, Alfven+, fast+) are the mirror-image branches of 2/1/0 and are
# deliberately NOT run: they exercise the same reconstruction/flux path with the sign of the
# eigenvalue flipped and add no independent information about the order.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/wp14_order
echo "WP-14 entropy ladder job $SLURM_JOB_ID $(date)"; md5sum $B

for N in 16 32 64 128; do
  G=$H/wf3_n${N}
  rm -rf $G; mkdir -p $G          # the previous wf3_* dirs hold the null-test runs; replace them
  M=16; [ $N -le 16 ] && M=8
  env -C $G $B -i $H/lw_mhd.in \
    problem/linear_wave/wave_flag=3 \
    problem/linear_wave/vflow=1.0 \
    parthenon/mesh/nx1=$((N*2)) parthenon/mesh/nx2=$N parthenon/mesh/nx3=$N \
    parthenon/meshblock/nx1=$M parthenon/meshblock/nx2=$M parthenon/meshblock/nx3=$M \
    > $G/run.log 2>&1
  echo "wf=3 n=$N exit=$? tlim=$(grep -m1 -o 'tlim=[^ ]*' $G/run.log) $(date +%H:%M:%S)"
done
echo "done $(date)"
