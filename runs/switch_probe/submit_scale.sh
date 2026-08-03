#!/bin/bash
#SBATCH --job-name=scaleprobe
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#
# B3 DECIDED DIRECTLY: absolute vs relative gravity-solver tolerance, as a function of SOURCE
# SCALE. This replaces three failed attempts at reaching collapse in a cheap run.
#
# WHY THE COLLAPSE ROUTE KEPT FAILING. The question is "what happens to the Poisson solve as
# rms(rhs) = rms(4 pi G rho) grows by decades", and every attempt to get there by actually
# collapsing gas hit a different wall:
#   attempt 1 (full physics, numlevel=5):  290 s/step, cycle 12 in 47 min -- could not reach it.
#   attempt 2 (adiabatic gamma=5/3):       pressure-supported core; nblk froze at 288 by t=0.98.
#   attempt 3 (tabulated EOS, no rad):     refused at init -- eos=hydrogen REQUIRES radiation.
#   attempt 4 (adiabatic gamma=1.1):       nblk ALSO froze at 288; grav-iters plateaued at 8/200.
#
# The mistake was insisting on a physically-driven source. The Poisson operator does not care how
# the gas got dense -- it only sees the source. And because the operator is LINEAR, scaling
# rho -> lambda*rho scales phi and the residual by exactly lambda. That makes the whole question
# answerable with a STATIC source of controlled amplitude:
#
#   * RELATIVE criterion (rms_res < tol * rms(rhs)):  invariant under lambda. Predicted:
#     iteration count and final residual/rms(rhs) IDENTICAL at every rho_in.
#   * ABSOLUTE criterion (rms_res < tol):  the threshold is effectively divided by lambda.
#     Predicted: iteration count RISES with rho_in and eventually pins at max_iterations,
#     after which the solve exits at an uncontrolled accuracy.
#
# runs/validation_ws5a/sphere_multipole.in is exactly the right vehicle: a static uniform sphere,
# multipole BCs (production's), `nlim = 1` so only the first Poisson solve runs, and `rho_in` as a
# free amplitude. Each leg is seconds.
#
# rho_in spans 1 -> 1e8. Production reaches rho/rho_ambient ~ 1e10 in code units by the first
# core, so this brackets the real range.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/switch_probe/scale
DECK=/beegfs/u/bbg6470/athenapk/runs/validation_ws5a/sphere_multipole.in
rm -rf $H; mkdir -p $H
echo "job $SLURM_JOB_ID $(date)"; md5sum $B

# ATTEMPT 1 (2448614) WAS INVALID and its result must not be used: it passed
# `relative_residual_tolerance=$3` on the command line in EVERY leg. `BiCGSTABParams`
# (bicgstab_solver.hpp:56-66) reads that key with GetOrAddReal, so supplying it kept the
# RELATIVE test live even in the legs labelled "absolute" -- and since the accept test is an OR
# (`rms_res < rel_tol || rms_res < abs_tol`, :411) the relative branch simply won. The tell was
# that abs6 and rel6 came out IDENTICAL to every printed digit, and that abs6 "converged" at
# rho_in=1e8 with grav-res = 4.43e-6, four times ABOVE its own 1e-6 absolute ceiling.
#
# The deck (`sphere_multipole.in`) is a second source of the same contamination: it sets
# `absolute_residual_tolerance = 1e-12`, `relative_residual_tolerance = 1e-8` and
# `relative_residual = true`. So BOTH keys must be pinned explicitly in every leg to reproduce
# what the code does when a deck omits them -- the unselected one has to be driven to 0.
#
# name  relative  tol  rel_tol  abs_tol
for CRIT in "abs6 false 1.0e-6 0.0 1.0e-6" "abs8 false 1.0e-8 0.0 1.0e-8" "rel6 true 1.0e-6 1.0e-6 0.0"; do
  set -- $CRIT
  for R in 1 1e2 1e4 1e6 1e8 1e10; do
    G=$H/$1_r$R; mkdir -p $G
    # The deck has NO hst block (its output0 is hdf5), so solver_diag's grav-iters/grav-res
    # columns had nowhere to go -- add one. And scale the pressure with rho_in: at fixed
    # pressure = 1 the rho_in >= 1e6 legs hit "Got negative pressure" in the hydro update
    # AFTER the Poisson solve, killing the run before any output. Scaling pressure keeps the
    # sphere at constant temperature, which is both survivable and the physically sane choice;
    # it does not touch the Poisson source, which is rho alone.
    ( cd $G && $B -i $DECK \
        problem/poisson_test/rho_in=$R \
        problem/poisson_test/pressure=$R \
        parthenon/output1/file_type=hst \
        parthenon/output1/dt=1.0e-30 \
        hydro/pfloor=1.0e-10 \
        self_gravity/solver_diag=true \
        self_gravity/solver_params/relative_residual=$2 \
        self_gravity/solver_params/residual_tolerance=$3 \
        self_gravity/solver_params/relative_residual_tolerance=$4 \
        self_gravity/solver_params/absolute_residual_tolerance=$5 \
        self_gravity/solver_params/max_iterations=200 \
        > run.log 2>&1 )
    echo "$1 rho_in=$R exit=$?"
  done
done
echo "done $(date)"
