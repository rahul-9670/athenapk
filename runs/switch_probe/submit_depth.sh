#!/bin/bash
#SBATCH --job-name=depthprobe
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#
# THE DECIDING EXPERIMENT for B3: absolute vs relative gravity-solver tolerance AT DEPTH.
#
# Both shallow probes (jobs 2448322, 2448494) agree with each other and settle nothing about
# production, because they only reach t = 1.11 with rms(rho) = 0.435:
#   * `rel` (relative 1e-6) and `tol8` (absolute 1e-8) are BYTE-IDENTICAL on every physics column.
#   * both differ from `base` (absolute 1e-6, current production) by <= 1.8e-5 relative.
#   * cost: 4.92 -> 7.50 mean BiCGSTAB iterations (+52 %) for an 83x tighter residual.
# So at shallow depth the choice does not matter. The whole question is what happens as the
# source grows.
#
# THE TWO FAILURE MODES, stated so the run can falsify one of them.
#   rms(res) at convergence scales with rms(rhs) = rms(4 pi G rho) = rms(rho) [four_pi_G = 1].
#   rho_crit = 1e-13 / 5.467e-19 = 1.83e5 in code units, and the criterion is evaluated over
#   GetTotalCells() -- under Jeans AMR most cells ARE core cells, so rms(rhs) climbs by decades.
#
#   ABSOLUTE 1e-6 (production today): a fixed absolute threshold against a growing source is an
#     ever-TIGHTER relative demand -- 2.3e-6 relative at the smoke scale, ~1e-11 at rho_crit.
#     Predicted failure: the solver stops hitting tolerance, burns all 200 iterations, and exits
#     at an uncontrolled accuracy. Signature: `grav-iters` -> 200 and `grav-nonconv` > 0.
#
#   RELATIVE 1e-6: a fixed FRACTIONAL accuracy in phi at every depth, which is the physically
#     meaningful demand since phi and its gradient scale with the source. The absolute error in
#     phi does grow -- but so does phi. Predicted behaviour: bounded iteration count, no
#     nonconvergence. Signature: `grav-iters` stays O(10), `grav-nonconv` stays 0.
#
# Whichever signature appears decides B3. This supersedes the earlier reasoning in
# submit_tol.sh, which argued from the absolute error alone and did not weigh it against the
# growing physical scale.
#
# CONFIG: same 32^3 base deck, but numlevel 2 -> 5 and nlim removed, so the core actually
# refines and rho climbs. `-t` bounds each leg; a leg that runs out of time still produces the
# history needed, since the question is answered on the way down, not at the end.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/switch_probe
DECK=/beegfs/u/bbg6470/athenapk/runs/b2b4_gate/fhc.in
echo "job $SLURM_JOB_ID $(date)"; md5sum $B

run () {  # name  relative  tol
  G=$H/$1; rm -rf $G; mkdir -p $G
  ( cd $G && mpirun -n 4 --oversubscribe $B -i $DECK -t 02:20:00 \
      hydro/cons_diag=true self_gravity/solver_diag=true \
      self_gravity/solver_params/relative_residual=$2 \
      self_gravity/solver_params/residual_tolerance=$3 \
      parthenon/mesh/numlevel=5 \
      parthenon/time/nlim=-1 parthenon/time/tlim=1.5 \
      > run.log 2>&1 )
  echo "$1 (relative=$2 tol=$3) exit=$?  last: $(grep '^cycle=' $G/run.log | tail -1)"
}

run deep_abs6  false 1.0e-6     # current production
run deep_abs8  false 1.0e-8     # the tightened-absolute candidate
run deep_rel6  true  1.0e-6     # the relative candidate
echo "done $(date)"
