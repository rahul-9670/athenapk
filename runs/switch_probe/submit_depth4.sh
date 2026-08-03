#!/bin/bash
#SBATCH --job-name=depth4
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#
# THE DECIDING EXPERIMENT for B3 (absolute vs relative gravity tolerance), attempt 2.
#
# ATTEMPT 1 (job 2448512) WAS CANCELLED, not failed: with the full production physics stack
# (M1 radiation + gow17 chemistry + Ohm/Hall/AD + tabulated EOS) on top of numlevel=5, the leg
# ran at 290 s/step and reached only cycle 12 / t = 0.55 in 47 minutes. It could not have reached
# collapse inside its 2h20m budget, so it would have produced no evidence about depth at all.
#
# ATTEMPT 2 WAS ALSO CANCELLED, for the OPPOSITE reason: trimming to `hydro/eos=adiabatic`
# gamma=5/3 removed the thermodynamics that lets the collapse RUN AWAY. gamma=5/3 heating builds a
# pressure-supported core almost immediately -- nblk froze at 288 from cycle 45 (t=0.98) onward,
# the Jeans criterion stopped triggering refinement, and dt decayed only 2.2x over the next 200
# cycles. rms(rhs) was never going to grow by decades.
#
# ATTEMPT 3 (keep the TABULATED eos=hydrogen, drop radiation) DIED AT INIT in 4 s:
#   "eos=hydrogen requires <physics> radiation=true (else the barotropic cooling overwrites the
#    table EOS thermal energy)"
# -- a hard code constraint, so the table and radiation cannot be separated.
#
# ATTEMPT 4 (this one) stops trying to be a small version of production and instead poses the
# minimal problem that PRODUCES THE CONDITION UNDER TEST: a soft adiabatic collapse, gamma = 1.1.
# For gamma < 4/3 pressure cannot halt self-gravitating collapse, so the central density runs away
# (Larson-Penston) until the grid caps it -- numlevel=6 buys 64x in linear resolution, i.e. many
# decades in rho and therefore in rms(rhs), which is the ONLY property of the state the Poisson
# solver's difficulty depends on. Physical fidelity is irrelevant here: the question is whether
# BiCGSTAB under an ABSOLUTE tolerance stops converging when its source grows by decades, and the
# operator does not care how the gas got dense.
#
# LESSON from attempts 2 and 3: when trimming a problem to isolate a solver, keep whatever drives
# the state variable the solver's difficulty depends on -- but it does not have to be the
# PRODUCTION mechanism that drives it, just a cheap one that drives it the same way.
#
# THE FIX IS TO TRIM THE PROBLEM TO THE QUESTION. The question is entirely about the Poisson
# solver: how does the residual at convergence behave as rms(rhs) = rms(4 pi G rho) grows by
# decades? Radiation, chemistry and non-ideal MHD do not enter that; they only set how fast the
# gas gets there. So they are switched OFF and an adiabatic EOS replaces the table. What is kept
# is exactly what matters: self-gravity with multipole BCs, AMR (numlevel=6, deeper than
# attempt 1), and the same BE-sphere IC.
#
# THE PREDICTION UNDER TEST, restated. The accept test (bicgstab_solver.hpp:405-411) is
#   rms_res < rel_tol || rms_res < abs_tol,  rel_tol = relative_residual_tolerance * rms(rhs),
# and the two are MUTUALLY EXCLUSIVE by construction (lines 56-66 zero whichever is unselected).
# rms_res at convergence scales with rms(rhs), so:
#   * absolute 1e-6  -> an ever-TIGHTER relative demand as rho grows. Predicted signature:
#     grav-iters climbs toward max_iterations=200 and grav-nonconv leaves 0.
#   * relative 1e-6  -> a FIXED fractional accuracy at every depth. Predicted signature:
#     grav-iters stays O(10), grav-nonconv stays 0.
# For a LINEAR operator this is close to a theorem -- scaling rho by lambda scales phi and the
# residual by lambda, leaving the relative criterion invariant and multiplying the absolute
# demand by lambda -- but the density PROFILE also sharpens during collapse, which changes the
# conditioning, and that part is not analytic. Hence the measurement.
#
# Shallow probes already established that the choice does not matter at t ~ 1.1 with
# rms(rho) = 0.435: `rel` and `tol8` (absolute 1e-8) agree with each other and both differ from
# production's absolute 1e-6 by only 1.8e-5 relative, at a cost of 4.92 -> 7.50 mean iterations.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/switch_probe
DECK=/beegfs/u/bbg6470/athenapk/runs/b2b4_gate/fhc.in
echo "job $SLURM_JOB_ID $(date)"; md5sum $B

# Gravity + hydro + AMR only. Everything that does not touch the Poisson solve is off.
TRIM="physics/radiation=false physics/chemistry=false physics/dust=false
      diffusion/ambipolar=none diffusion/resistivity=none diffusion/hall=none
      diffusion/integrator=unsplit
      hydro/eos=adiabatic hydro/gamma=1.1 hydro/nscalars=0
      parthenon/mesh/numlevel=6 parthenon/time/nlim=-1 parthenon/time/tlim=3.0"

run () {  # name  relative  tol
  G=$H/$1; rm -rf $G; mkdir -p $G
  ( cd $G && mpirun -n 8 $B -i $DECK -t 02:20:00 \
      hydro/cons_diag=true self_gravity/solver_diag=true \
      self_gravity/solver_params/relative_residual=$2 \
      self_gravity/solver_params/residual_tolerance=$3 \
      $TRIM > run.log 2>&1 )
  echo "$1 (relative=$2 tol=$3) exit=$?  last: $(grep '^cycle=' $G/run.log | tail -1)"
}

run d4_abs6 false 1.0e-6     # current production
run d4_abs8 false 1.0e-8     # tightened-absolute candidate
run d4_rel6 true  1.0e-6     # relative candidate
echo "done $(date)"
