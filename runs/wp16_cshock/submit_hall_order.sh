#!/bin/bash
#SBATCH --job-name=wp16hall
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# WP-16, part 3 — ORDER OF ACCURACY OF THE HALL OPERATOR.
#
# SCOPE CORRECTION. VALIDATION_PLAN.md's WP-16 reads "Complete the Hall C-shock test [...] Hall is
# the least-validated of the three non-ideal terms". Parts 1 and 2 of this WP addressed the
# AMBIPOLAR operator -- that was necessary (the C-shock deck `runs/validation/cshock_ad.in` is
# ambipolar and turned out to be invalid), but it did not answer the Hall half. This does.
#
# METHOD, identical to part 2 so the two are directly comparable: measure against an ANALYTIC
# solution rather than by self-convergence. `src/pgen/diffusion.cpp` iprob=60 initialises a
# circularly polarised Hall eigenmode (B_x = B0, B_y = amp cos kx, B_z = h amp sin kx) whose
# dispersion relation is closed-form, and reports
#     rel. error (omega) = |omega_meas - omega_analytic| / omega_analytic
# with omega_meas recovered from the phase, plus the in/out amplitude (the Hall term is
# NON-DISSIPATIVE, so a circular mode must conserve <By^2+Bz^2> -- a second, independent check
# that does not depend on the dispersion relation being right).
#
# Both helicity branches are run because they are physically different waves and the Hall term
# splits them: helicity -1 is the whistler (fast) branch, +1 the ion-cyclotron (slow) branch.
# A scheme could get one right and the other wrong.
#
# CONFIGURATION NOTES (these are constraints, not choices):
#   * Hall is DISPERSIVE => `diffusion/integrator = unsplit` only; it errors under rkl2.
#   * It needs an Ohmic floor (`hall_ohmic_floor_code = 0.05`) for stability on the cell-centred
#     grid. That floor is a real physical resistivity in the dispersion relation the pgen solves
#     (etaO enters `bc` at diffusion.cpp:385), so it does NOT bias the comparison -- but it does
#     mean this measures "Hall + a small Ohmic term", which is exactly the production combination.
#   * The mode needs 3D (2D CT freezes B_z), hence nx2 = nx3 = 4 in the deck.
#
# Expected: p ~ 2. Only nx1 is refined -- the transverse directions carry no structure.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
DECK=/beegfs/u/bbg6470/athenapk/inputs/hall_whistler_glm.in
H=/beegfs/u/bbg6470/athenapk/runs/wp16_cshock
echo "WP-16c Hall order ladder, job $SLURM_JOB_ID $(date)"; md5sum $B

for HEL in -1.0 1.0; do
  TAG=$( [ "$HEL" = "-1.0" ] && echo whi || echo ion )
  for N in 32 64 128 256 512; do
    G=$H/hall_${TAG}_n$N; rm -rf $G; mkdir -p $G
    env -C $G $B -i $DECK \
      problem/diffusion/helicity=$HEL \
      parthenon/mesh/nx1=$N parthenon/meshblock/nx1=$N \
      > $G/run.log 2>&1
    echo "hel=$HEL n=$N exit=$? $(grep 'rel. error' $G/run.log | tail -1)"
  done
done
echo "done $(date)"
