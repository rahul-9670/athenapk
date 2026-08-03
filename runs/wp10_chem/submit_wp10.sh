#!/bin/bash
#SBATCH --job-name=wp10chem
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#
# WP-10 — chemistry integrator: substep convergence and element conservation.
#
# WHAT IS ACTUALLY BEING TESTED. `src/chemistry/chemistry.cpp` subcycles the gow17_reduced network
# inside each hydro step under two knobs: `cfl_cool` (fraction of the shortest chemical/cooling
# timescale a substep may take; production 0.1) and `nsub_max` (hard cap on substeps per hydro
# step; production 400). Neither has ever been varied. Two questions follow, and they fail
# differently:
#
#   (a) CONVERGENCE. Is `cfl_cool = 0.1` in the converged regime? Halving it twice must move the
#       species negligibly. If it does not, the production chemistry is substep-limited and the
#       electron abundance x_e -- which feeds eta_A through `ambipolar_coeff = ionization_chem`,
#       i.e. straight into the fossil-field result -- is not converged.
#
#   (b) SATURATION. `nsub_max` is a CAP. If a cell wants more substeps than the cap allows, the
#       integrator stops early and the error is silent -- lowering cfl_cool would then produce NO
#       improvement, because the cap binds first. Raising nsub_max alongside cfl_cool
#       distinguishes "converged" from "capped", which the cfl_cool ladder alone cannot.
#
# Element conservation is the independent check: gow17_reduced carries carbon as C, C+ and CO, so
# x_C + x_Cplus + x_CO must equal `x_Ctot` = 1.6e-4 in every cell at all times, to round-off. It
# is conserved by the network's structure, not by any timestep choice, so a drift that GROWS as
# substeps shrink would indicate an integrator defect rather than a resolution issue.
#
# Same 32^3 smoke deck as the B2/B4 gate (chemistry + radiation + full non-ideal + self-gravity
# all live), nlim=12, so each leg is minutes on CPU.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/wp10_chem
DECK=/beegfs/u/bbg6470/athenapk/runs/b2b4_gate/fhc.in
echo "job $SLURM_JOB_ID $(date)"; md5sum $B

run () {  # name  cfl_cool  nsub_max
  G=$H/$1; rm -rf $G; mkdir -p $G
  ( cd $G && $B -i $DECK \
      chemistry/cfl_cool=$2 chemistry/nsub_max=$3 \
      > run.log 2>&1 )
  echo "$1 (cfl_cool=$2 nsub_max=$3) exit=$?"
}

# (a) the cfl_cool ladder at production nsub_max
run c100  0.1    400      # production
run c050  0.05   400
run c025  0.025  400
run c0125 0.0125 400
# (b) the same two finest rungs with the cap lifted 8x -- if these differ from the ones above,
#     nsub_max was binding and the ladder above was measuring the cap, not the tolerance.
run c025_big  0.025  3200
run c0125_big 0.0125 3200
echo "done $(date)"
