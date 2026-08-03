#!/bin/bash
#SBATCH --job-name=tolprobe
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# FOLLOW-UP to submit_probe.sh, which REJECTED `relative_residual = true` for production.
#
# WHY IT WAS REJECTED (read bicgstab_solver.hpp:56-66 and :405-411). The switch is not additive,
# it REPLACES the criterion: with relative_residual=false the code sets abs_tol=residual_tolerance
# and rel_tol=0; with relative_residual=true it sets rel_tol=residual_tolerance and
# abs_tol=**0**. The accept test is `rms_res < rel_tol || rms_res < abs_tol`, so exactly one of
# the two is ever live, and there is NO way to demand that both hold.
#   rel_tol_effective = 1e-6 * sqrt(rhs2/GetTotalCells()) = 1e-6 * rms(rho)   [four_pi_G = 1]
# Measured on the probe deck: rms(rho) = 0.435 at t=1.11, so relative is 2.3x TIGHTER *today*.
# But rho_crit = 1e-13/5.467e-19 = 1.83e5 in code units, and that rms is weighted by CELL COUNT,
# which under Jeans AMR is overwhelmingly core cells. At depth rel_tol therefore becomes
# 3-5 decades LOOSER than the absolute 1e-6 it replaced -- the criterion would silently relax
# exactly where the solve matters most. So B3's fix is NOT to flip this switch.
#
# THIS PROBE tests the alternative: keep the absolute criterion and TIGHTEN it. B3's evidence
# stands -- `grav-res` runs at 70-84% of the 1e-6 tolerance by cycle 12 on a 32^3 smoke, so the
# current setting is already marginal and will be breached as rhs grows.
#
#   tol6  = current production        (1e-6, absolute)   -- reference, = probe leg `base`
#   tol8  = tightened                 (1e-8, absolute)
#   tol10 = tightened further         (1e-10, absolute)
#   d8    = tol8 + diode BC           (the production candidate)
#
# What we need out of it: (i) the iteration cost of each tolerance, (ii) confirmation that the
# state shift is small (the `rel` leg, which reached ~1e-8, moved the state by only ~1e-5
# relative -- so tightening should be a near-no-op today while protecting the deep phase), and
# (iii) whether 1e-10 is affordable or hits max_iterations.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/switch_probe
DECK=/beegfs/u/bbg6470/athenapk/runs/b2b4_gate/fhc.in
echo "job $SLURM_JOB_ID $(date)"; md5sum $B

run () {  # name  tol  bc
  G=$H/$1; rm -rf $G; mkdir -p $G
  BCS=""; for f in ix1 ox1 ix2 ox2 ix3 ox3; do BCS="$BCS parthenon/mesh/${f}_bc=$3"; done
  env -C $G $B -i $DECK \
    hydro/cons_diag=true self_gravity/solver_diag=true \
    self_gravity/solver_params/relative_residual=false \
    self_gravity/solver_params/residual_tolerance=$2 \
    $BCS > $G/run.log 2>&1
  echo "$1 (tol=$2 bc=$3) exit=$?"
}

run tol8  1.0e-8  outflow
run tol10 1.0e-10 outflow
run d8    1.0e-8  diode
echo "done $(date)"
