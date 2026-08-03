#!/bin/bash
#SBATCH --job-name=wp16ad
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# WP-16, part 2 — ORDER OF ACCURACY OF THE AMBIPOLAR OPERATOR, measured against an ANALYTIC
# solution rather than by self-convergence.
#
# WHY THIS REPLACES THE C-SHOCK MEASUREMENT. runs/validation/cshock_ad.in cannot answer WP-16:
# its two states violate the steady 1D MHD jump conditions by 327 % (x-momentum), 200 % (y-mom)
# and 292 % (energy) -- only the mass and induction fluxes match -- so the IC is NOT a shock. It
# is an arbitrary discontinuity that blows apart into a Riemann fan, and with outflow BCs on both
# x-faces nothing sustains an upstream state, so the fan simply drains: measured at tlim=400 all
# three resolutions are spatially UNIFORM (max-min < 1e-6) at three DIFFERENT constant states
# (v_x = -1.021 / -1.227 / -1.401). The "self-convergence order" there compares three constants.
# The tlim=40 number (p ~ 0.35) was taken when the structure had already decayed to 14-27 % of
# its initial amplitude, at a resolution-dependent rate. The deck also cites a reference solution
# `runs/validation/cshock_ode.py` that DOES NOT EXIST.
#
# This test instead uses the existing quantitative AD eigenmode (src/pgen/diffusion.cpp,
# iprob=50): B_x = B0 uniform, B_y = amp*sin(kx), which evolves as a damped Alfven mode with a
# CLOSED-FORM amplitude. The pgen already prints `relative error` = |A_meas - A_pred|/|A_pred| at
# t_fin, so the ladder below gives a true (not self-) observed order p = log2(e(h)/e(h/2)).
# It isolates exactly the operator WP-16 is about, with no shock and no boundary influence
# (periodic), and it is cheap: 1D, 4 rungs, seconds each.
#
# Expected: p ~ 2. If the AD operator is 2nd order here, then any low order at a C-shock is a
# property of the DISCONTINUITY (1st order is correct for shock capturing), not of the AD term.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
DECK=/beegfs/u/bbg6470/athenapk/inputs/diffusion_ambipolar.in
H=/beegfs/u/bbg6470/athenapk/runs/wp16_cshock
echo "WP-16b AD eigenmode order ladder, job $SLURM_JOB_ID $(date)"; md5sum $B

for N in 32 64 128 256 512; do
  G=$H/ad_n$N; rm -rf $G; mkdir -p $G
  env -C $G $B -i $DECK \
    parthenon/mesh/nx1=$N parthenon/meshblock/nx1=$N \
    > $G/run.log 2>&1
  echo "n=$N exit=$? $(grep -A0 'relative error' $G/run.log | tail -1)"
done
echo "done $(date)"
