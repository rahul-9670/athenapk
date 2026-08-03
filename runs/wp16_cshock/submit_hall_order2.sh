#!/bin/bash
#SBATCH --job-name=wp16hall2
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# WP-16 part 3, ATTEMPT 2 — Hall order of accuracy, with the ladder specified correctly.
#
# ATTEMPT 1 (job 2449287) WAS MIS-SPECIFIED. It refined `nx1` only. But the deck's transverse
# directions are 4 cells across a FIXED extent of 0.03125, so dx2 = dx3 = 0.0078125 always, while
# dx1 = 1/N. The timestep is set by min(dx) over all directions, so for N <= 128 the TRANSVERSE
# spacing was the binding constraint and dt came out IDENTICAL (3.6621e-05) at N = 32, 64 and 128.
# A ladder in which dt does not change is not a temporal-convergence ladder at all.
#
# Note dx2 = dx1 exactly at N = 128 -- i.e. the stock deck is self-consistent at its own native
# resolution, and it was attempt 1 that broke that, not the deck.
#
# CORRECT LADDER: hold the transverse cells CUBIC by scaling the transverse extent with N,
#   x2min = x3min = -2/N,  x2max = x3max = +2/N   =>  dx2 = dx3 = 4/N/4 = 1/N = dx1.
# Then every direction refines together and dt falls as it should.
#
# SEPARATELY, ATTEMPT 1 RAISED A QUESTION THAT MUST BE ANSWERED FIRST. Its N = 128 leg IS the
# stock deck (helicity -1.0 is the deck default, meshblock/nx1 = 128 is the deck value), and it
# did not merely lose accuracy -- it went UNSTABLE:
#     N=32   omega err 1.81e-03   amplitude 1.00e-06 -> 9.07e-07  (-9.3 %, the Ohmic floor)
#     N=64   omega err 3.75e-03   amplitude 1.00e-06 -> 9.08e-07  (-9.25 %)
#     N=128  omega err 6.23e-01   amplitude 1.00e-06 -> 2.65e-01  (+5 DECADES)
# Hall is non-dissipative, so a circular mode must conserve <By^2+Bz^2>; growing it by five
# decades is a numerical instability, not an accuracy loss. Against a recorded historical result
# of 0.4 % for this test, that is either a REGRESSION or a configuration difference. Leg `stock`
# below runs all three stock decks completely unmodified to settle which.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/wp16_cshock
echo "WP-16c attempt 2, job $SLURM_JOB_ID $(date)"; md5sum $B

echo "### PART A — stock decks, ZERO overrides (is attempt 1's N=128 blow-up a regression?)"
for D in hall_whistler hall_whistler_glm hall_whistler_ct; do
  G=$H/stock_$D; rm -rf $G; mkdir -p $G
  env -C $G $B -i /beegfs/u/bbg6470/athenapk/inputs/$D.in > $G/run.log 2>&1
  echo "  $D exit=$?"
  grep -E "measured omega|rel. error|amplitude in/out" $G/run.log | sed 's/^/    /'
done

echo "### PART B — correct ladder: transverse extent scaled so dx2 = dx3 = dx1"
for HEL in -1.0 1.0; do
  TAG=$( [ "$HEL" = "-1.0" ] && echo whi || echo ion )
  for N in 32 64 128 256; do
    E=$(/beegfs/u/bbg6470/venvs/analysis_env/bin/python -c "print(2.0/$N)")
    G=$H/h2_${TAG}_n$N; rm -rf $G; mkdir -p $G
    env -C $G $B -i /beegfs/u/bbg6470/athenapk/inputs/hall_whistler_glm.in \
      problem/diffusion/helicity=$HEL \
      parthenon/mesh/nx1=$N parthenon/meshblock/nx1=$N \
      parthenon/mesh/x2min=-$E parthenon/mesh/x2max=$E \
      parthenon/mesh/x3min=-$E parthenon/mesh/x3max=$E \
      > $G/run.log 2>&1
    echo "  hel=$HEL n=$N exit=$? dt0=$(grep -m1 '^cycle=0' $G/run.log | grep -oE 'dt=[0-9.e+-]+') $(grep 'rel. error' $G/run.log | tail -1)"
    grep -oE "amplitude in/out.*" $G/run.log | sed 's/^/      /'
  done
done
echo "done $(date)"
