#!/bin/bash
#SBATCH --job-name=wp16hd
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# WP-16 part 3, DIAGNOSIS — why does the Hall whistler go unstable in 3D but not in 1D?
#
# WHAT IS ESTABLISHED (jobs 2449287, 2449288):
#
#   Stock decks, ZERO overrides, all at N = 128:
#     hall_whistler.in      (1D, nx2=nx3=1)  omega err 4.2e-03   amplitude 1e-6 -> 9.08e-07  OK
#     hall_whistler_glm.in  (3D, nx2=nx3=4)  omega err 6.2e-01   amplitude 1e-6 -> 2.65e-01  BLOWS UP
#     hall_whistler_ct.in   (3D, nx2=nx3=4)  omega err 6.2e-01   amplitude 1e-6 -> 2.65e-01  BLOWS UP
#
#   The two 3D variants fail IDENTICALLY to every digit, so it is not the divergence-control
#   choice (one is GLM, one is CT). The 1D deck's 4.2e-03 matches the historically recorded
#   "Hall whistler branch 0.4%" -- i.e. **Hall was only ever validated in 1D**, and this is a GAP
#   rather than a regression.
#
#   With the ladder corrected so all three directions refine together, dt scales exactly as dx^2
#   (5.86e-4 / 1.46e-4 / 3.66e-5 / 9.16e-6 for N = 32/64/128/256) -- the Hall timestep limiter is
#   working. The instability is therefore NOT a dt bug. It onsets between N = 32 and N = 64 and
#   worsens monotonically:
#     whistler  N=32 amp -9.3 %   N=64 amp +76 %   N=128 amp +2.65e5   N=256 crash
#     ion-cyc   N=32 amp -7.2 %   N=64 amp +17x    N=128 amp +2.60e5   N=256 crash
#
# THE HYPOTHESIS, and how this job tries to KILL it. Hall is dispersive: the whistler frequency
# goes as omega ~ eta_H k^2, so the fastest mode the grid supports is at k_max ~ pi/dx and grows
# like 1/dx^2. The Ohmic floor that stabilises it damps like eta_O k^2 -- the same power -- so a
# FIXED floor should in principle work at every resolution, and the observed 1/dx-ish onset would
# then NOT be explained by "the floor is too small". If raising the floor DOES restore stability
# with a threshold that scales predictably, the mechanism is a floor calibration and production
# needs a resolution-dependent floor. If raising it does NOT, the fault is in the 3D Hall EMF
# itself and no amount of floor will fix it. Those two outcomes are distinguishable and this job
# distinguishes them.
#
# PART A: 1D at high resolution -- is the failure 3D-SPECIFIC or just RESOLUTION-specific?
# PART B: 3D with the Ohmic floor raised -- does more floor buy stability, and how much is needed?
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/wp16_cshock
D1=/beegfs/u/bbg6470/athenapk/inputs/hall_whistler.in
D3=/beegfs/u/bbg6470/athenapk/inputs/hall_whistler_glm.in
echo "WP-16c diagnosis, job $SLURM_JOB_ID $(date)"; md5sum $B

report () { grep -E "rel. error \(omega\)|amplitude in/out" $1/run.log | sed 's/^/      /'; }

echo "### PART A — 1D ladder (is the blow-up 3D-specific, or does 1D fail too at high N?)"
for N in 128 256 512 1024; do
  G=$H/d1_n$N; rm -rf $G; mkdir -p $G
  env -C $G $B -i $D1 parthenon/mesh/nx1=$N parthenon/meshblock/nx1=$N > $G/run.log 2>&1
  echo "  1D n=$N exit=$?"; report $G
done

echo "### PART B — 3D, Ohmic floor raised (does more floor restore stability?)"
for N in 64 128; do
  E=$(/beegfs/u/bbg6470/venvs/analysis_env/bin/python -c "print(2.0/$N)")
  for FL in 0.05 0.1 0.2 0.4 0.8; do
    G=$H/d3_n${N}_f$FL; rm -rf $G; mkdir -p $G
    env -C $G $B -i $D3 \
      parthenon/mesh/nx1=$N parthenon/meshblock/nx1=$N \
      parthenon/mesh/x2min=-$E parthenon/mesh/x2max=$E \
      parthenon/mesh/x3min=-$E parthenon/mesh/x3max=$E \
      diffusion/hall_ohmic_floor_code=$FL \
      > $G/run.log 2>&1
    echo "  3D n=$N floor=$FL exit=$?"; report $G
  done
done
echo "done $(date)"
