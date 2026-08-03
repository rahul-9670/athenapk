#!/bin/bash
#SBATCH --job-name=wp16ht
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#
# WP-16 part 3, THRESHOLD — how much Ohmic floor does 3D Hall need, as a fraction of eta_H?
#
# ESTABLISHED SO FAR: 3D Hall is unstable at the production floor 0.05 with the stress value
# eta_H = Q_H*B/rho = 0.5, and raising the floor to 0.1 fixes it -- at BOTH N=64 and N=128, with
# results then agreeing to ~5 % between the two resolutions. So the requirement is a fixed RATIO
# eta_floor/eta_H, not a resolution-dependent number, and it sits between 0.1 and 0.2.
#
# This scan measures the threshold from the other side: hold the floor at production's 0.05 and
# vary Q_H. Prediction from the floor scan: the onset is at eta_floor/eta_H ~ 0.1-0.2, i.e.
# Q_H ~ 0.25-0.5. If instead it onsets somewhere unrelated to that ratio, the ratio picture is
# wrong and the mechanism is something else.
#
# WHY IT MATTERS: production sets a FIXED hall_ohmic_floor_code = 0.05 while eta_H is computed
# from the ionization model and varies over the domain. A fixed floor is only safe while eta_H
# stays below the threshold. Measured on prod_v9: max |eta_H| = 4.316e-02 over 78.7M cells, and
# 0.0000 % of cells exceed 0.05 -- so production sits below the stress value by 11.6x. This scan
# turns that into a quantified margin.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash
B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/wp16_cshock
D3=/beegfs/u/bbg6470/athenapk/inputs/hall_whistler_glm.in
echo "job $SLURM_JOB_ID $(date)"; md5sum $B
N=128; E=$(/beegfs/u/bbg6470/venvs/analysis_env/bin/python -c "print(2.0/$N)")
for QH in 0.5 0.4 0.3 0.25 0.2 0.15 0.1 0.05; do
  G=$H/th_q$QH; rm -rf $G; mkdir -p $G
  env -C $G $B -i $D3 \
    parthenon/mesh/nx1=$N parthenon/meshblock/nx1=$N \
    parthenon/mesh/x2min=-$E parthenon/mesh/x2max=$E \
    parthenon/mesh/x3min=-$E parthenon/mesh/x3max=$E \
    diffusion/hall_coeff_code=$QH \
    diffusion/hall_ohmic_floor_code=0.05 \
    > $G/run.log 2>&1
  R=$(/beegfs/u/bbg6470/venvs/analysis_env/bin/python -c "print(f'{0.05/$QH:.3f}')")
  echo "  Q_H=$QH  floor/eta_H=$R  exit=$?  $(grep -oE 'rel. error \(omega\).*' $G/run.log|tail -1)"
  grep -oE "amplitude in/out.*" $G/run.log | sed 's/^/      /'
done
echo "done $(date)"
