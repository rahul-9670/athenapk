#!/bin/bash
#SBATCH --job-name=b10ver
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#
# B10 VERIFICATION. The claim: `nsub_max` silently overrides `chemistry/cfl_cool`, and the
# existing counter could not see it because it reported only loop exhaustion, which the
# dt_floor makes unreachable. The instrument now reports floor engagement.
#
# Falsifiable predictions:
#   (1) production settings (cfl_cool=0.1, nsub_max=400) now EMIT the warning, where the old
#       binary printed nothing;
#   (2) the same deck with nsub_max raised far enough emits it LESS or not at all;
#   (3) abundances are UNCHANGED by the instrument itself -- it must be diagnostic-only.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=8 OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash
B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/wp10_chem
DECK=/beegfs/u/bbg6470/athenapk/runs/b2b4_gate/fhc.in
echo "job $SLURM_JOB_ID $(date)"; md5sum $B
run () { G=$H/$1; rm -rf $G; mkdir -p $G
  ( cd $G && $B -i $DECK chemistry/cfl_cool=$2 chemistry/nsub_max=$3 > run.log 2>&1 )
  echo "$1 (cfl_cool=$2 nsub_max=$3) exit=$? warn=$(grep -c 'WARNING Chemistry' $G/run.log)"
  grep -m1 -A1 'WARNING Chemistry' $G/run.log | head -2
}
run v_prod   0.1    400        # production settings
run v_big    0.1    100000     # cap effectively lifted
echo "done $(date)"
