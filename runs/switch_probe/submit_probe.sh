#!/bin/bash
#SBATCH --job-name=swprobe
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# Measure what the two production switches actually DO, before committing the GPU re-baseline
# to them. Four legs, one knob at a time, identical deck otherwise:
#
#   base   outflow + absolute tolerance   (current production)
#   diode  diode   + absolute
#   rel    outflow + relative
#   both   diode   + relative             (the proposed new production)
#
# WHY THIS IS NOT A FORMALITY. With relative_residual = true the criterion becomes
#   rms_res < 1e-6 * rms(rhs),  rhs = 4 pi G rho.
# Through collapse rms(rhs) grows by DECADES, so the relative tolerance gets LOOSER, not
# tighter -- the opposite of what B3's framing suggests. The absolute error in phi therefore
# grows with density, and phi's gradient is what drives the collapse. `grav-iters` and
# `grav-res` (solver_diag) are the direct evidence of which way it goes.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/switch_probe
DECK=/beegfs/u/bbg6470/athenapk/runs/b2b4_gate/fhc.in
echo "job $SLURM_JOB_ID $(date)"; md5sum $B

run () {  # name  bc  relative
  G=$H/$1; rm -rf $G; mkdir -p $G
  BCS=""; for f in ix1 ox1 ix2 ox2 ix3 ox3; do BCS="$BCS parthenon/mesh/${f}_bc=$2"; done
  env -C $G $B -i $DECK \
    hydro/cons_diag=true self_gravity/solver_diag=true \
    self_gravity/solver_params/relative_residual=$3 \
    $BCS > $G/run.log 2>&1
  echo "$1 (bc=$2 relative=$3) exit=$?"
}

run base  outflow false
run diode  diode  false
run rel   outflow true
run both   diode  true
echo "done $(date)"
