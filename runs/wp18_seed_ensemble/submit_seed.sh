#!/bin/bash
#SBATCH --job-name=wp18
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:h100:8
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
set -o pipefail
#
# WP-18 -- one turbulence-seed realization. Single slot, not self-chaining: the pre-collapse
# leg to t=1.0 takes ~1100 s at 256^3 on 8 GPUs (scaled from the measured 2040 s on 4, job
# 2433514) plus ~3 min one-off IC generation at turb_nmodes=2048.
#
# Env: RUNDIR = this seed's directory (staged by stage.sh)
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

: "${RUNDIR:?set RUNDIR}"
BIN="${BIN:-/beegfs/u/bbg6470/athenapk/build_gpu_wp20/bin/athenaPK}"
DECK=$RUNDIR/fhc_seed.in
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
cd $RUNDIR

echo "WP-18 seed run $(basename $RUNDIR) job $SLURM_JOB_ID $(date)"
echo "binary:"; md5sum $BIN
echo "seed:"; grep "^turb_seed" $DECK
grep -q "Driver completed" run.log 2>/dev/null && { echo "already done"; exit 0; }

stdbuf -oL -eL mpirun -n 8 $MCA $RUNDIR/wrap_mod.sh $BIN -i $DECK -t 01:45:00 \
  parthenon/mesh/nx1=256 parthenon/mesh/nx2=256 parthenon/mesh/nx3=256 \
  parthenon/mesh/do_coalesced_comms=true \
  diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
  diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
  diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
  diffusion/cap_diag=true hydro/mag_diag=true \
  > run.log 2>&1
echo "RUN_EXIT $? $(date)"

# Same IC guard as the root ladder: a binary without turb_ksample silently runs the OLD IC,
# which would make the whole ensemble measure the scatter of a superseded initial condition.
grep -q "k sampling *: k2" run.log \
  && echo "IC_CHECK ok: $(grep -m1 'E(k) ~' run.log | tr -s ' ')" \
  || echo "IC_CHECK **FAILED** -- k2 sampler not reported; this seed's result is VOID."
