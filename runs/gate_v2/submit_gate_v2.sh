#!/bin/bash
#SBATCH --job-name=gatev2
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# OFF-STATE GATE for the new GPU binary.
#
#   old = build_gpu/bin/athenaPK      md5 49d9c257   (current production; gravity fix + 4 diags)
#   new = build_gpu_v2/bin/athenaPK   md5 f181c0a1   (old + B1 diode BC + B2 solver-converged
#                                                     flag & unconditional warning + B3 startup
#                                                     tolerance notice + B4 coarse-EOS-table
#                                                     warning)
#
# Every one of those four is opt-in or diagnostic-only, so with the production deck unchanged the
# two binaries MUST produce bit-identical history. This is the same gate the CPU batch used, and
# it is what licenses swapping the production binary without re-baselining anything.
#
# It is run on the PRODUCTION deck (multipole self-gravity + ambipolar + RKL2 + AMR), not a smoke
# deck, because that is the code path the fixes touch. nlim=40 keeps it to a few minutes while
# still crossing the first AMR regrid.
#
# 2 GPUs, not 8: WP-12 measured decomposition invariance, and 4->8 scales at only 74 % efficiency,
# so a bit-for-bit comparison has no reason to burn 8 cards.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

H=/beegfs/u/bbg6470/athenapk/runs/gate_v2
DECK=/beegfs/u/bbg6470/athenapk/runs/root_ladder/fhc_rootladder.in
WRAP=$H/wrap_mod.sh
install -m 755 /beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh $WRAP
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

leg () {   # name  binary
  G=$H/$1; rm -rf $G; mkdir -p $G
  echo "=== $1 ==="; md5sum $2
  ( cd $G && stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $2 -i $DECK \
      parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
      parthenon/mesh/do_coalesced_comms=true \
      parthenon/time/nlim=40 \
      diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
      diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
      diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
      diffusion/cap_diag=true hydro/mag_diag=true \
      > run.log 2>&1 )
  echo "$1 exit=$?"
}

leg old /beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
leg new /beegfs/u/bbg6470/athenapk/build_gpu_v2/bin/athenaPK

echo "=== VERDICT ==="
if cmp -s $H/old/parthenon.out0.hst $H/new/parthenon.out0.hst; then
  echo "PASS: history files are BYTE-IDENTICAL -> f181c0a1 is an OFF-state no-op, safe to swap."
else
  echo "FAIL: history files DIFFER. Diff below; do NOT swap the production binary."
  diff $H/old/parthenon.out0.hst $H/new/parthenon.out0.hst | head -20
fi
echo "done $(date)"
