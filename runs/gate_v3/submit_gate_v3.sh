#!/bin/bash
#SBATCH --job-name=gatev3
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# OFF-STATE GATE for build_gpu_v3.
#
#   old = build_gpu_v2/bin/athenaPK   f181c0a1   (current production)
#   new = build_gpu_v3/bin/athenaPK              (= f181c0a1 + B10 + B11)
#
# B10 is a diagnostic counter and a warning; B11 adds `hall_ohmic_floor_ratio`, which defaults
# to 0.0 and then makes EffectiveOhmicFloor() return the old constant. So with the production
# deck unchanged the two binaries must produce byte-identical history.
#
# THIS IS NOT REDUNDANT WITH THE CPU GATE. runs/wp16_cshock/gate_b11_ratio.sh established CPU
# bit-identity, which does NOT imply GPU bit-identity: the B11 edit changed the CONTROL FLOW
# inside device lambdas (`need_eta_h` now gates the eta_H evaluation), and nvcc is free to
# reassociate, re-fuse multiply-adds, or spill differently when a branch is restructured. Only
# running it on the H100 settles it.
#
# It also exercises what the CPU gate could not: the real production stack at 128^3 with AMR,
# multipole self-gravity, M1 radiation, chemistry and the mixed RKL2 non-ideal path.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

H=/beegfs/u/bbg6470/athenapk/runs/gate_v3
DECK=/beegfs/u/bbg6470/athenapk/runs/root_ladder/fhc_rootladder.in
WRAP=$H/wrap_mod.sh
install -m 755 /beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh $WRAP
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

leg () {   # name  binary  [extra CLI...]
  local G=$H/$1; shift; local BIN=$1; shift
  rm -rf $G; mkdir -p $G
  echo "=== $(basename $G) ==="; md5sum $BIN
  ( cd $G && stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN -i $DECK \
      parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
      parthenon/mesh/do_coalesced_comms=true parthenon/time/nlim=40 \
      diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
      diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
      diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
      diffusion/cap_diag=true hydro/mag_diag=true "$@" \
      > run.log 2>&1 )
  echo "exit=$?"
}

leg old /beegfs/u/bbg6470/athenapk/build_gpu_v2/bin/athenaPK
leg new /beegfs/u/bbg6470/athenapk/build_gpu_v3/bin/athenaPK
# B11 ON at the production floor: 0.2*max|eta_H| = 8.6e-3 < 0.05, so max() must pick the
# absolute floor in every cell and this must ALSO be byte-identical. If it is not, the
# production margin measured on prod_v9 does not hold at 128^3 and the ratio is NOT free here.
leg new_r02 /beegfs/u/bbg6470/athenapk/build_gpu_v3/bin/athenaPK diffusion/hall_ohmic_floor_ratio=0.2

echo "=== VERDICT ==="
for L in new new_r02; do
  if cmp -s $H/old/parthenon.out0.hst $H/$L/parthenon.out0.hst; then
    echo "PASS ($L vs old): history BYTE-IDENTICAL"
  else
    echo "FAIL ($L vs old): history DIFFERS"; diff $H/old/parthenon.out0.hst $H/$L/parthenon.out0.hst | head -8
  fi
done
echo "--- B10/B11 banners in the new binary ---"
grep -h "Hall Ohmic stabilizer\|NOTE \[Hall\]\|WARNING \[Hall\]\|WARNING Chemistry" $H/new_r02/run.log | head -4
echo "done $(date)"
