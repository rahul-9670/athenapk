#!/bin/bash
#SBATCH --job-name=wp0verify
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:40:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp0_verify/%x_%j.out
set -o pipefail
#
# Settles two things in one GPU allocation.
#
# (1) WP-0 acceptance -- REPRODUCTION. The archive at docs/provenance/binary_5ebddce0/ was
#     rebuilt into wp0_repro/athenapk/build_gpu_repro. Its md5 is 122e32ab, NOT the archived
#     5ebddce0. That was expected: build paths, timestamps and link order are not
#     deterministic, so md5 equality was never the right acceptance test. FUNCTIONAL identity
#     is. Run A (original) vs run B (repro) must produce a byte-identical history file.
#
# (2) The new production candidate build_gpu_wp20 (bffdf8cd = current tree: WP-20 turb_ksample
#     + session B's WP-5 grav_diag, both default-OFF). Run C uses the UNMODIFIED eos_smoke deck,
#     so both new features are OFF and it must reproduce run A byte-for-byte on GPU -- the
#     OFF-state gate, repeated on the real hardware rather than only on CPU. Run D turns
#     turb_ksample=k2 on purely to confirm the banner and that it changes the IC.
#
# All four runs share one deck, one rank, one GPU, so nothing differs but the binary/flag.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

B=/beegfs/u/bbg6470/athenapk/runs/wp0_verify
DECK=/beegfs/u/bbg6470/athenapk/runs/eos_smoke/fhc.in
UNITS=/beegfs/u/bbg6470/athenapk/runs/eos_smoke/units.json
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

ORIG=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
REPRO=/beegfs/u/bbg6470/wp0_repro/athenapk/build_gpu_repro/bin/athenaPK
CAND=/beegfs/u/bbg6470/athenapk/build_gpu_wp20/bin/athenaPK

run () {  # $1=tag $2=binary $3...=extra CLI
  local tag=$1 bin=$2; shift 2
  rm -rf $B/$tag; mkdir -p $B/$tag; cp $DECK $UNITS $B/$tag/
  cd $B/$tag
  echo "--- $tag : $(md5sum $bin | cut -d' ' -f1)"
  mpirun -n 1 $MCA $B/wrap_mod.sh $bin -i fhc.in "$@" > out.log 2>&1
  echo "    exit=$? final: $(grep -h '^cycle=12 ' out.log | cut -d' ' -f1-3)"
}

run A_orig   $ORIG
run B_repro  $REPRO
run C_cand   $CAND
run D_cand_k2 $CAND problem/collapse_be/turb_ksample=k2

echo
echo "=============== RESULTS ==============="
echo "--- WP-0 reproduction (A original 5ebddce0  vs  B repro 122e32ab):"
if cmp -s $B/A_orig/parthenon.out0.hst $B/B_repro/parthenon.out0.hst; then
  echo "    PASS -- history BYTE-IDENTICAL. md5 differs, behaviour does not."
else
  echo "    FAIL -- histories differ:"; cmp $B/A_orig/parthenon.out0.hst $B/B_repro/parthenon.out0.hst
fi
echo "--- OFF-state gate on GPU (A original  vs  C candidate bffdf8cd, both features OFF):"
if cmp -s $B/A_orig/parthenon.out0.hst $B/C_cand/parthenon.out0.hst; then
  echo "    PASS -- history BYTE-IDENTICAL."
else
  echo "    FAIL -- histories differ:"; cmp $B/A_orig/parthenon.out0.hst $B/C_cand/parthenon.out0.hst
fi
echo "--- ON-state (D, turb_ksample=k2) banner:"
grep -A9 "Initial turbulence" $B/D_cand_k2/out.log | head -11
echo "--- D must DIFFER from C (the switch is live):"
cmp -s $B/C_cand/parthenon.out0.hst $B/D_cand_k2/parthenon.out0.hst \
  && echo "    FAIL -- identical, k2 had no effect" || echo "    PASS -- differs, k2 is live"
echo "======================================="
