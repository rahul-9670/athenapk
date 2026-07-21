#!/bin/bash
#SBATCH --job-name=fused3way
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:h100:4
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/dt_attrib/%x_%j.out
set -o pipefail
# FUSED-IONIZATION_CHEM 3-WAY PROBE (2026-07-12): prices/validates extending the fused
# non-ideal evaluator (FusedNonidealEval) to the production ionization_chem mix.
#   A: athenaPK_eos_v2 (current production binary)          -> baseline
#   B: athenaPK_eos_v3, diffusion/fused_nonideal_dt=false   -> isolates the fused
#      PrecomputeNonidealEta cache fill (expected bit-identical eta: same pure
#      functions, same arguments; dt trace must match A exactly)
#   C: athenaPK_eos_v3, fused on (default)                  -> adds the FusedMixed dt
#      estimator (expected small dt shift: EOS-table T vs p/rho in the per-term
#      estimators -- documented; wsec_step should drop ~1-3 s/cycle)
# 8 cycles each on the NEWEST production restart, picked at RUN time (read-only; never
# delete it -- live chain resume point). All runs use the production changeover config
# (rkl2 mixed + freeze + cap 1000).
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN2=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_eos_v2
BIN3=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_eos_v3
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
WDIR=/beegfs/u/bbg6470/athenapk/runs/dt_attrib
WRAP=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
cd $WDIR

RESTART=$(ls -t /beegfs/u/bbg6470/athenapk/runs/prod_t4_full/parthenon.out2.*.rhdf | head -1)
NCYC=$($PY -c "import h5py; print(int(h5py.File('$RESTART','r')['Info'].attrs['NCycle']))")
NLIM=$((NCYC + 8))
echo "restart: $RESTART (cycle $NCYC) -> nlim $NLIM"
echo "binaries:"; md5sum $BIN2 $BIN3

COMMON="parthenon/time/nlim=$NLIM parthenon/output0/dt=1e9 parthenon/output1/dn=10000000 \
parthenon/output2/dn=10000000 parthenon/mesh/do_coalesced_comms=true \
diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true"

run_seg () { # $1=tag $2=binary $3=extra overrides
  mkdir -p $WDIR/f3w_$1
  echo "=== $1 start $(date) ==="
  stdbuf -oL -eL mpirun -n 4 $MCA $WRAP $2 -r $RESTART -d $WDIR/f3w_$1 $COMMON $3 \
    > $WDIR/log_f3w_$1.log 2>&1
  echo "=== $1 exit $? $(date) ==="
}

run_seg A_v2_baseline  $BIN2 ""
run_seg B_v3_fusedoff  $BIN3 "diffusion/fused_nonideal_dt=false"
run_seg C_v3_fusedon   $BIN3 ""

echo "---- dt traces (cycle dt wsec_step) ----"
for t in A_v2_baseline B_v3_fusedoff C_v3_fusedon; do
  echo "[$t]"
  grep '^cycle=' $WDIR/log_f3w_$t.log | \
    grep -o 'cycle=[0-9]*\|dt=[0-9.e+-]*\|wsec_step=[0-9.e+]*' | paste - - -
done
echo "---- STS ratios ----"
for t in A_v2_baseline B_v3_fusedoff C_v3_fusedon; do
  echo "[$t]"; grep -o 'STS ratio: [0-9.e+]* Taking [0-9]* steps' $WDIR/log_f3w_$t.log | tail -2
done
echo "=== fused3way done $(date) ==="
