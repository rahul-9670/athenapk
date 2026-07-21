#!/bin/bash
#SBATCH --job-name=etacap
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:h100:4
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/dt_attrib/%x_%j.out
set -o pipefail
# ETA_O CAP PRICING PROBE (2026-07-12): prices diffusion/eta_ohm_cap_code at the live
# production state before it goes into the chain (user-approved physics change).
#   A: v4, no cap        -> must reproduce v3 behavior (cap default = disabled no-op)
#   B: v4, cap = 0.3     -> weak cap (price curve point)
#   C: v4, cap = 0.1     -> the approved production value
# 8 cycles each on the NEWEST production restart, picked at RUN time (read-only; never
# delete it -- live chain resume point). Production changeover config throughout
# (rkl2 mixed + freeze + cap 1000). Expect: dt_par ~ 1/eta_max rises with the cap; dt
# rises until the whistler/hydro constraint binds or the STS ratio unpins from 1000.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BINDIR=/beegfs/u/bbg6470/athenapk/build_gpu/bin
# Freeze the just-built eta-cap binary as v4 (build job 2344650 leaves it at bin/athenaPK).
cp -p $BINDIR/athenaPK $BINDIR/athenaPK_eos_v4
BIN4=$BINDIR/athenaPK_eos_v4
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
WDIR=/beegfs/u/bbg6470/athenapk/runs/dt_attrib
WRAP=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
cd $WDIR

RESTART=$(ls -t /beegfs/u/bbg6470/athenapk/runs/prod_t4_full/parthenon.out2.*.rhdf | head -1)
NCYC=$($PY -c "import h5py; print(int(h5py.File('$RESTART','r')['Info'].attrs['NCycle']))")
NLIM=$((NCYC + 8))
echo "restart: $RESTART (cycle $NCYC) -> nlim $NLIM"
echo "binary:"; md5sum $BIN4

COMMON="parthenon/time/nlim=$NLIM parthenon/output0/dt=1e9 parthenon/output1/dn=10000000 \
parthenon/output2/dn=10000000 parthenon/mesh/do_coalesced_comms=true \
diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true"

run_seg () { # $1=tag $2=extra overrides
  mkdir -p $WDIR/ecap_$1
  echo "=== $1 start $(date) ==="
  stdbuf -oL -eL mpirun -n 4 $MCA $WRAP $BIN4 -r $RESTART -d $WDIR/ecap_$1 $COMMON $2 \
    > $WDIR/log_ecap_$1.log 2>&1
  echo "=== $1 exit $? $(date) ==="
}

run_seg A_nocap ""
run_seg B_cap03 "diffusion/eta_ohm_cap_code=0.3"
run_seg C_cap01 "diffusion/eta_ohm_cap_code=0.1"

echo "---- cap banners ----"
grep -H "eta_O cap" $WDIR/log_ecap_*.log
echo "---- dt traces (cycle dt wsec_step) ----"
for t in A_nocap B_cap03 C_cap01; do
  echo "[$t]"
  grep '^cycle=' $WDIR/log_ecap_$t.log | \
    grep -o 'cycle=[0-9]*\|dt=[0-9.e+-]*\|wsec_step=[0-9.e+]*' | paste - - -
done
echo "---- STS ratios (last 3 each) ----"
for t in A_nocap B_cap03 C_cap01; do
  echo "[$t]"; grep -o 'STS ratio: [0-9.e+]* Taking [0-9]* steps' $WDIR/log_ecap_$t.log | tail -3
done
echo "---- NaN/error check ----"
grep -iH "nan\|fail\|abort" $WDIR/log_ecap_*.log | grep -v "GetOrAdd" | head
echo "=== etacap done $(date) ==="
