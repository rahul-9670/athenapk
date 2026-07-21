#!/bin/bash
#SBATCH --job-name=creduc
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/dt_attrib/%x_%j.out
set -o pipefail
# RSLA CONVERGENCE PROBE (2026-07-18): does creduc=1000 throttle the first core's
# radiative cooling at tau~20-40 (warm bias)? 40-cycle restarts of the corrected redo
# state (ckpt 00160, cycle ~75500, rho_c 2.7e-9, T_c 276 K), binary v6:
#   A_c1000: production creduc=1000 (chat_code ~1578)
#   B_c100 : creduc=100 (chat_code ~15779, 10x more RT subcycles)
# Discriminator: L(r)=4pi r^2 <F_r> luminosity profile + central T drift + Er/aT^4.
# If A ~= B => RSLA error negligible at this state; if B cools faster => quantify bias.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_eos_v6
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
WDIR=/beegfs/u/bbg6470/athenapk/runs/dt_attrib
WRAP=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
cd $WDIR

RESTART=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/parthenon.out2.00160.rhdf
NCYC=$($PY -c "import h5py; print(int(h5py.File('$RESTART','r')['Info'].attrs['NCycle']))")
NLIM=$((NCYC + 40))
echo "restart: $RESTART (cycle $NCYC) -> nlim $NLIM"
echo "binary:"; md5sum $BIN

COMMON="parthenon/time/nlim=$NLIM parthenon/output0/dt=1e9 parthenon/output1/dn=20 \
parthenon/output2/dn=10000000 parthenon/mesh/do_coalesced_comms=true \
diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
diffusion/eta_ohm_cap_code=0.1"

run_seg () { # $1=tag $2=extra overrides
  mkdir -p $WDIR/cr_$1
  echo "=== $1 start $(date) ==="
  stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN -r $RESTART -d $WDIR/cr_$1 $COMMON $2 \
    > $WDIR/log_cr_$1.log 2>&1
  echo "=== $1 exit $? $(date) ==="
}

run_seg A_c1000 ""
run_seg B_c100 "radiation/creduc=100"

echo "---- dt/wsec tails ----"
for t in A_c1000 B_c100; do
  echo "[$t]"; grep '^cycle=' $WDIR/log_cr_$t.log | tail -3
done
echo "---- NaN/error check ----"
grep -iH "nan\|fail\|abort" $WDIR/log_cr_*.log | grep -v "GetOrAdd" | head
echo "=== creduc probe done $(date) ==="
