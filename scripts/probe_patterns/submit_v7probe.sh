#!/bin/bash
#SBATCH --job-name=v7probe
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/probe_v7/%x_%j.out
set -o pipefail
# V7 PRE-SWAP PROBE (2026-07-20): pinned-restart A/B/C on the live production state
# (post-cleanup, chain slot 12 era), 30 cycles each, read-only on ckpt 00185.
#   A_v6ctl : production binary athenaPK_eos_v6 (03a299cf), exact production config.
#   B_v7def : candidate v7 (505491e8), same config, amr_check_interval UNSET (default 1).
#             The ONLY live code difference vs A on this config is audit fix F1
#             (chemistry reaction T from the EOS table under eos=hydrogen), so A-vs-B
#             *measures the F1 correction*, expected small at T_c~350 K (enters via
#             chem-x_e -> capped AD only). NOT expected bit-identical.
#   C_v7iv10: v7 + parthenon/mesh/amr_check_interval=10 -> prices the AMR-gating
#             speedup (wsec_AMR share) and its physics drift vs B.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/athenapk/runs/probe_v7/tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN_V6=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_eos_v6
BIN_V7=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
WDIR=/beegfs/u/bbg6470/athenapk/runs/probe_v7
WRAP=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
cd $WDIR

# Pinned EXPLICITLY (not newest-at-runtime): all three segments must read the same
# checkpoint while the chain keeps writing newer ones. Read-only; never delete.
RESTART=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/parthenon.out2.00185.rhdf
NCYC=$($PY -c "import h5py; print(int(h5py.File('$RESTART','r')['Info'].attrs['NCycle']))")
NLIM=$((NCYC + 30))
echo "restart: $RESTART (cycle $NCYC) -> nlim $NLIM"
echo "binaries:"; md5sum $BIN_V6 $BIN_V7

# Exact production CLI overrides (mirrors prod_t4_full/submit.sh), outputs disabled
# except the automatic final dump (the state-comparison artifact).
COMMON="parthenon/time/nlim=$NLIM parthenon/output0/dt=1e9 parthenon/output1/dn=10000000 \
parthenon/output2/dn=10000000 parthenon/mesh/do_coalesced_comms=true \
diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
diffusion/eta_ohm_cap_code=0.1"

run_seg () { # $1=tag $2=binary $3=extra overrides
  mkdir -p $WDIR/$1
  echo "=== $1 start $(date) ==="
  stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $2 -r $RESTART -d $WDIR/$1 $COMMON $3 \
    > $WDIR/log_$1.log 2>&1
  echo "=== $1 exit $? $(date) ==="
}

run_seg A_v6ctl  $BIN_V6 ""
run_seg B_v7def  $BIN_V7 ""
run_seg C_v7iv10 $BIN_V7 "parthenon/mesh/amr_check_interval=10"

echo "---- per-segment cycle stats (mean wsec_step / wsec_AMR over the 30 cycles) ----"
for t in A_v6ctl B_v7def C_v7iv10; do
  echo -n "[$t] "
  grep -oP 'wsec_step=\S+ zone-cycles/wsec=\S+ wsec_AMR=\S+' $WDIR/log_$t.log | \
    awk -F'[= ]' '{s+=$2; a+=$6; n++} END {if(n>0) printf "n=%d mean_step=%.2f mean_AMR=%.2f total=%.1f s\n", n, s/n, a/n, s+a; else print "NO CYCLES"}'
done
echo "---- final dt + time per segment ----"
for t in A_v6ctl B_v7def C_v7iv10; do
  echo -n "[$t] "; grep '^cycle=' $WDIR/log_$t.log | tail -1 | grep -oP 'cycle=\d+ time=\S+ dt=\S+'
done
echo "---- NaN/error check ----"
grep -iH "nan\|fail\|abort" $WDIR/log_*.log | grep -v "GetOrAdd" | head
echo "=== v7probe done $(date) ==="
