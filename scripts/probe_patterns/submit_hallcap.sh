#!/bin/bash
#SBATCH --job-name=hallcap
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/dt_attrib/%x_%j.out
set -o pipefail
# ETA_H CAP PRICING PROBE (2026-07-14): prices diffusion/eta_hall_cap_code (binary v5,
# md5 d145a7df) at the live second-collapse state (checkpoint 00163, cycle 77250).
# Companion to submit_hallprobe.sh (v4 A_ctl / B_nohall on the same restart).
#   A_v5_nocap: v5, cap disabled -> OFF-STATE CHECK: dt trace must be bit-identical
#               to hallprobe A_ctl (v4); any drift = the cap plumbing leaks.
#   B_cap03:    |eta_H| cap 0.3
#   C_cap01:    |eta_H| cap 0.1 (the eta_O-cap value; whistler dt ~ dx^2/|eta_H|)
# 8 cycles each on the SAME newest restart (read-only; chain resume point).
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN5=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_eos_v5
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
WDIR=/beegfs/u/bbg6470/athenapk/runs/dt_attrib
WRAP=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
cd $WDIR

# 2026-07-16: PINNED to the quarantined ckpt 00163 (cycle 77250). Post-bug outputs were
# quarantined (newest live restart is now 00142/cycle 71000), but the cap pricing and the
# off-state bit-identity check vs hallprobe (which already ran on 00163) need the SAME
# state. Binary stays v5 (bug-era STS, matching hallprobe's v4) so the off-state dt trace
# is comparable; the cap's RELATIVE dt gain transfers to v6.
RESTART=/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/quarantine_postbug_71000/parthenon.out2.00163.rhdf
NCYC=$($PY -c "import h5py; print(int(h5py.File('$RESTART','r')['Info'].attrs['NCycle']))")
NLIM=$((NCYC + 8))
echo "restart: $RESTART (cycle $NCYC) -> nlim $NLIM"
echo "binary:"; md5sum $BIN5

COMMON="parthenon/time/nlim=$NLIM parthenon/output0/dt=1e9 parthenon/output1/dn=10000000 \
parthenon/output2/dn=10000000 parthenon/mesh/do_coalesced_comms=true \
diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
diffusion/eta_ohm_cap_code=0.1"

run_seg () { # $1=tag $2=extra overrides
  mkdir -p $WDIR/hc_$1
  echo "=== $1 start $(date) ==="
  stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN5 -r $RESTART -d $WDIR/hc_$1 $COMMON $2 \
    > $WDIR/log_hc_$1.log 2>&1
  echo "=== $1 exit $? $(date) ==="
}

run_seg A_v5_nocap ""
run_seg B_cap03 "diffusion/eta_hall_cap_code=0.3"
run_seg C_cap01 "diffusion/eta_hall_cap_code=0.1"

echo "---- banners ----"
grep -H "eta_O cap\|eta_H| cap\|Hall effect" $WDIR/log_hc_*.log
echo "---- dt traces (cycle dt wsec_step) ----"
for t in A_v5_nocap B_cap03 C_cap01; do
  echo "[$t]"
  grep '^cycle=' $WDIR/log_hc_$t.log | \
    grep -o 'cycle=[0-9]*\|dt=[0-9.e+-]*\|wsec_step=[0-9.e+]*' | paste - - -
done
echo "---- STS ratios (last 3 each) ----"
for t in A_v5_nocap B_cap03 C_cap01; do
  echo "[$t]"; grep -o 'STS ratio: [0-9.e+]* Taking [0-9]* steps' $WDIR/log_hc_$t.log | tail -3
done
echo "---- NaN/error check ----"
grep -iH "nan\|fail\|abort" $WDIR/log_hc_*.log | grep -v "GetOrAdd" | head
echo "=== hallcap done $(date) ==="
