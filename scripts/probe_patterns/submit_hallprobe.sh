#!/bin/bash
#SBATCH --job-name=hallprobe
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/dt_attrib/%x_%j.out
set -o pipefail
# HALL DT-ATTRIBUTION PROBE (2026-07-14): measures the Hall/whistler share of the
# production dt at the live second-collapse state (rho ~ 1.1e6 x rhocrit, L14),
# before deciding on an eta_H cap (the eta_O-cap sequel). User directed: chain
# stopped at checkpoint 00163 for this retooling.
#   A_ctl:    exact production config (v4 + eta_O cap 0.1, mixed rkl2+Hall)
#   B_nohall: diffusion/hall=none  (whistler limit AND Ohmic floor both gone;
#             eta_O cap + AD unchanged) -> dt_B/dt_A = what Hall costs
# 8 cycles each on the NEWEST production restart, picked at RUN time (read-only;
# never delete it -- chain resume point). 5 GPUs: L14 footprint ~335 GiB no longer
# fits on 4 cards.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
export OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
BIN4=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_eos_v4
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
diffusion/integrator=rkl2 diffusion/rkl2_max_dt_ratio=1000 \
diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1"

run_seg () { # $1=tag $2=extra overrides
  mkdir -p $WDIR/hp_$1
  echo "=== $1 start $(date) ==="
  stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN4 -r $RESTART -d $WDIR/hp_$1 $COMMON $2 \
    > $WDIR/log_hp_$1.log 2>&1
  echo "=== $1 exit $? $(date) ==="
}

run_seg A_ctl    "diffusion/hall_floor_integrator=rkl2"
run_seg B_nohall "diffusion/hall=none"

echo "---- banners ----"
grep -H "eta_O cap\|Mixed diffusion\|Hall effect" $WDIR/log_hp_*.log
echo "---- dt traces (cycle dt wsec_step) ----"
for t in A_ctl B_nohall; do
  echo "[$t]"
  grep '^cycle=' $WDIR/log_hp_$t.log | \
    grep -o 'cycle=[0-9]*\|dt=[0-9.e+-]*\|wsec_step=[0-9.e+]*' | paste - - -
done
echo "---- STS ratios (last 3 each) ----"
for t in A_ctl B_nohall; do
  echo "[$t]"; grep -o 'STS ratio: [0-9.e+]* Taking [0-9]* steps' $WDIR/log_hp_$t.log | tail -3
done
echo "---- NaN/error check ----"
grep -iH "nan\|fail\|abort" $WDIR/log_hp_*.log | grep -v "GetOrAdd" | head
echo "=== hallprobe done $(date) ==="
