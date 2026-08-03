#!/bin/bash
#SBATCH --job-name=stsopt_ab
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=00:40:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/stsopt_ab/%x_%j.out
set -o pipefail
# A/B for the STS exchange-narrowing optimization. base = audit+units (no STS opt),
# cand = audit+units+STS opt. Both restart from a COPY of a prod_v8 first-core restart
# (cycle 1000, 1807 blocks) and run 8 cycles. Isolates the STS opt (units migration common
# to both). Result-preserving <=> base/cand agree to roundoff (a stale-ghost bug shows as a
# localized ~1e-2 boundary error). Also compares wsec_step.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
W=/beegfs/u/bbg6470/athenapk/runs/stsopt_ab; WRAP=$W/wrap_mod.sh
CLI="parthenon/mesh/do_coalesced_comms=true diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 parthenon/time/nlim=1008 parthenon/output1/dn=8 parthenon/output2/dn=1000000"
cd $W
for TAG in base cand; do
  BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_stsopt_$TAG
  OUT=$W/out_$TAG; rm -rf $OUT; mkdir -p $OUT
  echo "=== $TAG: $(md5sum $BIN) ==="
  stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN -r $W/seed.rhdf -d $OUT $CLI > $OUT/log.txt 2>&1
  echo "$TAG exit=$? $(grep '^cycle=' $OUT/log.txt | tail -1)"
  echo "$TAG NaN/abort: $(grep -ciE 'nan|abort|Inconsistent' $OUT/log.txt)"
  echo "$TAG med wsec_step: $(grep '^cycle=' $OUT/log.txt | grep -oE 'wsec_step=[0-9.e+-]+' | tail -6 | cut -d= -f2 | sort -n | sed -n 3p)"
done
echo "DONE $(date)"
