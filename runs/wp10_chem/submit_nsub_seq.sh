#!/bin/bash
#SBATCH --job-name=nsubseq
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:h100:2
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=%x_%j.out
#
# WP-10 follow-up — nsub_max cost + physics on the PRODUCTION config, ALL THREE LEGS IN ONE JOB.
#
# WHY SEQUENTIAL IN ONE JOB, and not three concurrent jobs (the first design, jobs 2450143-45).
# The quantity being measured is COST (`wsec_step`), and three jobs landing on three different
# nodes at different times measure the machine as much as the parameter. WP-3 demonstrated exactly
# that failure: its cfl legs reported 1.02e6 vs 1.76e6 zone-cycles/wallsecond for runs that should
# be identical per zone-cycle, purely from node and contention differences -- which made their
# wall times unusable as a cost ratio. Running the three legs BACK TO BACK ON THE SAME TWO GPUs
# removes that confound entirely: same silicon, same neighbours, minutes apart.
#
# WHY 2 GPUs. It is what is actually free (g001 has 2 idle GPUs and 16 idle cores = exactly
# 2 x cpus-per-task; g004's 7 are unusable, the node is DRAINED with a RestrictedCoresPerGPU
# config fault). Rank count does not affect the physics (WP-12: decomposition-invariant to
# 0.000e+00 over 179 rows) and here it does not affect the COMPARISON either, because all three
# legs share it -- which is the entire point of putting them in one job.
#
# The earlier 3-GPU legs are superseded and cancelled. The completed `gpu_ns400` is KEPT: it
# doubled as an OFF-state gate (byte-identical to r128_sw on build_gpu_v2 over the full 89-cycle
# run), which is a result in its own right and is unaffected by this rerun.
#
# TIME LIMIT is deliberately 1.5 h, not the 6 h the root-ladder script defaults to. `ns400` took
# 4 m 27 s on 3 GPUs, so three legs on 2 GPUs is ~30 min even if nsub_max = 40000 costs 3x. An
# 80x over-request is not free: SLURM backfill can only place a job into a window as long as its
# limit, so a 6 h request sits behind gaps a 1.5 h request would drop straight into.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

H=/beegfs/u/bbg6470/athenapk/runs/wp10_chem
DECK=/beegfs/u/bbg6470/athenapk/runs/root_ladder/fhc_rootladder.in
BIN=/beegfs/u/bbg6470/athenapk/build_gpu_v3/bin/athenaPK
WRAP=$H/wrap_mod.sh
install -m 755 /beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh $WRAP
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
echo "job $SLURM_JOB_ID $(date)"; echo "node: $SLURM_JOB_NODELIST"
md5sum $BIN; md5sum $DECK

for NS in 400 4000 40000; do
  G=$H/seq_ns$NS; rm -rf $G; mkdir -p $G
  echo "=== nsub_max=$NS start $(date +%H:%M:%S) ==="
  ( cd $G && stdbuf -oL -eL mpirun -n 2 $MCA $WRAP $BIN -i $DECK -t 01:10:00 \
      parthenon/mesh/nx1=128 parthenon/mesh/nx2=128 parthenon/mesh/nx3=128 \
      parthenon/mesh/do_coalesced_comms=true \
      diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
      diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
      diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
      diffusion/cap_diag=true hydro/mag_diag=true \
      chemistry/nsub_max=$NS \
      > run.log 2>&1 )
  echo "nsub_max=$NS exit=$? $(date +%H:%M:%S)  $(grep -c '^cycle=' $G/run.log) cycles"
  grep -m1 -oE "WARNING Chemistry: [0-9]+ cell" $G/run.log | sed 's/^/    /'
done
echo "done $(date)"
