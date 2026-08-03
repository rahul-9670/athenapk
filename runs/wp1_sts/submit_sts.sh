#!/bin/bash
#SBATCH --job-name=wp1sts
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gres=gpu:h100:8
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp1_sts/slurm_%x_%j.out
set -o pipefail
#
# 2026-08-02 TAKEOVER FIX. The first four legs (jobs 2446363/64/65/66) ran ~55x too slow
# (wsec_step ~400 s, 4.5e4 zone-cycles/wsec) against root_ladder/r256's 5-12 s and 1.4-3.3e6
# on the BYTE-IDENTICAL deck. GPU load was 2-4% with CPU at 100% => host-bound, GPUs idle.
# Cause was NOT the binary (a bisect over three binaries returned the same 4.3e4): this script
# omitted `parthenon/mesh/do_coalesced_comms=true`, whose Parthenon default is FALSE
# (external/parthenon/src/mesh/mesh.cpp:96-97). Every other production submit script passes it
# (submit_root.sh:73, submit_flagship_test.sh:52, submit_conv.sh:67, submit_seed.sh:37).
# Measured A/B: runs/wp1_sts/attrib/RESULT_COALESCE (job 2446454).
# Also raised 4 -> 8 GPUs to match root_ladder/r256, so WP-1 legs are directly comparable to it.
source ~/athenapk_env.sh
module load cuda/12.5.1
R=/beegfs/u/bbg6470/athenapk/runs/wp1_sts/$LEG
B=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
W=/beegfs/u/bbg6470/athenapk/runs/wp1_sts/wrap_alloc.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"
export OMP_NUM_THREADS=1 PMIX_MCA_gds=hash OMPI_MCA_io=romio341
echo "LEG=$LEG OVERRIDES=$OV  BINARY: $(md5sum $B)"
# The binary MUST carry the WP-13 gravity stale-density fix, else the leg is measured on the
# superseded Poisson operator and is not comparable with anything else in the re-baselined
# campaign. `grav_rho_floor` is the fix's marker symbol.
# `grep -qa` on the file, not `strings | grep -q`: this script sets `-o pipefail`, and grep -q
# exiting early gives strings a SIGPIPE that would make the check fail on a GOOD binary.
if ! grep -qa grav_rho_floor $B; then
  echo "ABORT: $B predates the WP-13 gravity fix (no grav_rho_floor). Rebuild build_gpu."; exit 1
fi
# self-resuming: pick the newest restart if one exists
LATEST=$(ls -1t $R/*.rhdf 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then RA="-r $LATEST"; echo "resuming from $LATEST"; else RA="-i $R/fhc.in"; fi
cd $R && stdbuf -oL -eL mpirun -n 8 $MCA $W $B $RA \
  parthenon/mesh/do_coalesced_comms=true $OV > $R/run_$SLURM_JOB_ID.log 2>&1
echo "exit=$?" >> $R/status
# GUARD: the corrected k^2 turbulence sampler must have been used (WP-20).
grep -q "k sampling *: k2" $R/run_$SLURM_JOB_ID.log \
  || echo "IC_CHECK FAILED: k2 sampler absent -- results void" >> $R/status
