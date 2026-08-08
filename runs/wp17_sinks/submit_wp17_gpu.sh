#!/bin/bash
#SBATCH --job-name=wp17gpu
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:h100:4
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp17_sinks/wp17gpu_%j.out
set -o pipefail
#
# WS-1 increment 6 — the sink ACCRETION path on GPU.
#
# WHY. WP-17 closed the accretion gate with `Mdot / Mdot_B = 1.01 +/- 0.02` at N=128 (job
# 2453375), but on `build_cpu`, 32 MPI ranks. The GPU sinks path is described everywhere in this
# repo as "smoke-tested" only -- WP-17 itself says at line 188 that "the full-physics accretion
# path on GPU is not validated (WS-1 inc6), so that is not a free swap", and that is what blocks
# the WP-17 N=512 rung. This runs the IDENTICAL deck on GPU so the two answers are comparable.
#
# WHAT MAKES IT NON-VACUOUS. Attempt 1 of WP-17 (job 2452595) exited 0 on all three legs while
# testing nothing, because `sinks/accretion` defaults to false and the sink `mass` was absent from
# the dumps. Both are fixed in bondi_wp17.in, and the positive checks below are kept verbatim from
# submit_wp17b.sh: the sink mass must MOVE off its 2.5 seed value, and `mass` must be present in
# the dump. A GPU leg that runs clean but leaves the sink mass at 2.5 is a FAILURE, not a pass.
#
# CONFIGURATION, matched to the CPU leg it is being compared against:
#   nx = 128, meshblock = 32  -> 64 blocks (>= 4 ranks; Parthenon needs nblocks >= nranks)
#   racc_cells = 4            -> r_acc = 0.5, the value the CPU N=128 answer was measured at
#   tlim = 6.0                -> from the deck; 2.4 settling times (r_B/c = 2.5)
# 64 blocks over 4 ranks = 16 blocks/rank; at the measured 0.159 GiB/block (D1) that is ~2.5 GiB
# of 79.2 per card, so memory is not a constraint here.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

H=/beegfs/u/bbg6470/athenapk/runs/wp17_sinks
DECK=$H/bondi_wp17.in
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_PRESERVED_84a6d248
WRAP=$H/wrap_mod.sh
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

echo "job $SLURM_JOB_ID $(date) host=$(hostname)"; md5sum $BIN; md5sum $DECK
echo "analytic target: Mdot_B(seed M=2.5) = 4 pi * 0.625 * (G M)^2 rho / c^3"

N=128; M=32
G=$H/g$N; rm -rf $G; mkdir -p $G
echo "=== GPU nx1=$N meshblock=$M racc_cells=4 start $(date +%H:%M:%S) ==="
( cd $G && stdbuf -oL -eL mpirun -n 4 $MCA $WRAP $BIN -i $DECK \
    parthenon/mesh/nx1=$N parthenon/mesh/nx2=$N parthenon/mesh/nx3=$N \
    parthenon/meshblock/nx1=$M parthenon/meshblock/nx2=$M parthenon/meshblock/nx3=$M \
    sinks/racc_cells=4 \
    > run.log 2>&1 )
echo "gpu n=$N exit=$?  $(grep -c '^cycle=' $G/run.log) cycles  $(date +%H:%M:%S)"

# POSITIVE CHECKS -- verbatim from submit_wp17b.sh. Exit 0 is not evidence the test ran.
/beegfs/u/bbg6470/venvs/analysis_env/bin/python - "$G" <<'PY'
import sys, glob, os, h5py, numpy as np
fs = sorted(glob.glob(os.path.join(sys.argv[1], "parthenon.out0.*.phdf")))
if not fs: print("    NO DUMPS"); raise SystemExit
def m(fn):
    with h5py.File(fn,'r') as f:
        sv = f['sinks/SwarmVars']
        if 'mass' not in sv: return None, float(f['Info'].attrs['Time'])
        return float(np.array(sv['mass'])[0]), float(f['Info'].attrs['Time'])
m0,t0 = m(fs[0]); m1,t1 = m(fs[-1])
if m0 is None:
    print("    **`mass` ABSENT from dump -- swarm_variables did not take; rate unreadable.**")
else:
    print(f"    sink mass {m0:.6f} (t={t0:.2f}) -> {m1:.6f} (t={t1:.2f});"
          f"  mean dM/dt = {(m1-m0)/max(t1-t0,1e-30):.4f}"
          + ("   **NO ACCRETION -- knob inert**" if abs(m1-m0) < 1e-12 else ""))
PY

echo "=== GPU vs CPU comparison ==="
/beegfs/u/bbg6470/venvs/analysis_env/bin/python $H/wp17_gpu_compare.py
echo "done $(date)"
