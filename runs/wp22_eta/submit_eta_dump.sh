#!/bin/bash
#SBATCH --job-name=wp22dump
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#
# WP-22 — write out the APPLIED ambipolar diffusivity at first-core density.
#
# WHY THIS IS NEEDED. `wp22_eta_physical.py` reconstructs eta_A offline from (rho, B, x_e) using
# `Ionization::AmbipolarEtaFromXe`, and gets eta_A/eta_num ~ 1e4-1e5 at the first core. But that
# is only the CHEMISTRY branch. What the flux kernel actually applies is
#   eta_A = min(eta_chem, eta_eq, eta_ad_cap)          [diffusion.hpp:310-341]
# where eta_eq is the equilibrium NICIL/Wardle value used as a self-calibrating ceiling, added
# precisely because "in dense gas the chemistry x_e collapses (C->CO ionization minimum,
# x_e ~ 1e-13) which makes the single-fluid eta_A blow up". Dense gas IS the first core, so the
# ceiling may well bind exactly where the WP-22 claim lives. `eta_ad_cap_code` is NOT set in the
# production deck and defaults to numeric_limits<Real>::max() (hydro.cpp:1085), so that third term
# is inactive -- but eta_eq is not, and it cannot be reproduced offline without porting the whole
# grain + Saha charge solve.
#
# So the offline number is an UPPER BOUND on the applied eta_A, and an upper bound is the wrong
# direction for a claim that physical eta DOMINATES numerical eta. This job replaces it with the
# real thing.
#
# HOW. The code already caches the applied (eta_O, eta_H, eta_A) per cell in the Parthenon field
# `nonideal_eta` (hydro.cpp:1391-1395, components eta_ohmic/eta_hall/eta_ambipolar), registered
# whenever `diffusion/eta_cache` is on -- which it is in production. It is simply never written.
# Adding it to the output1 variable list makes the applied diffusivity a phdf field.
#
# SAFETY. prod_v9 is HELD (newest restart 2026-07-23) and is READ ONLY here: `-r` points at its
# restart, but the process runs with CWD in this directory, so every new file lands here and
# nothing in runs/prod_v9 is touched. Binary is prod_v9's OWN `athenaPK_eos_v9` (md5 17af621a),
# not the current production binary -- a restart must be read by the code that wrote it.
# nlim is set 2 cycles past the restart so it dumps and stops.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

H=/beegfs/u/bbg6470/athenapk/runs/wp22_eta
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_eos_v9
RST=/beegfs/u/bbg6470/athenapk/runs/prod_v9/parthenon.out2.00019.rhdf
WRAP=$H/wrap_mod.sh
install -m 755 /beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh $WRAP
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

G=$H/dump; rm -rf $G; mkdir -p $G
echo "job $SLURM_JOB_ID $(date)"; md5sum $BIN; ls -la $RST

# The restart's own cycle count is unknown here, so bound by WALL time (-t) instead of nlim and
# force a dump every cycle; the first dump is the one we want and anything after is harmless.
#
# WALL LIMIT MUST STAY SHORT. With `output1/dn=1` this writes a **5.67 GB** phdf EVERY CYCLE
# (2402 blocks x 32^3 x 18 components). Attempt 2 (2448571) was cancelled after 4 minutes having
# already written 24 GB; left to its 35-minute `-t` it would have produced ~350 GB and eaten the
# BeeGFS quota. Only the FIRST dump is needed. If this is ever re-run, drop `-t` to ~00:03:00 and
# delete the extras immediately (explicit per-file paths -- see CLAUDE.md, no tree-wide patterns).
#
# ATTEMPT 1 (2448542) ABORTED IN 38 s: `dt and dn are enabled for the same output block`
# (outputs.cpp:166). The deck drives BOTH hdf5 blocks by CYCLE COUNT -- output1 `dn = 50`,
# output2 `dn = 250` -- so adding a `dt` override left both set. Use `dn` for both instead of
# mixing the two. (The run DID get far enough to print the eta-cache banner, confirming
# `nonideal_eta` is registered on this binary.)
( cd $G && stdbuf -oL -eL mpirun -n 5 $MCA $WRAP $BIN -r $RST -t 00:03:00 \
    parthenon/output1/variables=prim,grav.phi,nonideal_eta \
    parthenon/output1/dn=1 \
    parthenon/output2/dn=1000000 \
    diffusion/cap_diag=true \
    > run.log 2>&1 )
echo "RUN_EXIT $? $(date)"
ls -la $G/*.phdf 2>/dev/null
grep -m1 "Non-ideal eta cache" $G/run.log || echo "WARNING: eta cache banner not found -- nonideal_eta may not exist"
