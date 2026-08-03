#!/bin/bash
# WP-7 RE-BASELINE onto the WP-13 gravity fix. Staged 2026-08-02; launch with ./stage_gfix.sh go
#
# WHY. The first root ladder (r128/r256/r512, jobs 2434918/2434917/2438414) ran on binary
# `bffdf8cd`, which carries the stale-density Poisson RHS that WP-13 root-caused and fixed:
# SelfGravity::FillPoissonRHS read `prim` (Hydro's FillDerived output) with no guaranteed
# inter-package ordering, so gravity was solved from a one-stage-lagged density on EVERY step of
# EVERY self-gravity run. The fix (read `cons`, floored) is result-changing and is validated
# against the analytic Jeans growth rate: 2.870% -> 2.242% error in the best-resolved window.
#
# WP-7 is a convergence claim about a gravity-driven collapse, so it must be measured on the
# corrected operator. This re-run is ALSO the campaign's calibration: the size of the r256 shift
# between `bffdf8cd` and `49d9c257` decides whether the other stale-gravity results (WP-8 njeans,
# WP-18 sigma) need redoing or can stand as paired comparisons. Judge against WP-18's ~16% sigma.
#
# Compare against the OLD legs in r128/ r256/ r512/ -- they are kept, not overwritten.
set -euo pipefail
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK      # 49d9c257: gravity fix + all diags
HERE=/beegfs/u/bbg6470/athenapk/runs/root_ladder
WRAPSRC=/beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh   # md5 097f81b6, identical to what r128/
                                                      # r256/r512/r256_rank4 all ran with
[ -r "$WRAPSRC" ] || { echo "ABORT: GPU-pinning wrapper missing: $WRAPSRC"; exit 1; }

# `grep -qa` on the file, NOT `strings ... | grep -q`: under `set -o pipefail` grep -q exits on
# the first match, strings takes SIGPIPE, and the pipeline reports failure even on success.
# (Plain `strings` without -a also misses this symbol entirely in a 350 MB CUDA binary.)
grep -qa grav_rho_floor $BIN \
  || { echo "ABORT: $BIN predates the WP-13 gravity fix"; exit 1; }
echo "binary: $(md5sum $BIN)"

# NX  RUNDIR       NRANK  MBLK   (512^3 needs 64^3 blocks: 32^3 projects to 82.8 GiB > 80 GiB)
#
# RANK COUNT IS A FREE PARAMETER, packed to fit the free capacity on ONE node (MaxNodes=1).
# It is free because decomposition is measured not to affect results: WP-12 (r256 8-rank vs
# r256_rank4 4-rank, same deck) agreed to all 7 printed figures on every history quantity across
# 179 rows, and cycle-1 dt agrees to 11 significant figures. So there is nothing special about 8,
# or about powers of two -- 3 and 6 are equally valid and start NOW instead of queueing.
# Memory at 5.09 KiB/cell (conservative: that anchor includes AMR+multigrid overhead this uniform
# grid does not carry):  128^3 = 10.2 GiB total,  256^3 = 81.4 GiB,  512^3 = 651.5 GiB.
#   r128_gfix @ 3 ranks ->  3.4 GiB/GPU
#   r256_gfix @ 6 ranks -> 13.6 GiB/GPU
#   r512_gfix @ 8 ranks -> 81.4 GiB/GPU at 32^3 blocks (does NOT fit) -> MBLK=64 -> ~69.8 GiB.
# 512^3 is the ONE leg that cannot shrink: it barely fits on a full node and MaxNodes=1 forbids
# going wider. It will queue until a node frees; that is expected, not a fault.
# Block counts bound the rank count (256^3/32^3 = 512 blocks, 128^3/32^3 = 64); all fine here.
LEGS="128 r128_gfix 3 -
      256 r256_gfix 6 -
      512 r512_gfix 8 64"

echo "$LEGS" | while read NX DIR NRANK MBLK; do
  [ -z "${NX:-}" ] && continue
  mkdir -p $HERE/$DIR
  # The GPU-pinning wrapper lives at runs/wrap_mod.sh, NOT in root_ladder/ -- and it is stored
  # mode 0644 there, so it must be installed executable. The previous
  # `cp -n $HERE/wrap_mod.sh ... 2>/dev/null || true` pointed at a path that does not exist and
  # swallowed the error; mpirun then failed with "lacked permissions to execute" (jobs
  # 2446540/2446541, RUN_EXIT 183). `install -m 755` + a hard failure, so it can never be silent.
  install -m 755 $WRAPSRC $HERE/$DIR/wrap_mod.sh
  printf '%-12s NX=%-4s NRANK=%s MBLK=%s\n' "$DIR" "$NX" "$NRANK" "$MBLK"
  if [ "${1:-}" = "go" ]; then
    ENVS="NX=$NX,RUNDIR=$HERE/$DIR,NRANK=$NRANK,BIN=$BIN"
    [ "$MBLK" != "-" ] && ENVS="$ENVS,MBLK=$MBLK"
    # --ntasks/--gres MUST be passed here: submit_root.sh's #SBATCH header sets neither, so a
    # bare `sbatch submit_root.sh` requests ZERO GPUs and mpirun -n $NRANK then has no device.
    # The original ladder (jobs 2434917/2434918/2438414) supplied them on the sbatch CLI.
    sbatch --job-name=root${NX}g --nodes=1 --ntasks=$NRANK --gres=gpu:h100:$NRANK \
           --export=ALL,$ENVS $HERE/submit_root.sh
  fi
done

[ "${1:-}" = "go" ] || echo $'\nSTAGED ONLY. Re-run with:  ./stage_gfix.sh go'
