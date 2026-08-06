#!/bin/bash
#SBATCH --job-name=d1leg
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# D1 MEASUREMENT, ONE LEG. Submit with the rank count on the command line, e.g.
#     sbatch --ntasks=4 --gres=gpu:h100:4 --job-name=d1r4 submit_d1_leg.sh
#     sbatch --ntasks=5 --gres=gpu:h100:5 --job-name=d1r5 submit_d1_leg.sh
# The rank count is taken from SLURM_NTASKS (see below), so --export is NOT required and there is
# no wrapper script. (An earlier version of this header referred to a `drive_d1.sh` that does the
# submitting; no such file has ever existed. Corrected 2026-08-06.)
#
# WHY ONE LEG PER JOB rather than the previous both-legs-in-one-4h-reservation
# (submit_d1_meminfo.sh, superseded). The two legs are independent measurements that get
# compared afterwards by d1_analyse.py -- nothing in leg 5 depends on leg 4 having run. Packing
# them together bought nothing and cost queue position twice over: it forced the WIDER of the two
# reservations (5 GPUs) onto BOTH legs, and it forced a 4 h limit when a leg needs ~1.5 h. Split,
# the 4-GPU leg backfills on a narrower hole and neither leg blocks the other.
#
# THE TIME LIMIT IS MEASURED, NOT GUESSED. b7_closure2 (job 2454734), the same restart at 5 ranks,
# went cycle=250 -> 369 with wsec_total 2.55 -> 5.34e+03, i.e. the OOM arrives 89 min after start.
# --time=02:00:00 covers that with ~35 % margin; athenaPK gets -t 01:45:00 so it stops itself
# cleanly before SLURM kills the step.
#
# ---------------------------------------------------------------------------------------------
# THE MEASUREMENT (unchanged from submit_d1_meminfo.sh -- see docs/validation/
# D1_gpu_memory_imbalance.md). Deep-AMR runs die on a few-MiB `bnd_flux::*` allocation while other
# cards have tens of GiB free: 49.3 / 56.6 / 79.2 GiB on IDENTICAL 198-block-per-rank assignments.
# Three hypotheses are FALSIFIED and must not be re-derived:
#   1. block-distribution imbalance     -- Parthenon reports exactly equal counts.
#   2. the pinning wrapper aliases GPUs -- probe job 2454732: ranks 0-3 select devices 0-3.
#   3. "more ranks buys headroom"       -- 4 and 5 ranks OOM at cycle 369, TO THE CYCLE.
#
# ONE HYPOTHESIS SURVIVES: the exhausted memory is AMR coarse/prolongation buffer space, sized by
# the mesh's COARSE-FINE BOUNDARY COUNT rather than by each rank's block count. `<hydro>
# d1_meminfo=true` (default false => bit-identical when off) makes each rank print, whenever ITS
# block count changes (every regrid/rebalance):
#     [D1] cycle=N rank=R nblocks=... coarse_fine_nbrs=... same_level_nbrs=... dev_free_GiB=...
#
# THE DISCRIMINATOR needs BOTH legs. Within one leg nblocks and coarse_fine_nbrs are correlated
# across ranks and cannot be separated. Changing the RANK COUNT breaks the degeneracy: blocks/rank
# falls ~20 % from 4 to 5, while the hierarchy's total coarse-fine boundary count is a property of
# the MESH and barely moves. So 5 GPUs is the FLOOR for the second leg -- shrinking it is not a
# scheduling tweak, it deletes the experiment.
#
# nlim=380 deliberately runs PAST the known cycle-369 failure. An OOM here is DATA, not a failed
# job: the point is per-rank memory logged right up to it.
#
# nlim IS ABSOLUTE, AND THIS IS WHERE THE PREVIOUS ATTEMPT DIED A SECOND WAY. Parthenon's
# KeepGoing() is `(time < tlim) && (nlim < 0 || ncycle < nlim)` (basic_types.hpp:329) and a
# restart RESTORES ncycle -- here 250. submit_d1_meminfo.sh passed nlim=130, i.e. 250 < 130 is
# false, so it would have executed ZERO cycles and reported zero [D1] samples even after the
# dust-guard fix that killed job 2461421. The working b7_closure runs used nlim=500. Any change
# to this number must be >= 250 + (cycles wanted).
#
# THE BINARY IS PINNED TO A PRESERVED HARD LINK, not build_gpu/bin/athenaPK: that path is the
# scratch build slot and the next `make -C build_gpu` would silently replace the binary this
# measurement ran on. It must NOT be build_gpu_v5 (6d5b9895): v5 carries the original N3
# PARTHENON_REQUIRE_THROWS on the retired <dust> unit keys, and $RST (written 2026-08-04 by a
# pre-N3 binary) carries those keys because GetOrAddReal injects its defaults into the
# ParameterInput and Parthenon writes the whole list into every restart. Job 2461421 died in 9 s,
# both legs, with SLURM reporting COMPLETED 0:0.
set -u
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

H=/beegfs/u/bbg6470/athenapk/runs/deep_amr
RST=$H/run/parthenon.out2.00001.rhdf
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_PRESERVED_84a6d248
WRAP=$H/wrap_mod.sh
# Rank count comes from SLURM's own --ntasks, which is guaranteed to be set in a batch job and
# always agrees with the --gres the job actually holds. NR remains an override for interactive use.
# It used to be `NR=${NR:?...}`, i.e. wholly dependent on `--export=ALL,NR=N` having landed -- an
# assumption that cannot be checked from outside the job (scontrol does not print the submit
# environment), and whose failure mode is a 9-second death that costs a whole scheduling day.
# Deriving it from SLURM_NTASKS removes the assumption instead of documenting it.
NR=${NR:-${SLURM_NTASKS:?neither NR nor SLURM_NTASKS is set}}
if [ -n "${SLURM_NTASKS:-}" ] && [ "$NR" != "$SLURM_NTASKS" ]; then
  echo "REFUSING: NR=$NR disagrees with SLURM_NTASKS=$SLURM_NTASKS (GPUs are allocated per task)." >&2
  exit 1
fi
install -m 755 /beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh $WRAP
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

[ -x "$BIN" ] || { echo "MISSING $BIN -- aborting"; exit 1; }
[ -f "$RST" ] || { echo "MISSING restart $RST -- aborting"; exit 1; }
echo "job $SLURM_JOB_ID  nranks=$NR  $(date)"; md5sum $BIN

G=$H/d1_r$NR; rm -rf $G; mkdir -p $G
( cd $G && stdbuf -oL -eL mpirun -n $NR $MCA $WRAP $BIN -r $RST -t 01:45:00 \
    parthenon/time/nlim=380 \
    parthenon/mesh/do_coalesced_comms=true \
    diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 \
    diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true \
    diffusion/eta_ohm_cap_code=0.1 diffusion/ion_zeta=1.0e-16 \
    diffusion/hall_ohmic_floor_ratio=0.2 \
    hydro/d1_meminfo=true \
    hydro/boundary_flux_clamp=false \
    parthenon/output2/dn=1000000 \
    parthenon/output0/dt=-1.0 parthenon/output0/dn=1 \
    > run.log 2>&1 )
RC=$?
LASTCYC=$(grep -oE '^cycle=[0-9]+' $G/run.log | tail -1 | cut -d= -f2)
ND1=$(grep -c '^\[D1\]' $G/run.log)
echo "exit=$RC last_cycle=${LASTCYC:-none} d1_samples=$ND1"
# POSITIVE CHECK. exit=0 is worthless here -- job 2461421 was SLURM-COMPLETED 0:0 after dying at
# parse in 9 s, and only this check caught it. The leg is usable ONLY if the instrument fired.
if [ "$ND1" -eq 0 ]; then
  echo "  **VOID: instrument produced no output. First error --**"
  grep -A4 -iE "PARTHENON ERROR|Kokkos ERROR|what\(\):" $G/run.log | head -8 | sed 's/^/     /'
  # A run that stops at the restart cycle executed nothing: nlim/tlim already satisfied.
  [ "${LASTCYC:-0}" = "250" ] && echo "     ^ stopped AT the restart cycle 250 => nlim/tlim gate, not a crash."
fi
grep -iE "Kokkos ERROR|failed to allocate" $G/run.log | head -3 | sed 's/^/  OOM: /'
echo "done $(date)"
