#!/bin/bash
#SBATCH --job-name=d1meas
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --gres=gpu:h100:5
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#
# ##############################################################################################
# SUPERSEDED 2026-08-06 by submit_d1_leg.sh. DO NOT SUBMIT THIS SCRIPT -- it refuses to run.
#
# It is kept only as the record of job 2461421. It carries TWO independent defects, either of
# which alone produces a run with zero [D1] samples:
#   1. `parthenon/time/nlim=130`. nlim is ABSOLUTE and a restart RESTORES ncycle -- here 250.
#      KeepGoing() is `(time < tlim) && (nlim < 0 || ncycle < nlim)` (parthenon
#      basic_types.hpp:329), so 250 < 130 is false and the driver executes ZERO cycles. The
#      b7_closure runs that reached the OOM used nlim=500; submit_d1_leg.sh uses 380.
#   2. It pointed at build_gpu_v5, whose N3 dust guard aborts on this restart (see below).
# It also packs both legs into one 5-GPU / 4 h reservation, which forces the wider reservation
# and the longer limit onto the 4-rank leg for no reason -- the legs are independent.
# ##############################################################################################
echo "REFUSING TO RUN: submit_d1_meminfo.sh is superseded by submit_d1_leg.sh (see header)." >&2
echo "  NR=4 sbatch --ntasks=4 --gres=gpu:h100:4 --job-name=d1r4 --export=ALL,NR=4 submit_d1_leg.sh" >&2
echo "  NR=5 sbatch --ntasks=5 --gres=gpu:h100:5 --job-name=d1r5 --export=ALL,NR=5 submit_d1_leg.sh" >&2
exit 1
#
# D1 MEASUREMENT — kill or confirm the last standing hypothesis for the deep-AMR GPU OOM.
#
# THE STATE OF D1 (docs/validation/D1_gpu_memory_imbalance.md). Deep-AMR runs die on a few-MiB
# `bnd_flux::*` allocation while other cards have tens of GiB free: 49.3 / 56.6 / **79.2** GiB on
# IDENTICAL 198-block-per-rank assignments. Three hypotheses are already FALSIFIED and must not
# be re-derived:
#   1. block-distribution imbalance   -- Parthenon reports exactly equal counts.
#   2. the pinning wrapper aliases GPUs -- probe job 2454732: ranks 0-3 select devices 0-3.
#   3. "more ranks buys headroom"     -- 4 ranks and 5 ranks OOM at cycle 120, TO THE CYCLE.
#
# ONE HYPOTHESIS SURVIVES: the exhausted memory is AMR coarse/prolongation buffer space, sized by
# the mesh's COARSE-FINE BOUNDARY COUNT rather than by each rank's block count. It explains both
# facts at once -- the per-rank spread on equal block counts (Z-order intervals straddle different
# numbers of refinement boundaries) and the immovable failure cycle (rank count changes how blocks
# are distributed, not how many coarse-fine boundaries the hierarchy has). The failing label on
# the 5-rank run was literally `bnd_flux::cons.coarse`.
#
# THE MEASUREMENT. `<hydro> d1_meminfo=true` (default false => bit-identical when off) makes each
# rank print, whenever ITS block count changes, i.e. at every regrid and rebalance:
#     [D1] cycle=N rank=R nblocks=... coarse_fine_nbrs=... same_level_nbrs=... dev_free_GiB=...
# `coarse_fine_nbrs` counts that rank's neighbours at a DIFFERENT refinement level.
#
# THE DISCRIMINATOR, and why this job runs TWO legs. Within one leg, nblocks and coarse_fine_nbrs
# are correlated across ranks and a single leg cannot separate them. Changing the RANK COUNT
# breaks the degeneracy: blocks/rank falls by 1/nranks from 4 to 5, while the hierarchy's total
# coarse-fine boundary count is a property of the MESH and barely moves. So:
#     if consumed memory tracks nblocks           -> it must drop ~20 % per rank from 4 to 5
#     if it tracks coarse_fine_nbrs               -> it stays put, and so does the OOM cycle
# The second is what the falsified hypothesis 3 already implies; this job measures it directly
# instead of inferring it from a crash.
#
# THIS REQUIRES A POST-AUDIT-BATCH BINARY. The instrument landed in the 2026-08-05 audit batch
# (9f406ce) and exists in NO earlier binary -- which is why D1 was stuck: round 2 implemented the
# instrument but could not run it, because nothing had compiled the audit batch for GPU.
#
# 2026-08-06: MOVED OFF build_gpu_v5 (6d5b9895), which CANNOT run this job. v5 carries the
# original N3 guard -- a PARTHENON_REQUIRE_THROWS on the retired <dust> rho_unit_cgs/t_unit_cgs/
# T_unit_K keys existing. Those keys are not in any deck; `GetOrAddReal` INJECTS its default into
# the ParameterInput and Parthenon writes the whole list into every restart, so $RST (written
# 2026-08-04 by a pre-N3 binary) carries them tagged "# Default value added at run time" and the
# guard fires on a file no deck ever wrote. Job 2461421 died in 9 s, both legs, SLURM reporting
# COMPLETED 0:0 -- only the d1_samples==0 check below caught it. The guard is now a warning
# (src/dust/dust.cpp) and this points at that build. Do not point it back at v5.
#
# nlim=130 deliberately runs PAST the known cycle-120 failure. An OOM here is DATA, not a
# failure of the job: the whole point is to have per-rank memory logged right up to it.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1 OMPI_MCA_io=romio341
export TMPDIR=/beegfs/u/bbg6470/.chem_tmp; mkdir -p "$TMPDIR"
export LD_LIBRARY_PATH=/sw/env/gcc-13.3.0_openmpi-5.0.7/pkgsrc/2025Q1/lib:$LD_LIBRARY_PATH

H=/beegfs/u/bbg6470/athenapk/runs/deep_amr
RST=$H/run/parthenon.out2.00001.rhdf
# Pinned to the PRESERVED hard link, not build_gpu/bin/athenaPK: that path is the scratch build
# slot and the next `make -C build_gpu` would silently replace the binary this measurement ran on.
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK_PRESERVED_84a6d248
WRAP=$H/wrap_mod.sh
install -m 755 /beegfs/u/bbg6470/athenapk/runs/wrap_mod.sh $WRAP
MCA="--mca mtl ^psm2 --mca btl tcp,self,sm -x LD_LIBRARY_PATH -x PMIX_MCA_gds -x OMP_NUM_THREADS -x OMPI_MCA_io -x TMPDIR"

[ -x "$BIN" ] || { echo "MISSING $BIN -- build_gpu_v5 not built; aborting"; exit 1; }
[ -f "$RST" ] || { echo "MISSING restart $RST -- aborting"; exit 1; }
echo "job $SLURM_JOB_ID $(date)"; md5sum $BIN

leg () {   # nranks
  local N=$1
  local G=$H/d1_r$N; rm -rf $G; mkdir -p $G
  echo "=== leg nranks=$N ==="
  ( cd $G && stdbuf -oL -eL mpirun -n $N $MCA $WRAP $BIN -r $RST -t 01:45:00 \
      parthenon/time/nlim=130 \
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
  local RC=$?
  local NCYC=$(grep -c '^cycle=' $G/run.log)
  local ND1=$(grep -c '^\[D1\]' $G/run.log)
  echo "exit=$RC cycles=$NCYC d1_samples=$ND1"
  # POSITIVE CHECK. exit=0 is worthless here (this campaign has repeatedly seen SLURM report
  # COMPLETED for runs that died in seconds). The leg is only usable if the instrument fired.
  if [ "$ND1" -eq 0 ]; then
    echo "  **VOID: instrument produced no output. First error --**"
    grep -A4 -iE "PARTHENON ERROR|Kokkos ERROR|what\(\):" $G/run.log | head -8 | sed 's/^/     /'
  fi
  grep -iE "Kokkos ERROR|failed to allocate" $G/run.log | head -3 | sed 's/^/  OOM: /'
}

leg 4
leg 5

echo "=== D1 ANALYSIS ==="
/beegfs/u/bbg6470/venvs/analysis_env/bin/python \
  /beegfs/u/bbg6470/athenapk/runs/deep_amr/d1_analyse.py
echo "done $(date)"
