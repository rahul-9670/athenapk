#!/bin/bash
# GPU pinning that uses the allocation SLURM actually gave us.
#
# THE BUG THIS REPLACES. The stock wrap_mod.sh does:
#     NGPU=$(nvidia-smi -L | wc -l); export CUDA_VISIBLE_DEVICES=$(( local_rank % NGPU ))
# i.e. it derives the device from what is VISIBLE, and overwrites the correct value SLURM had
# already exported. `nvidia-smi` is not reliably restricted to the job's allocation on these
# nodes: on g004 it listed SIX devices for a FIVE-GPU allocation. NGPU is then 6, ranks 0-4 map
# to devices 0,1,2,3,4, and a rank can land on a card the job does not own.
#
# Measured, job 2495385 on g004: rank 2 aborted, and the pre-run table showed device index 2
# holding 35,975 MiB of another job's memory -- Kokkos then failed to allocate 17.09 MiB. Job
# 2495384 on g002 failed the other way, rank 4 with cudaGetDeviceCount ->
# cudaErrorInitializationError, i.e. a device we were not permitted to touch at all.
# Probe 2493185 confirms SLURM does export the right list: every rank saw
# CUDA_VISIBLE_DEVICES='0,1,2,3' for a 4-GPU allocation BEFORE the wrapper replaced it.
#
# So: index into SLURM's list, never recompute it. If SLURM gave us nothing (interactive /
# non-SLURM use), fall back to the old behaviour rather than failing.
LR=${OMPI_COMM_WORLD_LOCAL_RANK:-0}
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then
  IFS=',' read -ra DEVS <<< "$CUDA_VISIBLE_DEVICES"
  N=${#DEVS[@]}
  [ "$N" -gt 0 ] && export CUDA_VISIBLE_DEVICES="${DEVS[$(( LR % N ))]}"
else
  N=$(nvidia-smi -L 2>/dev/null | wc -l); [ "$N" -lt 1 ] && N=1
  export CUDA_VISIBLE_DEVICES=$(( LR % N ))
fi
[ "${WRAP_VERBOSE:-0}" = "1" ] && echo "rank $LR -> CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES" >&2
exec "$@"
