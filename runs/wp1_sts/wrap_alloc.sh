#!/bin/bash
# GPU pinning that respects the SLURM ALLOCATION.
#
# The production wrap_mod.sh does:  CUDA_VISIBLE_DEVICES = local_rank % (nvidia-smi -L | wc -l)
# nvidia-smi ignores CUDA_VISIBLE_DEVICES and always reports all 8 physical GPUs, so on a
# partially-allocated node that formula indexes GPUs the job does NOT own -> the job dies with
# cudaErrorInitializationError (job 2432856). That is why every GPU job so far had to request
# the whole node. This version indexes into the device list SLURM actually gave us, so a
# 2-GPU job is safe on a shared node.
LR=${OMPI_COMM_WORLD_LOCAL_RANK:-0}
if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
  IFS=',' read -ra DEVS <<< "$CUDA_VISIBLE_DEVICES"
else
  mapfile -t DEVS < <(nvidia-smi -L 2>/dev/null | awk '{print NR-1}')
fi
N=${#DEVS[@]}
if [ "$N" -lt 1 ]; then DEVS=(0); N=1; fi
export CUDA_VISIBLE_DEVICES=${DEVS[$(( LR % N ))]}
echo "rank(local)=$LR -> CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES (allocated: $N GPU(s))" >&2
exec "$@"
