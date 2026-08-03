#!/bin/bash
# GPU pinning wrapper (canonical copy -- WP-12, 2026-07-31).
#
# BUG FIXED HERE. The previous version was unconditionally:
#     NGPU=$(nvidia-smi -L | wc -l)
#     export CUDA_VISIBLE_DEVICES=$(( ${OMPI_COMM_WORLD_LOCAL_RANK:-0} % NGPU ))
# which OVERWRITES the device list SLURM already assigned to the task. That is
# only harmless when the job owns every GPU on the node starting at device 0.
# On a PARTIALLY allocated node the modulo indexes GPUs the job does not own:
# job 2432856 (4 of g002's 8 GPUs, another job holding 2) died in 5 s with
#     cudaGetDeviceCount(&count) error( cudaErrorInitializationError )
# from Kokkos_Core.cpp:135.
#
# Correct behaviour: SLURM already tells us which devices we may use, via
# CUDA_VISIBLE_DEVICES. Pick this local rank's entry FROM THAT LIST. Two cases,
# both handled:
#   * per-JOB binding (--gres=gpu:N, no per-task binding): the list holds all N
#     allocated devices, and LR % N gives each rank a distinct one.
#   * per-TASK binding (--gpus-per-task): the list holds a single device, and
#     LR % 1 == 0 selects it.
# Only when SLURM provided nothing do we fall back to the old node-wide count.
LR=${OMPI_COMM_WORLD_LOCAL_RANK:-${SLURM_LOCALID:-0}}
if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
  IFS=',' read -r -a DEVS <<< "$CUDA_VISIBLE_DEVICES"
  N=${#DEVS[@]}
  [ "$N" -lt 1 ] && N=1
  export CUDA_VISIBLE_DEVICES=${DEVS[$(( LR % N ))]}
else
  N=$(nvidia-smi -L 2>/dev/null | wc -l); [ "$N" -lt 1 ] && N=1
  export CUDA_VISIBLE_DEVICES=$(( LR % N ))
fi
exec "$@"
