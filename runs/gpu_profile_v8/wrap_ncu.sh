#!/bin/bash
# Per-rank GPU pin; run ncu ONLY on rank 0 (serialized kernel replay), others run bare.
# Other ranks will block at MPI collectives while rank0 is slow under ncu -- expected.
NGPU=$(nvidia-smi -L 2>/dev/null | wc -l); [ "$NGPU" -lt 1 ] && NGPU=1
export CUDA_VISIBLE_DEVICES=$(( ${OMPI_COMM_WORLD_LOCAL_RANK:-0} % NGPU ))
R=${OMPI_COMM_WORLD_RANK:-0}
if [ "$R" = "0" ]; then
  OUT=/beegfs/u/bbg6470/athenapk/runs/gpu_profile_v8/ncu/ncu_rep
  exec ncu \
    --target-processes application-only \
    --nvtx \
    --set basic \
    --launch-skip ${NCU_SKIP:-300} \
    --launch-count ${NCU_COUNT:-40} \
    --force-overwrite \
    -o "$OUT" \
    "$@"
else
  exec "$@"
fi
