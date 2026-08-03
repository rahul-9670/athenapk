#!/bin/bash
# Per-rank GPU pin (mirror of prod wrap_mod.sh) + nsys profile wrapper.
NGPU=$(nvidia-smi -L 2>/dev/null | wc -l); [ "$NGPU" -lt 1 ] && NGPU=1
export CUDA_VISIBLE_DEVICES=$(( ${OMPI_COMM_WORLD_LOCAL_RANK:-0} % NGPU ))
R=${OMPI_COMM_WORLD_RANK:-0}
OUT=/beegfs/u/bbg6470/athenapk/runs/gpu_profile_v8/nsys/profile_r${R}
exec nsys profile \
  --trace=cuda,nvtx,mpi \
  --sample=none --cpuctxsw=none \
  --cuda-memory-usage=false \
  --force-overwrite=true \
  -o "$OUT" \
  "$@"
