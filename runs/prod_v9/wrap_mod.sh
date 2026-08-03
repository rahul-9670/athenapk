#!/bin/bash
NGPU=$(nvidia-smi -L 2>/dev/null | wc -l); [ "$NGPU" -lt 1 ] && NGPU=1
export CUDA_VISIBLE_DEVICES=$(( ${OMPI_COMM_WORLD_LOCAL_RANK:-0} % NGPU ))
exec "$@"
