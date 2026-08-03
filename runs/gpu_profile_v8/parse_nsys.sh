#!/bin/bash
# Run AFTER submit_nsys.sh (job 2377882) completes and writes nsys/profile_r*.nsys-rep.
# Emits the measured tables the GPU_OPTIMIZATION.md "PENDING" boxes need.
source ~/athenapk_env.sh; module load cuda/12.5.1
cd /beegfs/u/bbg6470/athenapk/runs/gpu_profile_v8/nsys || exit 1
for rep in profile_r*.nsys-rep; do
  [ -f "$rep" ] || continue
  base=${rep%.nsys-rep}
  echo "################ $rep ################"
  # Per-kernel GPU time + call counts (NVTX labels appear via the kp_nvtx connector)
  echo "==== gpukernsum (top kernels by GPU time) ===="
  nsys stats -q --report cuda_gpu_kern_sum --format table "$rep" 2>/dev/null | head -40
  # NVTX range summary: physics-labeled totals (STS/rad/gravity/eta) + subcycle counts
  echo "==== nvtxsum (labeled range totals + instance counts) ===="
  nsys stats -q --report nvtx_sum --format table "$rep" 2>/dev/null | head -60
  # CUDA API: launches, device/stream/event syncs, malloc/free churn
  echo "==== cudaapisum (API time: launches, syncs, malloc/free) ===="
  nsys stats -q --report cuda_api_sum --format table "$rep" 2>/dev/null | head -40
  # GPU memory ops (H2D/D2H, memset)
  echo "==== gpumemtimesum ===="
  nsys stats -q --report cuda_gpu_mem_time_sum --format table "$rep" 2>/dev/null | head -20
done
echo "### radiation subcycle count = (CalculateRadFluxes instances)/(#hydro steps profiled)"
echo "### per-rank load balance = compare total kernel GPU time across profile_r0..r4"
