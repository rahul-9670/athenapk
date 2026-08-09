#!/bin/bash
#SBATCH --job-name=gpumap
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:h100:4
#SBATCH --cpus-per-task=8
#SBATCH --time=00:05:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/wp13b_ab_deep/gpumap_%j.out
#
# WHY THIS EXISTS. Leg 2492769 (AB_LEG=new) died at 38 s with 78.1 GiB on ONE of its four H100s
# while two others sat at 9.9 GiB and the fourth at 0.0. I called that "an occupied card at init"
# and resubmitted with --exclusive. That diagnosis is NOT supported by the evidence I have:
#   * the RRZ epilogue is a NODE-level device table, so it cannot separate "another job's memory
#     on a shared card" from "two of my own ranks landed on the same card";
#   * a 4-way rank collision would not leave device 5 at exactly 0.0 GiB, so the collision story
#     does not fit either;
#   * the run.log holding the actual CUDA error was deleted when I cleaned $D for the resubmit.
# So the cause is UNKNOWN, and --exclusive is a fix aimed at one of two untested hypotheses. It
# also costs ~9 h of queue: 2492773 is Reason=Resources purely because no node is fully free,
# while g001 and g004 each have 5 idle H100s right now.
#
# WHAT THIS SETTLES, in ~1 min of 4 GPUs:
#   (a) does the SLURM cgroup MASK the GPUs? i.e. does `nvidia-smi -L` inside the allocation list
#       4 devices or all 8? wrap_mod.sh computes NGPU from exactly that command and then sets
#       CUDA_VISIBLE_DEVICES=local_rank % NGPU. If masking works, NGPU=4 and the four ranks get
#       four DISTINCT logical devices -- the wrapper is correct. If nvidia-smi reports 8, the
#       wrapper hands ranks physical ids 0..3, of which some are not ours: that is the bug.
#   (b) are the allocated cards ALREADY holding another job's memory at t=0? memory.used per
#       device, read before this job allocates anything, answers it directly.
# (a) and (b) are the two hypotheses, and they have opposite fixes: (a) is a wrapper bug that
# --exclusive would only paper over; (b) is real card sharing that --exclusive genuinely fixes.
source ~/athenapk_env.sh; module load cuda/12.5.1
export PMIX_MCA_gds=hash OMP_NUM_THREADS=1
R=/beegfs/u/bbg6470/athenapk/runs/wp13b_ab_deep
cat > $R/_probe_rank.sh <<'RANK'
#!/bin/bash
LR=${OMPI_COMM_WORLD_LOCAL_RANK:-0}
NGPU=$(nvidia-smi -L 2>/dev/null | wc -l)
# exactly what wrap_mod.sh would choose, without acting on it
WOULD=$(( LR % (NGPU<1?1:NGPU) ))
echo "rank=$LR  SLURM_CVD_in='${CUDA_VISIBLE_DEVICES:-<unset>}'  nvidia-smi-L-count=$NGPU  wrapper_would_set=$WOULD"
[ "$LR" = "0" ] && { echo "--- nvidia-smi -L as seen inside the allocation:"; nvidia-smi -L | sed 's/^/    /';
  echo "--- per-visible-device memory.used BEFORE we allocate anything:";
  nvidia-smi --query-gpu=index,uuid,memory.used,memory.total --format=csv,noheader | sed 's/^/    /'; }
RANK
chmod +x $R/_probe_rank.sh
echo "=== node $(hostname)  job $SLURM_JOB_ID  $(date) ==="
echo "SLURM_JOB_GPUS=${SLURM_JOB_GPUS:-<unset>}  SLURM_STEP_GPUS=${SLURM_STEP_GPUS:-<unset>}  CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"
mpirun -n 4 --mca mtl ^psm2 --mca btl tcp,self,sm $R/_probe_rank.sh
echo "=== node-wide view (outside any mask) via scontrol ==="
scontrol show job $SLURM_JOB_ID | grep -iE "gres|nodelist" 
echo PROBE_DONE
