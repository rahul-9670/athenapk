#!/bin/bash
#SBATCH --job-name=ct_gpu_smoke
#SBATCH --account=banerjee_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:20:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/ct_tests/ct_gpu_smoke_%j.out
# Device validation of the CT path (Phase 2 increments 1-3): runs the 3D field loop
# (iprob=2, all E1/E2/E3 active) and the static-AMR field loop with GS05 CT on one H100,
# and prints the absolute face divergence. Expect round-off (~1e-16..1e-18) on-device.
set -o pipefail
source ~/athenapk_env.sh
module load cuda/12.5.1
export OMPI_MCA_mtl=^psm2 OMPI_MCA_btl=tcp,self,sm
BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
HERE=/beegfs/u/bbg6470/athenapk/runs/ct_tests
md5sum "$BIN"

maxcol(){ awk -v want="$2" '/^# \[/{for(i=1;i<=NF;i++)if($i~("\\]="want"$")){s=index($i,"[");e=index($i,"]");ci=substr($i,s+1,e-s-1)+0}next}/^#/{next}{if(ci>0&&($ci+0)>m)m=$ci+0}END{printf"%.4e",m}' $1; }

for cfg in "3d:field_loop_ct_3d.in:100" "amr:field_loop_ct_amr.in:140"; do
  name=${cfg%%:*}; rest=${cfg#*:}; inp=${rest%%:*}; nl=${rest##*:}
  dir="$HERE/gpu_${name}"; rm -rf "$dir"; mkdir -p "$dir"; cd "$dir"
  mpirun -n 1 "$BIN" -i "$HERE/$inp" parthenon/time/nlim="$nl" parthenon/output0/dt=1000 \
    >run.log 2>&1
  echo "GPU $name: EXIT $?  nb=$(awk '/^#/{next}{print $4;exit}' *.hst 2>/dev/null)  absDivB=$(maxcol *.hst ct_maxAbsDivB)"
done
echo "CT_GPU_SMOKE_DONE"
