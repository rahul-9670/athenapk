#!/bin/bash
#SBATCH --job-name=ct_ot_amr
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:25:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/ct_coalesced_fix/%x_%j.out
set -o pipefail
source ~/athenapk_env.sh
export OMP_NUM_THREADS=1 PMIX_MCA_gds=hash OMPI_MCA_io=romio341
BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
WD=/beegfs/u/bbg6470/athenapk/runs/ct_coalesced_fix; cd $WD
echo "=== OFF (ref) $(date) ==="; rm -rf ot_off; mkdir ot_off
srun --mpi=pmix -n 8 $BIN -i ot_ct_amr.in parthenon/mesh/do_coalesced_comms=false -d ot_off > ot_off/log 2>&1
echo "OFF exit=$? cyc=$(grep -cE '^cycle=' ot_off/log) driver=$(grep -c 'Driver completed' ot_off/log)"
echo "=== ON (reproduce?) $(date) ==="; rm -rf ot_on; mkdir ot_on
srun --mpi=pmix -n 8 $BIN -i ot_ct_amr.in parthenon/mesh/do_coalesced_comms=true -d ot_on > ot_on/log 2>&1
echo "ON exit=$? cyc=$(grep -cE '^cycle=' ot_on/log) recv-buf=$(grep -c 'receiving buffer' ot_on/log) driver=$(grep -c 'Driver completed' ot_on/log)"
echo "=== DONE $(date) ==="
