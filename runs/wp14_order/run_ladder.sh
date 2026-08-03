#!/bin/bash
# WP-14 ladder: for each wave family, refine 16->32->64->128 (nx1; nx2=nx3=nx1/2) and read
# the L1 error that linear_wave_mhd::UserWorkAfterLoop writes. Observed order p between
# successive rungs is log2(L1_coarse / L1_fine).
source ~/athenapk_env.sh >/dev/null 2>&1
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false PMIX_MCA_gds=hash
B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/wp14_order
for WF in ${WAVES:-0 1 2 3}; do
  for N in 16 32 64 128; do
    G=$H/wf${WF}_n${N}; rm -rf $G; mkdir -p $G
    M=16; [ $N -le 16 ] && M=8
    env -C $G $B -i $H/lw_mhd.in \
      problem/linear_wave/wave_flag=$WF \
      parthenon/mesh/nx1=$((N*2)) parthenon/mesh/nx2=$N parthenon/mesh/nx3=$N \
      parthenon/meshblock/nx1=$M parthenon/meshblock/nx2=$M parthenon/meshblock/nx3=$M \
      > $G/run.log 2>&1
    echo "wf=$WF n=$N exit=$?"
  done
done
