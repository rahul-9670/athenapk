#!/bin/bash
source ~/athenapk_env.sh >/dev/null 2>&1
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false PMIX_MCA_gds=hash
B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
R=/beegfs/u/bbg6470/athenapk/runs/wp13_rad
env -C $R/straight $B -i $R/straight/rad_shock.in > $R/straight/run.log 2>&1
echo "straight exit=$?" >> $R/status
env -C $R/split $B -i $R/split/rad_shock.in parthenon/time/nlim=6 > $R/split/run1.log 2>&1
echo "split1 exit=$?" >> $R/status
RST=$R/split/parthenon.out9.00001.rhdf
[ -f "$RST" ] || RST=$(ls -1 $R/split/*.rhdf | tail -1)
echo "restart_from=$RST" >> $R/status
env -C $R/split $B -r $RST parthenon/time/nlim=12 > $R/split/run2.log 2>&1
echo "split2 exit=$?" >> $R/status
echo DONE >> $R/status
