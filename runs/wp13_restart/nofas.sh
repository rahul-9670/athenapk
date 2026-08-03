#!/bin/bash
source ~/athenapk_env.sh >/dev/null 2>&1
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false PMIX_MCA_gds=hash
B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
R=/beegfs/u/bbg6470/athenapk/runs/wp13_restart/nofas
rm -rf $R; mkdir -p $R/{straight,split}
for d in straight split; do cp /beegfs/u/bbg6470/athenapk/runs/wp13_restart/straight/{fhc.in,units.json} $R/$d/; done
MC="self_gravity/solver_params/do_FAS=false"
env -C $R/straight $B -i $R/straight/fhc.in $MC parthenon/time/nlim=8 > $R/straight/run.log 2>&1
echo "straight exit=$?" >> $R/status
env -C $R/split $B -i $R/split/fhc.in $MC parthenon/time/nlim=6 > $R/split/run1.log 2>&1
echo "split1 exit=$?" >> $R/status
cp $R/split/parthenon.out2.final.rhdf $R/rst.rhdf
env -C $R/split $B -r $R/rst.rhdf $MC parthenon/time/nlim=8 > $R/split/run2.log 2>&1
echo "split2 exit=$?" >> $R/status
echo DONE >> $R/status
