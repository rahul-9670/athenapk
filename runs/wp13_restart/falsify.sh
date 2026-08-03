#!/bin/bash
source ~/athenapk_env.sh >/dev/null 2>&1
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false PMIX_MCA_gds=hash
B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
R=/beegfs/u/bbg6470/athenapk/runs/wp13_restart/falsify
rm -rf $R; mkdir -p $R/{fresh,rst}
cp /beegfs/u/bbg6470/athenapk/runs/wp13_restart/straight/{fhc.in,units.json} $R/fresh/
# radiation OFF: if the restart STILL diverges, radiation is not the cause -- the lost
# barotropic cooling is. If it matches, radiation was the proximate blow-up path.
env -C $R/fresh $B -i $R/fresh/fhc.in physics/radiation=false parthenon/time/nlim=7 \
    parthenon/output2/dn=6 > $R/fresh/run.log 2>&1
echo "fresh exit=$?" >> $R/status
cp $R/fresh/parthenon.out2.00001.rhdf $R/ 2>/dev/null
env -C $R/rst $B -r $R/parthenon.out2.00001.rhdf physics/radiation=false \
    parthenon/time/nlim=7 > $R/rst/run.log 2>&1
echo "rst exit=$?" >> $R/status
echo DONE >> $R/status
