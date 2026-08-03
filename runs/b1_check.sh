#!/bin/bash
source ~/athenapk_env.sh >/dev/null 2>&1
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false PMIX_MCA_gds=hash
B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
R=/beegfs/u/bbg6470/athenapk/runs
for MODE in outflow diode; do
  G=$R/b1chk_$MODE; rm -rf $G; mkdir -p $G; cp $R/b2b4_gate/fhc.in $G/
  BCS=""; for f in ix1 ox1 ix2 ox2 ix3 ox3; do BCS="$BCS parthenon/mesh/${f}_bc=$MODE"; done
  env -C $G $B -i $G/fhc.in parthenon/time/nlim=3 hydro/cons_diag=true $BCS > $G/run.log 2>&1
  echo "$MODE exit=$?"
done
