#!/bin/bash
cd /beegfs/u/bbg6470/athenapk/runs/flagship_integration
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python; JOB=2416619
for i in $(seq 1 240); do
  st=$(squeue -j $JOB -h -o %T 2>/dev/null)
  reached=$(awk 'NR>2 && $1>=1.102905{print $1; exit}' probe_glmadon/parthenon.out0.hst 2>/dev/null)
  if [ -z "$st" ] || [ -n "$reached" ]; then
    NEW=$(ls -t probe_glmadon/parthenon.out1.*.phdf probe_glmadon/parthenon.out2.*.rhdf 2>/dev/null | head -1)
    R=$($PY - "$NEW" 2>/dev/null <<'PYEOF'
import h5py,numpy as np,sys
f=sys.argv[1]
with h5py.File(f,"r") as h:
  t=float(h["Info"].attrs["Time"])
  k="prim" if "prim" in h else "cons"; c=np.array(h[k])
  if k=="prim": rho=c[:,0];P=c[:,4];b=c[:,5]**2+c[:,6]**2+c[:,7]**2;v=c[:,1]**2+c[:,2]**2+c[:,3]**2; KE=0.5*rho*v; eint=P/0.666667
  else: rho=c[:,0];b=c[:,5]**2+c[:,6]**2+c[:,7]**2;m=c[:,1]**2+c[:,2]**2+c[:,3]**2;KE=0.5*m/np.maximum(rho,1e-30);E4=c[:,4];eint=E4-KE-0.5*b
  ME=0.5*b;E=eint+KE+ME;mef=ME/np.maximum(E,1e-30)
  print("t=%.5f maxME/E=%.3f"%(t,mef.max()))
PYEOF
)
    echo "VERDICT GLM-AD-ON (clean same-start): $R  (compare CT-AD-ON=0.689, GLM=0.107 @ t=1.10291). AD-off~0.689=>AD-under-CT BROKEN; AD-off>>0.689=>ideal-CT amplifies"
    exit 0
  fi
  sleep 20
done
echo "ctadoff watcher timeout"
