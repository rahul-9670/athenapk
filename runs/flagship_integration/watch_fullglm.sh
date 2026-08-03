#!/bin/bash
cd /beegfs/u/bbg6470/athenapk/runs/flagship_integration
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
for i in $(seq 1 720); do   # up to ~6h, checking every 30s
  # blowup?
  nan=$(awk 'NR>2 && ($6=="-nan"){print $1; exit}' fc128glm/parthenon.out0.hst 2>/dev/null)
  [ -n "$nan" ] && { echo "FULL-GLM: NaN at t=$nan (GLM also blows => NOT just CT)"; exit 0; }
  # first core?
  [ -f STOP_128glm ] && { echo "FULL-GLM: STOP -> $(cat STOP_128glm)"; exit 0; }
  # progress + ME/E health from newest phdf
  NEW=$(ls -t fc128glm/parthenon.out1.*.phdf 2>/dev/null | head -1)
  if [ -n "$NEW" ]; then
    R=$($PY - "$NEW" <<'PYEOF'
import h5py,numpy as np,sys
try:
 with h5py.File(sys.argv[1],"r") as h:
  t=float(h["Info"].attrs["Time"]);p=h["prim"]
  rho=np.array(p[:,0]);P=np.array(p[:,4]);b=np.array(p[:,5])**2+np.array(p[:,6])**2+np.array(p[:,7])**2
  v=np.array(p[:,1])**2+np.array(p[:,2])**2+np.array(p[:,3])**2
  ME=0.5*b;KE=0.5*rho*v;E=P/0.666667+KE+ME;mef=(ME/np.maximum(E,1e-30))
  print("t=%.4f rho_max=%.2ecgs maxME/E=%.3f"%(t,rho.max()*5.467e-19,mef.max()))
except Exception as e: print("pending",e)
PYEOF
)
    echo "[$(date +%H:%M)] FULL-GLM $R"
  fi
  # if this member ended and no successor and not first core, note it
  q=$(squeue -u bbg6470 -h -o "%j" 2>/dev/null | grep -c fc128glm)
  [ "$q" = "0" ] && [ -f fc128glm/parthenon.out0.hst ] && { echo "FULL-GLM: no job in queue (chain ended); last $(tail -1 fc128glm/parthenon.out0.hst 2>/dev/null|awk '{print \"t=\"$1\" cyc=\"$3}')"; exit 0; }
  sleep 30
done
echo "FULL-GLM watcher timeout (still running)"
