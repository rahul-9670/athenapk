#!/bin/bash
cd /beegfs/u/bbg6470/athenapk/runs/flagship_integration
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
BUILD=2416947; BIN=/beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK; before=d1d809256341f720a988389554522272
# 1) wait for build
for i in $(seq 1 240); do [ -z "$(squeue -j $BUILD -h -o %T 2>/dev/null)" ] && break; sleep 15; done
md=$(md5sum "$BIN" 2>/dev/null|awk '{print $1}')
[ "$md" = "$before" ] || [ -z "$md" ] && { echo "BUILD FAILED/unchanged md5=$md"; exit 0; }
err=$(ls -t /beegfs/u/bbg6470/athenapk/runs/*build_gpu*.out 2>/dev/null|head -1)
echo "BUILD OK new md5=$md"; grep -ciE "error:|Error 1" $err 2>/dev/null | sed 's/^/  build errors: /'
# 2) launch confirm
rm -rf confirm_ctfix; JID=$(sbatch --parsable confirm_ctfix.sh); echo "confirm job=$JID"
# 3) track ME/E(t) trajectory (ref: OLD CT 0.393@1.100 -> 0.689@1.10291; GLM ~0.428@same)
last=""
for i in $(seq 1 600); do
  # blowup?
  nan=$(awk 'NR>2 && ($6=="-nan"){print $1; exit}' confirm_ctfix/parthenon.out0.hst 2>/dev/null)
  [ -n "$nan" ] && { echo "CONFIRM: NaN at t=$nan"; }
  NEW=$(ls -t confirm_ctfix/parthenon.out1.*.phdf 2>/dev/null|head -1)
  if [ -n "$NEW" ] && [ "$NEW" != "$last" ]; then
    last=$NEW
    $PY - "$NEW" <<PYEOF
import h5py,numpy as np,sys
with h5py.File("$NEW","r") as h:
    t=float(h["Info"].attrs["Time"]);p=h["prim"]
    rho=np.array(p[:,0]);P=np.array(p[:,4]);b=np.array(p[:,5])**2+np.array(p[:,6])**2+np.array(p[:,7])**2
    v=np.array(p[:,1])**2+np.array(p[:,2])**2+np.array(p[:,3])**2
    ME=0.5*b;KE=0.5*rho*v;E=P/0.666667+KE+ME;mef=ME/np.maximum(E,1e-30)
    print("  CONFIRM t=%.5f maxME/E=%.3f  [old-CT ref: 0.393@1.100->0.689@1.10291; GLM 0.428]"%(t,mef.max()))
PYEOF
  fi
  st=$(squeue -j $JID -h -o %T 2>/dev/null)
  # done if reached tlim or job ended
  reached=$(awk 'NR>2 && $1>=1.102905{print 1;exit}' confirm_ctfix/parthenon.out0.hst 2>/dev/null)
  { [ -z "$st" ] || [ -n "$reached" ] || [ -n "$nan" ]; } && { echo "CONFIRM DONE: job_state=${st:-ended} reached_tlim=${reached:-no} nan=${nan:-no}"; break; }
  sleep 30
done
