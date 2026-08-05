#!/bin/bash
# Bit-identity regression harness for the 2026-08-04 audit fix batch.
#
#   ./regress.sh <label>            run every deck, write results/<label>/
#   ./regress.sh <label> A_multipole   run one deck only
#
# Gate: results/<label>/HASHES must match the baseline byte for byte.
# The .hst is plain text and deterministic; the phdf is hashed through h5py on
# DATASET CONTENT ONLY, because the file carries wall-clock metadata that legitimately
# differs between runs.
set -u
cd "$(dirname "$0")"
ROOT=$PWD
source ~/athenapk_env.sh >/dev/null 2>&1

BIN=${BIN:-/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK}
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
LABEL=${1:?usage: regress.sh <label> [deck ...]}
shift
DECKS=${*:-"A_multipole B_swindle"}

OUT=results/$LABEL
mkdir -p "$OUT"
echo "binary : $BIN" | tee "$OUT/PROVENANCE"
md5sum "$BIN" | tee -a "$OUT/PROVENANCE"
(cd /beegfs/u/bbg6470/athenapk && git rev-parse HEAD) | tee -a "$OUT/PROVENANCE"

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=false
# A bare (non-mpirun) launch does not inherit the MCA flags ~/athenapk_env.sh puts on
# mpirun, and OpenMPI's psm2 MTL aborts in MPI_Init on the front-end. Same effect, via env.
export OMPI_MCA_mtl='^psm2'
export OMPI_MCA_btl='tcp,self'  # /dev/shm is full on the front-end; no shared-memory BTL
export OMPI_MCA_pml='ob1'
export OMPI_MCA_io=romio341   # OMPIO cannot open files on BeeGFS; same workaround the submit scripts use

rc=0
for d in $DECKS; do
  wd="$OUT/$d"
  rm -rf "$wd"; mkdir -p "$wd"
  ( cd "$wd" && "$BIN" -i "$ROOT/$d.in" ) > "$wd/run.log" 2>&1
  if [ $? -ne 0 ]; then
    echo "FAIL  $d : non-zero exit"; tail -25 "$wd/run.log"; rc=1; continue
  fi
  if grep -qiE "nan|error|abort" "$wd/run.log"; then
    echo "WARN  $d : suspicious token in log"; grep -inE "nan|error|abort" "$wd/run.log" | head -5
  fi
  echo "ok    $d"
done

# --- content hashes -------------------------------------------------------
$PY - "$OUT" $DECKS <<'EOF' | tee "$OUT/HASHES"
import sys, os, hashlib, glob
import numpy as np, h5py
out = sys.argv[1]
for d in sys.argv[2:]:
    wd = os.path.join(out, d)
    h = hashlib.sha256()
    hst = sorted(glob.glob(os.path.join(wd, "*.hst")))
    for f in hst:
        # drop comment lines: the header is stable but cheap insurance
        with open(f, "rb") as fh:
            for line in fh:
                if not line.startswith(b"#"):
                    h.update(line)
    print("%-14s hst      %s" % (d, h.hexdigest()[:32]))
    for f in sorted(glob.glob(os.path.join(wd, "*.phdf"))):
        hh = hashlib.sha256()
        with h5py.File(f, "r") as g:
            def walk(name, obj):
                if isinstance(obj, h5py.Dataset):
                    a = np.asarray(obj[()])
                    hh.update(name.encode())
                    hh.update(a.tobytes() if a.dtype != object else repr(a).encode())
            g.visititems(walk)
        print("%-14s %-8s %s" % (d, os.path.basename(f).split(".")[-2], hh.hexdigest()[:32]))
EOF
exit $rc
