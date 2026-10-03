#!/bin/bash
# Build the high resolution EOS table used for protostellar collapse runs:
#   400 x 1399 (rho, esp) nodes and 400 x 1370 (rho, T) nodes, 17.8 MB,
#   rho = 1e-20 .. 1 g/cm^3, T = 0.046 K .. 3e5 K, esp down to ~3.4e6 erg/g.
# It is the 400 x 1000 x 920 table with 399 esp nodes and 450 T nodes added below the
# axes (at the same spacing), so cold gas is still inside the table.
# Takes about 7 minutes on one core and needs only numpy.
#
# usage: build_hires_v5_cold.sh [output file]
#        (default: eos_table_hires_v5_cold.bin in the current directory)
set -euo pipefail
PYTHON=${PYTHON:-python3}
OUT=${1:-eos_table_hires_v5_cold.bin}
HERE=$(cd "$(dirname "$0")" && pwd)
"$PYTHON" "$HERE/gen_eos_table.py" --out "$OUT" --nr 400 --ne 1000 --nT 920 \
  --esp-pad 399 --T-pad 450
# Reference build (numpy 2.x, x86_64): md5 f712d6a806cb82889214fec95a6c2681
md5sum "$OUT"
