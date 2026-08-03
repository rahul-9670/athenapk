#!/bin/bash
# Phase 2 Constrained Transport validation driver (increments 1 + 2).
#
# Runs the field-loop and Orszag-Tang tests with divergence_control = ct and = glm,
# plus a STATIC-AMR field loop (increment 2) whose coarse-fine boundary bisects the
# loop, and prints the divergence diagnostics. Single-rank CPU; front-end safe (needs
# the romio341 MPI-IO workaround + psm2 disabled, see the build-run-env note).
#
# Usage:  bash run_ct_tests.sh          # runs on the front-end (serial)
# For multi-rank / larger grids submit through SLURM std partition instead.
set -e
source ~/athenapk_env.sh >/dev/null 2>&1

BIN=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
HERE=/beegfs/u/bbg6470/athenapk/runs/ct_tests
export OMP_NUM_THREADS=1
export OMPI_MCA_mtl=^psm2 OMPI_MCA_btl=tcp,self OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341

run() {  # <subdir> <input> <nlim>
  local dir="$HERE/$1"; rm -rf "$dir"; mkdir -p "$dir"; cd "$dir"
  # dt=1000 suppresses phdf dumps during the divB sweep (every input now carries an
  # output0 hdf5 block so this override is always valid).
  "$BIN" -i "$HERE/$2" parthenon/time/nlim="$3" parthenon/output0/dt=1000 >run.log 2>&1
}

# maxcol <hst> <column-name>  -> max over the run of the named history column, found by
# parsing the "# [i]=name" header so the index never has to be hard-coded.
maxcol() {
  awk -v want="$2" '
    /^# \[/ { for (i=1;i<=NF;i++) if ($i ~ ("\\]=" want "$")) { s=index($i,"["); e=index($i,"]");
              ci=substr($i,s+1,e-s-1)+0 } next }
    /^#/    { next }
    { if (ci>0 && ($ci+0)>m) m=$ci+0 }
    END { printf "%.4e", m }' "$1"
}

# NOTE on the metric: ct_maxAbsDivB = max|div B|_face*dx is the true CT invariant
# (div(curl E)=0 identically -> round-off for any single-valued edge EMF). ct_maxRelDivB
# additionally divides by |B_cc|, which blows up in near-zero-field regions (e.g. the
# field-loop ambient), so it can read ~1e-9 from pure round-off there -- PASS/FAIL uses
# the absolute metric.
echo "== field loop (single level) =="
run fl_ct  field_loop_ct.in  140
run fl_glm field_loop_glm.in 140
echo "  CT  max |divB|*dx (abs) = $(maxcol "$HERE"/fl_ct/*.hst ct_maxAbsDivB)   (expect ~1e-18, round-off)"
echo "  GLM max CC-divB (rel)   = $(maxcol "$HERE"/fl_glm/*.hst UserRelDivB)   (grows: GLM cleans advectively)"

echo "== field loop (STATIC AMR, C-F boundary bisects loop) -- increment 2 =="
run fl_ct_amr field_loop_ct_amr.in 140
NB=$(awk '/^#/{next}{print $4; exit}' "$HERE"/fl_ct_amr/*.hst)
CT_AMR=$(maxcol "$HERE"/fl_ct_amr/*.hst ct_maxAbsDivB)
echo "  nbtotal=$NB  CT max |divB|*dx across C-F = $CT_AMR"
awk -v v="$CT_AMR" 'BEGIN{ if (v+0 < 1e-12) print "  RESULT: PASS (div-free preserved across coarse-fine boundary)";
                            else print "  RESULT: FAIL (edge-EMF reflux broken -> divB spike at C-F)" }'

echo "== orszag-tang (single level) =="
run ot_ct  orszag_tang_ct.in  200
run ot_glm orszag_tang_glm.in 200
echo "  CT  max |divB|*dx (abs) = $(maxcol "$HERE"/ot_ct/*.hst ct_maxAbsDivB)   (expect ~1e-18, round-off)"
echo "  CT  max face-divB (rel) = $(maxcol "$HERE"/ot_ct/*.hst ct_maxRelDivB)   (B~O(1) so rel is clean too)"
