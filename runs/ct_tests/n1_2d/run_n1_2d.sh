#!/bin/bash
# AUDIT N1 VALIDATION -- 2D CT must now evolve B_z, and must agree with the validated 3D path.
#
# THE TEST. For each physics configuration, run the SAME z-invariant whistler eigenmode twice:
# once on a 2D mesh (nx3 = 1) and once on a 3D mesh (nx3 = 4). The initial condition has no
# z-dependence and x3 is periodic, so the exact solution stays z-invariant forever and the two
# legs must agree cell for cell up to floating-point reassociation. The 3D leg is the reference
# because that path is already validated (whistler dispersion 0.4%, DEV_LOG 2026-06).
#
# WHY THIS IS THE RIGHT TEST rather than an analytic one: it exercises the ideal E1/E2 assembly,
# the F3 curl, and the non-ideal Ohmic/ambipolar/Hall E1/E2 kernels simultaneously, against an
# answer nobody in this session derived. A sign error in any 2D branch shows up immediately.
#
# dt_ceil FORCES BOTH LEGS ONTO THE SAME TIMESTEP. The hydro CFL sums 1/dx over the ACTIVE
# directions, so an unconstrained 3D leg would take a smaller dt than the 2D leg and the two
# would separate at truncation-error level -- a real difference that says nothing about CT.
# Capping both below either CFL limit removes that confound; then any disagreement is the scheme.
#
# THE FALSIFIER. Leg `frozen_check` re-runs the ideal configuration and reports how far B_z
# actually moves. If the 2D fix were inert (B_z still frozen), max|B_z - B_z(t=0)| would be 0
# and the whole comparison would be passing for the wrong reason -- two frozen fields agree
# trivially. That number must be NON-ZERO and O(amp).
set -u
cd "$(dirname "$0")"
ROOT=$PWD
source ~/athenapk_env.sh >/dev/null 2>&1
BIN=${BIN:-/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK}
PY=/beegfs/u/bbg6470/venvs/analysis_env/bin/python
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false
export OMPI_MCA_mtl='^psm2' OMPI_MCA_btl='tcp,self' OMPI_MCA_pml='ob1' OMPI_MCA_io=romio341

# dt below either leg's CFL limit (dx = 7.8125e-3, and Hall shrinks it much further).
DTC=${DTC:-2.0e-5}

run_leg () {         # name  nx3  extra-args...
  local name=$1 nx3=$2; shift 2
  local wd=out/$name
  rm -rf "$wd"; mkdir -p "$wd"
  local x3half=0.015625
  local mbnx3=$nx3
  ( cd "$wd" && "$BIN" -i "$ROOT/base.in" \
      parthenon/mesh/nx3=$nx3 parthenon/meshblock/nx3=$mbnx3 \
      parthenon/mesh/x3min=-$x3half parthenon/mesh/x3max=$x3half \
      parthenon/time/dt_ceil=$DTC \
      "$@" ) > "$wd/run.log" 2>&1
  local rc=$?
  # Report the LAST CYCLE NUMBER, not the count of printed lines: ncycle_out = 1000 means the
  # log carries only cycle=0 and cycle=1000, so a line count reads "2" for a 1000-step run and
  # invites the conclusion that the leg barely ran. (It did: dt is pinned at dt_ceil, so this is
  # tlim/dt_ceil = 1000 steps.)
  local lastcyc=$(grep -oE '^cycle=[0-9]+' "$wd/run.log" | tail -1 | cut -d= -f2)
  local nprint=$(grep -c '^cycle=' "$wd/run.log")
  # POSITIVE CHECK: exit 0 is not evidence -- this cluster has produced exit-0 runs that died
  # in seconds. A leg counts only if it cycled and wrote its output.
  if [ $rc -ne 0 ] || [ "$nprint" -eq 0 ]; then
    echo "  VOID $name (exit=$rc, no cycle lines):"
    grep -A4 -iE "PARTHENON ERROR|Kokkos ERROR|what\(\):|FATAL" "$wd/run.log" | head -8 |
      sed 's/^/     /'
    return 1
  fi
  echo "  ok   $name  last cycle=${lastcyc:-?}"
  return 0
}

# "name;extra args" -- the args go to BOTH legs of that pair. `arithmetic` is the
# Balsara-Spicer cell-centred average (hydro.cpp:657); `gs05` is the upwind corner EMF the
# flagship uses. Both have their own 2D branch and both are covered.
# The `ideal_*` legs still carry hall_coeff_code in base.in -- see the note there; it only shapes
# the initial eigenmode. With `diffusion/hall` unset the evolution is pure ideal MHD.
# The non-ideal legs each need their `*_coeff = fixed` selector as well as the coefficient value;
# hydro.cpp:1246/1340 hard-fail if the term is enabled without one.
declare -a PAIRS=(
  "ideal_gs05;hydro/ct_emf=gs05"
  "ideal_bs;hydro/ct_emf=arithmetic"
  "ohmic;hydro/ct_emf=gs05 diffusion/resistivity=ohmic diffusion/resistivity_coeff=fixed diffusion/ohm_diff_coeff_code=0.01"
  "ambipolar;hydro/ct_emf=gs05 diffusion/ambipolar=ambipolar diffusion/ambipolar_coeff=fixed diffusion/ambipolar_coeff_code=0.05"
  "hall;hydro/ct_emf=gs05 diffusion/hall=hall diffusion/hall_coeff=fixed diffusion/hall_coeff_code=0.5 diffusion/hall_ohmic_floor_code=0.05"
  # Multi-block: 4 blocks in x1 forces the F3/E1/E2 ghost exchange to actually run in 2D. A
  # single-block test cannot see a broken face-field boundary.
  "ideal_mb;hydro/ct_emf=gs05 parthenon/meshblock/nx1=32"
)

echo "binary : $BIN"; md5sum "$BIN"
(cd /beegfs/u/bbg6470/athenapk && git rev-parse HEAD)
echo "dt_ceil: $DTC"
echo

# ONLY3D=1 runs just the 3D legs and hashes them. That is the N1 REGRESSION gate, distinct from
# the validation gate: every N1 edit is either inside an `ndim == 2` branch or an offset that is
# -1 in 3D, so 3D must be BIT-IDENTICAL across the change. Run it once with the pre-N1 binary,
# once after, and diff the two HASHES files. A structural argument is not a measurement.
ONLY3D=${ONLY3D:-0}
OUTTAG=${OUTTAG:-}

FAIL=0
for spec in "${PAIRS[@]}"; do
  name=${spec%%;*}
  read -r -a xargs <<< "${spec#*;}"
  echo "=== $name ==="
  run_leg "${name}_3d" 4 "${xargs[@]}" || { FAIL=1; continue; }
  [ "$ONLY3D" = "1" ] && continue
  run_leg "${name}_2d" 1 "${xargs[@]}" || { FAIL=1; continue; }
  $PY compare_2d3d.py "out/${name}_2d" "out/${name}_3d" "$name" || FAIL=1
done

if [ "$ONLY3D" = "1" ]; then
  echo
  echo "=== 3D CONTENT HASHES${OUTTAG:+ ($OUTTAG)} ==="
  $PY hash3d.py out | tee "HASHES3D${OUTTAG:+_$OUTTAG}"
  exit $FAIL
fi

echo
echo "=== FALSIFIER: is B_z actually moving in 2D? ==="
$PY frozen_check.py out/ideal_gs05_2d || FAIL=1

echo
[ $FAIL -eq 0 ] && echo "N1 2D CT VALIDATION: PASS" || echo "N1 2D CT VALIDATION: FAIL"
exit $FAIL
