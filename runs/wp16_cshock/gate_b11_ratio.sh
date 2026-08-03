#!/bin/bash
#SBATCH --job-name=b11gate
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#
# FULL GATE for B11 — `diffusion/hall_ohmic_floor_ratio`, the per-cell Hall Ohmic stabilizer.
#
# WHAT CHANGED. `hall_ohmic_floor_code` is an ABSOLUTE stabilizer applied regardless of the
# local eta_H, and WP-16 part 3 measured that 3D Hall AMPLIFIES rather than damps once
# eta_floor/|eta_H| falls below ~0.11. The new key makes the applied stabilizer
#     eta_floor_cell = max(hall_ohmic_floor_code, hall_ohmic_floor_ratio * |eta_H_cell|)
# so it tracks the term it stabilizes. Default 0.0 = disabled = the old constant.
#
# EVERY CONSUMER OF THE FLOOR WAS UPDATED, not just the EMF -- an EMF-only change would have
# produced a stable-but-wrongly-stepped run:
#   hall.cpp   HallDiffFluxIsoFixed x3 (X1/X2/X3 EMF, GLM path)
#   hall.cpp   EstimateHallTimestep    (the parabolic dx^2/eta_floor constraint)
#   ct.cpp     CT_AddHallEMF x3        (E1/E2/E3 edge EMFs, CT path)
#   diffusion.cpp  fused unsplit dt estimator  (eta_O_tot = eta_O + floor)
#   diffusion.cpp  fused mixed dt estimator    (eta_floor_par / eta_floor_strict)
# and in the mixed rkl2 mode eta_H is now evaluated on the FLOOR-ONLY kernel call too
# (`need_eta_h`), without which the fix would have been silently inert in exactly the
# integrator production uses.
#
# FOUR GATES:
#   A  OFF-state: ratio unset must be BYTE-IDENTICAL to the pre-B11 binary 2ddea223 on both
#      the GLM and CT Hall decks. Controls are that binary's own recorded history files.
#   B  ON-state: ratio=0.2 must FIX the 3D instability at N=128 (which amplifies 2.65e5) and
#      N=256 (which crashes).
#   C  PRODUCTION EQUIVALENCE: production runs floor=0.05 with max|eta_H| = 4.3e-2 measured
#      over 78.7M cells of prod_v9, so 0.2*|eta_H| <= 8.6e-3 < 0.05 in EVERY cell and the
#      max() must select the absolute floor everywhere => enabling the ratio must be
#      byte-identical on the production-physics deck. This is the "nothing else breaks" test.
#   D  Same as C but through the MIXED RKL2 path (production's actual integrator), which is a
#      different code path from the unsplit one C exercises.
set -o pipefail
source ~/athenapk_env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OMP_PROC_BIND=spread OMP_PLACES=threads
export OMPI_MCA_pml=ob1 OMPI_MCA_io=romio341 FI_PROVIDER=tcp PMIX_MCA_gds=hash

B=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK
H=/beegfs/u/bbg6470/athenapk/runs/wp16_cshock
I=/beegfs/u/bbg6470/athenapk/inputs
FHC=/beegfs/u/bbg6470/athenapk/runs/b2b4_gate/fhc.in
echo "B11 gate, job $SLURM_JOB_ID $(date)"; md5sum $B

run () { # dir  deck  extra-args...
  local G=$H/$1; shift; local D=$1; shift
  rm -rf $G; mkdir -p $G
  ( cd $G && $B -i $D "$@" > run.log 2>&1 )
  echo "  exit=$?"
}
cmp_hst () { # label  fileA  fileB
  if [ ! -f "$2" ] || [ ! -f "$3" ]; then echo "  $1: MISSING ($2 | $3)"; return; fi
  if cmp -s "$2" "$3"; then echo "  $1: PASS byte-identical"; else
    echo "  $1: **FAIL** differs"; diff "$2" "$3" | head -6; fi
}

echo "=== GATE A — OFF-state (ratio unset) vs the pre-B11 binary 2ddea223 ==="
run a_glm $I/hall_whistler_glm.in
cmp_hst "A1 GLM unsplit Hall" $H/gate_b11_new/parthenon.out1.hst $H/a_glm/parthenon.out1.hst
run a_ct  $I/hall_whistler_ct.in
cmp_hst "A2 CT Hall        " $H/stock_hall_whistler_ct/parthenon.out1.hst $H/a_ct/parthenon.out1.hst

echo "=== GATE B — ON-state: does ratio=0.2 fix the 3D instability? ==="
for N in 128 256; do
  E=$(/beegfs/u/bbg6470/venvs/analysis_env/bin/python -c "print(2.0/$N)")
  for R in 0.0 0.2; do
    run b_n${N}_r${R} $I/hall_whistler_glm.in \
      parthenon/mesh/nx1=$N parthenon/meshblock/nx1=$N \
      parthenon/mesh/x2min=-$E parthenon/mesh/x2max=$E \
      parthenon/mesh/x3min=-$E parthenon/mesh/x3max=$E \
      diffusion/hall_ohmic_floor_ratio=$R
    echo "  N=$N ratio=$R: $(grep -oE 'rel. error \(omega\).*' $H/b_n${N}_r${R}/run.log | tail -1)"
    grep -oE "amplitude in/out.*" $H/b_n${N}_r${R}/run.log | sed 's/^/      /'
  done
done

echo "=== GATE C — production equivalence, UNSPLIT (floor 0.05 dominates 0.2*|eta_H|) ==="
run c_off $FHC
run c_on  $FHC diffusion/hall_ohmic_floor_ratio=0.2
cmp_hst "C production unsplit" $H/c_off/parthenon.out0.hst $H/c_on/parthenon.out0.hst

echo "=== GATE D — same, through the MIXED RKL2 path (production's integrator) ==="
RKL="diffusion/integrator=rkl2 diffusion/hall_floor_integrator=rkl2 diffusion/rkl2_max_dt_ratio=1000 diffusion/rkl2_freeze_eta=true"
run d_off $FHC $RKL
run d_on  $FHC $RKL diffusion/hall_ohmic_floor_ratio=0.2
cmp_hst "D production rkl2  " $H/d_off/parthenon.out0.hst $H/d_on/parthenon.out0.hst

echo "=== banner check (the new startup notice) ==="
grep -h "Hall Ohmic stabilizer\|NOTE \[Hall\]\|WARNING \[Hall\]" $H/c_on/run.log $H/c_off/run.log | sed 's/^/  /' | head -4
echo "done $(date)"
