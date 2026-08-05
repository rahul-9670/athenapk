#!/bin/bash
# Post-A2/A3/N1 gate chain. Runs unattended; every step's verdict is printed.
#
# ORDER MATTERS. Step 1 runs the 3D CT reference legs with the CURRENT build_cpu binary, which
# has N13/N14/N6 but NOT A2/A3/N1. That is the only chance to capture a genuine "before" for the
# N1 bit-identity claim -- once the rebuild happens the old binary is gone. A structural argument
# ("every edit is inside an ndim==2 branch") is not a measurement.
set -u
R=/beegfs/u/bbg6470/athenapk
CT=$R/runs/ct_tests/n1_2d
AF=$R/runs/audit_fix_regress
source ~/athenapk_env.sh >/dev/null 2>&1

# STEP 1 (capture the pre-N1 3D reference) is run SEPARATELY and must already have produced
# HASHES3D_pre. It is not re-run here on purpose: it has to execute against the OLD binary, and
# the rebuild below destroys that binary. On 2026-08-05 an earlier version of this script ran
# step 1 inline, four of six legs died on deck errors, and the script cheerfully proceeded to
# rebuild -- which would have made the bit-identity claim untestable. Hence the hard gate.
echo "######## STEP 0: pre-N1 3D reference must already exist ########"
if [ ! -s "$CT/HASHES3D_pre" ]; then
  echo "MISSING or empty $CT/HASHES3D_pre -- run"
  echo "  ONLY3D=1 OUTTAG=pre $CT/run_n1_2d.sh"
  echo "with the PRE-A2/A3/N1 binary FIRST. Refusing to rebuild and lose the reference."
  exit 1
fi
NPRE=$(grep -c '_3d ' "$CT/HASHES3D_pre")
echo "  HASHES3D_pre has $NPRE entries"
if [ "$NPRE" -lt 12 ]; then
  echo "  EXPECTED 12 (6 legs x 2 outputs). A short file means legs died -- fix them and re-run"
  echo "  step 1 before rebuilding."
  exit 1
fi

echo
echo "######## STEP 2: rebuild build_cpu with A2/A3/N1 ########"
git -C $R status --short
make -C $R/build_cpu athenaPK -j8 2>&1 | tail -5
RC=${PIPESTATUS[0]}
if [ "$RC" != "0" ]; then echo "BUILD FAILED (rc=$RC) -- stopping"; exit 1; fi
md5sum $R/build_cpu/bin/athenaPK
git -C $R status --short

echo
echo "######## STEP 3: round6 regression (A2/A3 must be inert; all decks are 3D/GLM) ########"
"$AF/regress.sh" round6_a2a3n1 A_multipole B_swindle C_idealeos D_sinkhydro E_sinkmhd \
    H_dust_ideal I_dust_hydrogen 2>&1 | tail -20

echo
echo "######## STEP 4: 3D CT legs on the POST binary -- must be BIT-IDENTICAL to step 1 ########"
ONLY3D=1 OUTTAG=post "$CT/run_n1_2d.sh" || echo "  (step 4 reported a non-zero exit; see above)"
echo "--- diff HASHES3D_pre vs HASHES3D_post ---"
if diff "$CT/HASHES3D_pre" "$CT/HASHES3D_post"; then
  echo "  => N1 3D BIT-IDENTITY: PASS"
else
  echo "  => N1 3D BIT-IDENTITY: FAIL (3D changed; N1 was supposed to touch 2D only)"
fi

echo
echo "######## STEP 5: N1 2D validation ladder (2D vs 3D + falsifier) ########"
"$CT/run_n1_2d.sh"

echo
echo "######## STEP 6: G_ct2d_guard must now RUN instead of aborting ########"
G=$AF/results/round6_ct2d; rm -rf $G; mkdir -p $G
export OMP_NUM_THREADS=1 OMPI_MCA_mtl='^psm2' OMPI_MCA_btl='tcp,self' OMPI_MCA_pml='ob1' \
       OMPI_MCA_io=romio341
( cd $G && $R/build_cpu/bin/athenaPK -i $AF/G_ct2d_guard.in ) > $G/run.log 2>&1
echo "  exit=$? cycles=$(grep -c '^cycle=' $G/run.log)"
grep -E "CT 2D out-of-plane check" $G/run.log | sed 's/^/  /'
grep -iE "PARTHENON ERROR|what\(\):" $G/run.log | head -3 | sed 's/^/  ERR: /'

echo
echo "######## DONE ########"
