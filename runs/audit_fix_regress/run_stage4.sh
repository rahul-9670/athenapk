#!/bin/bash
# Stage-4 verification driver (A1 + A4/A6). Serial: the front-end is a 1-CPU cgroup.
#
#   gates           -- (default) re-print the verdicts from whatever is in results/
#   phase1          -- capture D/E baselines with the PRE-FIX binary (BEFORE rebuilding)
#   phase2          -- run every deck with the POST-FIX binary
#   phase_e_old     -- A4 A/B: E with only the sink energy line reverted
#   phase_a1discrim -- A1 A/B: F with only jeans.cpp reverted (needs $A1OLD)
#
# $OLD below is a session scratchpad path that no longer exists; phase1 is kept for the
# record. Regenerate a pre-fix binary the same way phase_a1discrim documents.
set -u
cd "$(dirname "$0")"
OLD=/beegfs/tmp/hummel2-front1/bbg6470/1783349135.889008.252851392/claude-7086470/-beegfs-u-bbg6470/11ff3eb5-c165-4c94-b9d9-04aa97ede850/scratchpad/athenaPK.baseline
NEW=/beegfs/u/bbg6470/athenapk/build_cpu/bin/athenaPK

case "${1:-gates}" in
  gates) ;;   # re-print the verdicts from whatever is already in results/
  phase1)
    # D only: the OLD binary predates the b0z option in poisson_test, so an E baseline taken
    # here would run with B = 0 and would compare two changes at once. E's A/B is done
    # separately (phase_e), against a build that reverts ONLY the one energy line.
    BIN=$OLD ./regress.sh baseline_sinkhydro D_sinkhydro
    ;;
  phase2)
    BIN=$NEW ./regress.sh stage4_ideal     C_idealeos
    BIN=$NEW ./regress.sh stage4_sinkhydro D_sinkhydro
    BIN=$NEW ./regress.sh stage4_sinkmhd   E_sinkmhd
    BIN=$NEW ./regress.sh stage4_eos       A_multipole B_swindle
    ;;
  phase_e_old)   # run E with a binary whose ONLY difference is the reverted energy line
    BIN=$NEW ./regress.sh oldenergy_sinkmhd E_sinkmhd
    ;;
  phase_a1discrim)
    # A1 POSITIVE test. A_multipole/B_swindle cannot discriminate the A1 fix: they sit at
    # numlevel = 2 with njeans = 8 and every block is already at the top level, so nbtotal is
    # pinned at 64 for the whole run and the refinement criterion never gets to act. Their
    # "identical" result is a NULL result, not evidence.
    #
    # F_a1discrim.in is built to discriminate. From the t=0 snapshot (scratchpad scan, same
    # bilinear/bisection reproduction as measure_A1.py, against src/eos/eos_table.bin):
    #     min over blocks of n_J(ideal gamma=1.4) = 4.537   (8 blocks below 4.7)
    #     min over blocks of n_J(EOS table)       = 4.925   (0 blocks below 4.7)
    # so njeans = 4.7 with numlevel = 3 splits the two criteria:
    #     OLD -> 8 blocks refine -> 64 - 8 + 8*8 = 120 blocks
    #     NEW -> 0 blocks refine ->                  64 blocks
    # $A1OLD must be a build with ONLY src/refinement/jeans.cpp reverted:
    #     git show <pre-fix-sha>:src/refinement/jeans.cpp > src/refinement/jeans.cpp
    #     make -C build_cpu athenaPK -j     (the fix landed in 6e3a55c, so <pre-fix-sha> = 6e3a55c~1)
    BIN=${A1OLD:?set A1OLD to a build with only jeans.cpp reverted} \
        ./regress.sh a1discrim_a1old F_a1discrim
    BIN=$NEW ./regress.sh a1discrim_a1new F_a1discrim
    ;;
esac

echo "=============================== GATES ==============================="
gate () { # name baseline candidate expectation
  local n=$1 b=$2 c=$3 e=$4
  if [ ! -s "$b" ] || [ ! -s "$c" ]; then echo "$n : SKIP (missing hashes)"; return; fi
  if diff -q "$b" "$c" >/dev/null; then r=identical; else r=different; fi
  if [ "$r" = "$e" ]; then echo "$n : PASS ($r, expected $e)";
  else echo "$n : FAIL ($r, expected $e)"; diff "$b" "$c" | head -12; fi
}
gate "A1 gate A  eos=adiabatic must be untouched " results/baseline_ideal/HASHES results/stage4_ideal/HASHES identical
gate "A4 gate A  sinks, pure hydro: ME=0 => no-op" results/baseline_sinkhydro/HASHES results/stage4_sinkhydro/HASHES identical
gate "A4 positive sinks, MHD: must CHANGE         " results/oldenergy_sinkmhd/HASHES results/stage4_sinkmhd/HASHES different
# NULL by construction -- see phase_a1discrim above. These two decks never exercise the
# refinement criterion (nbtotal is 64 for every cycle of every run), and on the t=0 state the
# A4/A6 sink paths are off, so identical is the CORRECT expectation, not a passing accident.
gate "A1/A4      eos=hydrogen decks (null test)  " results/baseline/HASHES     results/stage4_eos/HASHES   identical

# The A1 positive test: predicted block counts, not just "the answer changed".
a1_gate () {
  local o=results/a1discrim_a1old/parthenon.out0.hst
  local n=results/a1discrim_a1new/parthenon.out0.hst
  if [ ! -s "$o" ] || [ ! -s "$n" ]; then echo "A1 positive  discriminating mesh : SKIP (run phase_a1discrim)"; return; fi
  local nbo nbn
  nbo=$(awk 'NR==3{print $4}' "$o"); nbn=$(awk 'NR==3{print $4}' "$n")
  if [ "$nbo" = "120" ] && [ "$nbn" = "64" ]; then
    echo "A1 positive  discriminating mesh : PASS (old $nbo, new $nbn -- both as predicted)"
  else
    echo "A1 positive  discriminating mesh : FAIL (old $nbo expected 120, new $nbn expected 64)"
  fi
}
a1_gate
