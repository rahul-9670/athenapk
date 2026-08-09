#!/usr/bin/env python
"""Gate 3 for the WP-08 density-split restart: are the new columns REAL?

Three checks, in order of what would embarrass us most if skipped:

1. PRESENT -- the six WP-8 columns exist at all. If mag_diag_rho_split did not survive the
   restart's stored deck, they are simply absent and everything below is moot.
2. NON-DEGENERATE -- neither bin is identically zero. A split where one bin holds everything is
   the failure mode BOTH previous WP-8 splits died of (the Jsq density split put 97-99 % in one
   bin; the current-sheet split at thresh=0.3 left the other bin EMPTY). A column of zeros would
   look like a clean measurement in a plot and mean nothing.
3. RECONSTRUCTS -- hi + lo == global to round-off. This is the falsifier for the columns
   themselves: the two bins partition the same cells the global sums, so if they do not add up,
   the kernel's masking is wrong and no convergence claim built on them is safe. Checked against
   the global mag-dissO / mag-dissA in the SAME row, so no cross-file or cross-epoch mismatch can
   creep in.

Reports the volume fractions too: mag-Vhi vs the box volume says how much of the domain the
"core" bin actually is, which is the number that tells us whether this bin is a physical region
or another point sample.
"""
import sys
import re
import numpy as np

path = sys.argv[1] if len(sys.argv) > 1 else "parthenon.out0.hst"
# Box volume in CODE units, for f_eff. The njeans ladder runs a cube [-8, 8]^3 (fhc_ladder.in
# <parthenon/mesh>), so V_box = 16^3 = 4096. AMR does not change it -- refinement subdivides the
# same domain. Override as argv[2] if this is ever pointed at a different deck.
VBOX = float(sys.argv[2]) if len(sys.argv) > 2 else 16.0 ** 3

try:
    cols = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") and "[" in line:
                cols = re.findall(r"\[\d+\]=(\S+)", line)
                break
    d = np.loadtxt(path, comments="#")
except Exception as e:
    print(f"    CANNOT READ {path}: {e}")
    sys.exit(1)
if d.ndim == 1:
    d = d[None, :]

# NB: the concentration probes are mag-dissOsq / mag-dissAsq, NOT mag-dissOV / mag-dissAV.
# The first version of this checker demanded the -OV/-AV names and reported a FALSE failure on a
# run that was completely fine (job 2495356). It got them from the doc block at the top of
# mag_diag.hpp, which still advertised the old names -- but the enum right below it records that
# a literal "volume where q>0" was implemented first and is USELESS (eta and J are nonzero almost
# everywhere, so it just returns the box volume: measured 1.40608e5 on an L=52 deck, i.e. exactly
# 52^3) and was replaced by int q^2 dV, from which analysis forms
#     f_eff = (int q dV)^2 / (V_box * int q^2 dV).
# The stale header has been corrected. Lesson: verify column names against the REGISTRATION in
# hydro.cpp, not against prose.
need = ["mag-dissO-hi", "mag-dissO-lo", "mag-dissA-hi", "mag-dissA-lo",
        "mag-Vhi", "mag-dissOsq", "mag-dissAsq",
        # WP-8 round 3
        "mag-dissO-sheet", "mag-dissO-smooth", "mag-dissA-sheet", "mag-dissA-smooth",
        "mag-dissOhisq", "mag-dissOlosq", "mag-dissAhisq", "mag-dissAlosq", "mag-Vlo"]
missing = [c for c in need if c not in cols]
if missing:
    print(f"    *** SPLIT COLUMNS MISSING: {missing}")
    print("    *** Check (a) that the binary CONTAINS them -- grep -qa the key straight at the")
    print("    ***       binary, never `strings | grep -q` (SIGPIPE + pipefail => false ABSENT);")
    print("    ***   and (b) that hydro/mag_diag_rho_split reached the run.")
    print("    *** Do NOT launch nj8/nj16 until this passes.")
    sys.exit(2)

g = lambda n: d[:, cols.index(n)]
print(f"    rows={d.shape[0]}  t={d[0,0]:.6f} .. {d[-1,0]:.6f}")

ok = True
for base in ("dissO", "dissA"):
    tot, hi, lo = g(f"mag-{base}"), g(f"mag-{base}-hi"), g(f"mag-{base}-lo")
    # last row: the run has settled past any restart transient
    T, H, L = tot[-1], hi[-1], lo[-1]
    rec = H + L
    rel = abs(rec - T) / abs(T) if T else (0.0 if rec == 0 else float("inf"))
    frac = H / T if T else float("nan")
    zero = "  ** A BIN IS IDENTICALLY ZERO" if (np.all(hi == 0) or np.all(lo == 0)) else ""
    # Tolerance 1e-5, NOT 1e-10. The .hst is a text file printed to SIX significant figures, so
    # hi + lo can only agree with the global to ~1e-6/1e-7 relative no matter how exact the
    # kernel is. The first version used 1e-10 and duly reported "DOES NOT RECONSTRUCT" at
    # rel = 2.8e-07 on a perfectly correct run -- reading the file's own precision as a physics
    # defect. Anything above 1e-5 is a genuine masking error; below it is the print format.
    bad = "  ** DOES NOT RECONSTRUCT" if rel > 1e-5 else ""
    if zero or bad:
        ok = False
    print(f"    mag-{base:5s} global={T:12.5e}  hi={H:12.5e} ({frac*100:6.2f}%)  "
          f"lo={L:12.5e}  |hi+lo-glob|/glob={rel:.2e}{bad}{zero}")

vhi = g("mag-Vhi")[-1]
print(f"    mag-Vhi={vhi:.5e}   mag-dissOsq={g('mag-dissOsq')[-1]:.5e}   "
      f"mag-dissAsq={g('mag-dissAsq')[-1]:.5e}")
# f_eff = (int q dV)^2 / (V_box * int q^2 dV): the fraction of the box that would carry the whole
# integral if q were uniform. ~1 = genuinely volume-filling; ~1e-7 = a point sample of a few
# cells, which cannot converge under refinement. This is the number WP-8 exists to expose.
for base in ("dissO", "dissA"):
    I, Isq = g(f"mag-{base}")[-1], g(f"mag-{base}sq")[-1]
    if Isq > 0:
        print(f"    f_eff({base}) = {I*I/(VBOX*Isq):.3e}"
              + ("   ** POINT SAMPLE" if I*I/(VBOX*Isq) < 1e-3 else ""))

# WP-8 round 3: the sheet split must ALSO reconstruct the global, and per-bin f_eff is now
# formable because each bin has its own sq column and carrying volume.
for base in ("dissO", "dissA"):
    T = g(f"mag-{base}")[-1]
    sh, sm = g(f"mag-{base}-sheet")[-1], g(f"mag-{base}-smooth")[-1]
    rel = abs(sh + sm - T) / abs(T) if T else 0.0
    bad = "  ** DOES NOT RECONSTRUCT" if rel > 1e-5 else ""
    if bad:
        ok = False
    print(f"    mag-{base:5s} sheet={sh:12.5e} ({sh/T*100 if T else float('nan'):6.2f}%)  "
          f"smooth={sm:12.5e}  |sheet+smooth-glob|/glob={rel:.2e}{bad}")
vlo = g("mag-Vlo")[-1]
for base in ("dissO", "dissA"):
    for tag, V in (("hi", vhi), ("lo", vlo)):
        I = g(f"mag-{base}-{tag}")[-1]; Isq = g(f"mag-{base}{tag}sq")[-1]
        if Isq > 0 and V > 0:
            fe = I * I / (V * Isq)
            print(f"    f_eff({base}-{tag}) = {fe:.3e}"
                  + ("   ** POINT SAMPLE within its own bin" if fe < 1e-3 else
                     "   (resolved within its bin)"))

print("    GATE 3: " + ("PASS -- columns present, both bins populated, hi+lo == global"
                        if ok else "FAIL -- see markers above"))
sys.exit(0 if ok else 3)
