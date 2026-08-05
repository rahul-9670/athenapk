#!/usr/bin/env python
"""Find an njeans value that DISCRIMINATES the A1 fix on the regression deck.

nj_table >= nj_ideal (the table c_s is larger at T < theta_rot), so any njeans with
    nj_ideal_min  <  njeans  <=  nj_table_min
makes the OLD code tag `refine` and the NEW code tag `same` for that block.
Prints the per-block windows and the widest common choice.
"""
import sys
import numpy as np
import h5py

sys.path.insert(0, "/beegfs/u/bbg6470/athenapk/runs/audit_fix_regress")
from measure_A1 import read_eos, bilin, asq_from_rho_pres  # noqa: E402

SNAP = sys.argv[1]
EOS_BIN = sys.argv[2] if len(sys.argv) > 2 else \
    "/beegfs/u/bbg6470/athenapk/src/eos/eos_table.bin"
GAMMA = 1.4

E = read_eos(EOS_BIN)
with h5py.File(SNAP, "r") as g:
    prim = g["prim"]
    nb = prim.shape[0]
    levels = g["Levels"][:]
    xmin = float(g["Info"].attrs["RootGridDomain"][0])
    xmax = float(g["Info"].attrs["RootGridDomain"][1])
    nx = int(g["Info"].attrs["RootGridSize"][0])
    dx_root = (xmax - xmin) / nx
    print("snap %s\n  blocks %d  levels %d..%d  dx_root %.5g"
          % (SNAP, nb, levels.min(), levels.max(), dx_root))
    nj_i = np.empty(nb)
    nj_t = np.empty(nb)
    for b in range(nb):
        blk = prim[b]
        rho = blk[0].astype(np.float64)
        pre = blk[4].astype(np.float64)
        bsq = blk[5].astype(np.float64)**2 + blk[6].astype(np.float64)**2 + \
            blk[7].astype(np.float64)**2
        r, p, bb = rho.ravel(), pre.ravel(), bsq.ravel()
        cs_i = np.sqrt(GAMMA * p / r)
        a2, _ = asq_from_rho_pres(E, r, p)
        cs_t = np.sqrt(np.maximum(a2, 1e-300))
        va = np.sqrt(bb / r)
        fac = 2.0 * np.pi / (dx_root / 2.0 ** int(levels[b]))
        nj_i[b] = fac * np.min((cs_i + va) / np.sqrt(r))
        nj_t[b] = fac * np.min((cs_t + va) / np.sqrt(r))

lo, hi = nj_i.min(), nj_t.min()
print("  nj_ideal : min %.6g  max %.6g" % (nj_i.min(), nj_i.max()))
print("  nj_table : min %.6g  max %.6g" % (nj_t.min(), nj_t.max()))
print("  nj_table >= nj_ideal everywhere: %s" % bool(np.all(nj_t >= nj_i - 1e-12)))
print("  per-block ratio nj_table/nj_ideal: min %.4f med %.4f max %.4f"
      % ((nj_t/nj_i).min(), np.median(nj_t/nj_i), (nj_t/nj_i).max()))

# Widest window: any njeans in (nj_i[b], nj_t[b]] flips block b from refine -> same.
width = nj_t - nj_i
b = int(np.argmax(width))
print("\n  widest single-block window: block %d  (%.6g, %.6g]  width %.4g"
      % (b, nj_i[b], nj_t[b], width[b]))
cand = 0.5 * (nj_i[b] + nj_t[b])
n_old = int((nj_i < cand).sum())
n_new = int((nj_t < cand).sum())
print("  candidate njeans = %.6g  ->  OLD refines %d blocks, NEW refines %d blocks"
      % (cand, n_old, n_new))

# Also scan a grid for the choice that maximises the tag disagreement.
grid = np.unique(np.concatenate([nj_i, nj_t]))
best = None
for k in range(len(grid) - 1):
    c = 0.5 * (grid[k] + grid[k+1])
    d = int(((nj_i < c) != (nj_t < c)).sum())
    dd = int(((nj_i > 2.5*c) != (nj_t > 2.5*c)).sum())
    if best is None or d + dd > best[1]:
        best = (c, d + dd, d, dd)
print("  best njeans = %.6g  -> %d blocks disagree (%d refine, %d derefine)" % best)
