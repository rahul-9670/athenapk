#!/usr/bin/env python
"""Bit-level dataset comparison of two parthenon HDF5 outputs.

Byte-comparing .phdf files is useless: parthenon embeds wall-clock/walltime metadata,
so two bit-identical physics runs always differ as bytes. This walks every dataset and
compares the raw bits (via .tobytes()), so +0.0/-0.0 and NaN payloads are not glossed over
the way `==` would.
"""
import sys, h5py, numpy as np

def walk(g, prefix=""):
    for k in g:
        item = g[k]
        path = f"{prefix}/{k}"
        if isinstance(item, h5py.Group):
            yield from walk(item, path)
        else:
            yield path, item

a, b = sys.argv[1], sys.argv[2]
# metadata that legitimately differs between two runs of the same physics
IGNORE = ("walltime", "WallTime", "Time/wall")

fa, fb = h5py.File(a, "r"), h5py.File(b, "r")
da = dict(walk(fa)); db = dict(walk(fb))

only_a = sorted(set(da) - set(db)); only_b = sorted(set(db) - set(da))
if only_a: print(f"  ONLY IN A: {only_a}")
if only_b: print(f"  ONLY IN B: {only_b}")

ndiff = nsame = nskip = 0
for path in sorted(set(da) & set(db)):
    if any(t in path for t in IGNORE):
        nskip += 1; continue
    x, y = da[path][()], db[path][()]
    xa, ya = np.asarray(x), np.asarray(y)
    if xa.shape != ya.shape:
        print(f"  SHAPE DIFF {path}: {xa.shape} vs {ya.shape}"); ndiff += 1; continue
    if xa.dtype.kind in "fiub" and ya.dtype.kind in "fiub":
        same = xa.tobytes() == ya.tobytes()          # bitwise, not ==
    else:
        same = np.array_equal(xa, ya)
    if same:
        nsame += 1
    else:
        ndiff += 1
        if xa.dtype.kind == "f":
            d = np.abs(xa.astype(float) - ya.astype(float))
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = np.nanmax(d / np.maximum(np.abs(xa.astype(float)), 1e-300))
            print(f"  DIFF {path}: max|abs|={np.nanmax(d):.6e} max|rel|={rel:.6e}")
        else:
            print(f"  DIFF {path}: (non-float)")

print(f"  -> {nsame} datasets bit-identical, {ndiff} differ, {nskip} skipped (walltime metadata)")
sys.exit(1 if ndiff else 0)
