#!/usr/bin/env python
"""Content hashes of the 3D CT legs, for the N1 bit-identity regression.

Hashes DATASET CONTENT ONLY -- the phdf carries wall-clock metadata that legitimately differs
between runs of the same binary, so hashing the file bytes would report a spurious difference.
Same scheme as runs/audit_fix_regress/regress.sh.
"""
import sys, os, glob, hashlib
import numpy as np
import h5py

root = sys.argv[1] if len(sys.argv) > 1 else "out"
for d in sorted(glob.glob(os.path.join(root, "*_3d"))):
    name = os.path.basename(d)
    for f in sorted(glob.glob(os.path.join(d, "*.phdf"))):
        h = hashlib.sha256()
        with h5py.File(f, "r") as g:
            def walk(n, obj):
                if isinstance(obj, h5py.Dataset):
                    a = np.asarray(obj[()])
                    h.update(n.encode())
                    h.update(a.tobytes() if a.dtype != object else repr(a).encode())
            g.visititems(walk)
        print("%-16s %-10s %s" % (name, os.path.basename(f).split(".")[-2],
                                  h.hexdigest()[:32]))
