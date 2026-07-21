#!/usr/bin/env python
"""Supplement: sample the isothermal era (files 10..140) + find rhocrit crossing bracket."""
import h5py, numpy as np, glob, json, time
RUN = '/beegfs/u/bbg6470/athenapk/runs/prod_t4_full'
OUT = '/home/bbg6470/.claude/jobs/873b4541/tmp'
files = sorted(glob.glob(f'{RUN}/parthenon.out1.[0-9]*.phdf'))
rows = []
for i in [10, 25, 40, 55, 70, 85, 100, 110, 120, 130, 140, 145]:
    t0 = time.time()
    with h5py.File(files[i], 'r') as f:
        t = float(f['Info'].attrs['Time']); c = int(f['Info'].attrs['NCycle'])
        lev = int(f['Levels'][:].max()); nb = int(f['Info'].attrs['NumMeshBlocks'])
        r = float(f['prim'][:, 0, :, :, :].max())
    rows.append(dict(t=t, rho_code=r, cyc=c, lev=lev, nb=nb))
    print(f'  [{i:4d}] t={t:.6f} cyc={c} rho={r:.4e} lev={lev} nb={nb} ({time.time()-t0:.0f}s)', flush=True)
with open(f'{OUT}/t4_early.json', 'w') as fo:
    json.dump(rows, fo, indent=1)
print('DONE')
