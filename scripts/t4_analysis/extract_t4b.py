#!/usr/bin/env python
"""Lean extraction for the flow-diagram artifact (v2: fewer samples, streaming final read)."""
import h5py, numpy as np, glob, json, time

RUN = '/beegfs/u/bbg6470/athenapk/runs/prod_t4_full'
OUT = '/home/bbg6470/.claude/jobs/873b4541/tmp'
RHO0, V0, L0 = 5.467e-19, 1.9e4, 2.81e16
AU = 1.495979e13
MU, MH, KB = 2.33, 1.6726e-24, 1.380649e-16
TCOEF = V0*V0*MU*MH/KB

files = sorted(glob.glob(f'{RUN}/parthenon.out1.[0-9]*.phdf'))
n = len(files)
idx = [0, 150, 300, 450, 600, 750, 900, 1050, 1200, 1350, 1440,
       1480, 1492, 1501, 1510, 1520, 1530, 1536, 1541, 1545]
idx = [i for i in idx if i < n]
ts, rhomax, cyc, maxlev, nblk = [], [], [], [], []
for i in idx:
    t0 = time.time()
    with h5py.File(files[i], 'r') as f:
        t = float(f['Info'].attrs['Time']); c = int(f['Info'].attrs['NCycle'])
        lev = int(f['Levels'][:].max()); nb = int(f['Info'].attrs['NumMeshBlocks'])
        r = float(f['prim'][:, 0, :, :, :].max())
    ts.append(t); rhomax.append(r); cyc.append(c); maxlev.append(lev); nblk.append(nb)
    print(f'  [{i:5d}] t={t:.6f} cyc={c} rho={r:.4e} lev={lev} nb={nb} ({time.time()-t0:.1f}s)', flush=True)
ts = np.array(ts); rhomax = np.array(rhomax)

# ---- final snapshot deep dive: stream the whole prim once ----
fin = files[-1]
print('streaming final prim ...', flush=True)
t0 = time.time()
with h5py.File(fin, 'r') as f:
    prim = f['prim'][:]          # [blk, comp, k, j, i] float32, ~3.7 GB
    Er = f['rad.Er'][:]
    xg = f['Locations/x'][:]; yg = f['Locations/y'][:]; zg = f['Locations/z'][:]
    tfin = float(f['Info'].attrs['Time']); cfin = int(f['Info'].attrs['NCycle'])
print(f'  prim {prim.shape} {prim.dtype} in {time.time()-t0:.0f}s', flush=True)
Er = Er.reshape(Er.shape[0], *prim.shape[2:]) if Er.ndim == 5 else Er

rho = prim[:, 0]; vx, vy, vz = prim[:, 1], prim[:, 2], prim[:, 3]
P = prim[:, 4]; Bx, By, Bz = prim[:, 5], prim[:, 6], prim[:, 7]
xe = prim[:, 13]
ib, ik, ij, ii = np.unravel_index(np.argmax(rho), rho.shape)
xc = 0.5*(xg[:, 1:]+xg[:, :-1]); yc = 0.5*(yg[:, 1:]+yg[:, :-1]); zc = 0.5*(zg[:, 1:]+zg[:, :-1])
peak = (float(xc[ib, ii]), float(yc[ib, ij]), float(zc[ib, ik]))
dx_blk = (xg[:, 1] - xg[:, 0]).astype(np.float64)

central = dict(
    rho_code=float(rho[ib, ik, ij, ii]), P_code=float(P[ib, ik, ij, ii]),
    T_est_K=float(P[ib, ik, ij, ii]/rho[ib, ik, ij, ii]*TCOEF),
    Bmag_code=float(np.sqrt(Bx[ib, ik, ij, ii]**2+By[ib, ik, ij, ii]**2+Bz[ib, ik, ij, ii]**2)),
    xe=float(xe[ib, ik, ij, ii]), Er_code=float(Er[ib, ik, ij, ii]),
    peak_xyz=peak, dx_finest_code=float(dx_blk.min()))
print('central:', central, flush=True)

XX = np.broadcast_to(xc[:, None, None, :], rho.shape)
YY = np.broadcast_to(yc[:, None, :, None], rho.shape)
ZZ = np.broadcast_to(zc[:, :, None, None], rho.shape)
R = np.sqrt((XX-peak[0])**2+(YY-peak[1])**2+(ZZ-peak[2])**2)
Vol = np.broadcast_to((dx_blk**3)[:, None, None, None], rho.shape)
eps = 1e-30
r_au = (R*L0/AU).ravel()
rho_f = rho.ravel().astype(np.float64); P_f = P.ravel().astype(np.float64)
w = Vol.ravel()
B_f = np.sqrt(Bx.astype(np.float64)**2+By**2+Bz**2).ravel()
vr_f = (((XX-peak[0])*vx + (YY-peak[1])*vy + (ZZ-peak[2])*vz)/(R+eps)).ravel().astype(np.float64)
xe_f = xe.ravel().astype(np.float64)

rbins = np.logspace(np.log10(0.02), np.log10(5e4), 60)
ridx = np.digitize(r_au, rbins)
prof = {}
sw = np.bincount(ridx, weights=w, minlength=len(rbins)+1)
for k, q in (('rho', rho_f), ('T', P_f/np.maximum(rho_f, eps)*TCOEF),
             ('B', B_f), ('vr', vr_f), ('xe', xe_f)):
    s = np.bincount(ridx, weights=q*w, minlength=len(rbins)+1)
    prof[k] = (s/np.maximum(sw, 1e-300))
rmid = np.sqrt(rbins[:-1]*rbins[1:])
print('profiles done', flush=True)

lr = np.log10(np.maximum(rho_f, 1e-10))
rb = np.linspace(-1, float(lr.max()), 50)
bidx = np.digitize(lr, rb)
mw = rho_f*w
Bs = np.bincount(bidx, weights=B_f*mw, minlength=len(rb)+1)
Ws = np.bincount(bidx, weights=mw, minlength=len(rb)+1)
Bmean = Bs/np.maximum(Ws, 1e-300)
rbmid = 0.5*(rb[:-1]+rb[1:])
print('B-rho done', flush=True)

# midplane slices: per-BLOCK painting (fast), finest painted last
zlo_b = zg[:, :-1]; zhi_b = zg[:, 1:]
slices = {}
for tag, hw in (('full', 26.0), ('mid', 0.5), ('core', 0.004)):
    npx = 480
    gx = np.linspace(peak[0]-hw, peak[0]+hw, npx)
    gy = np.linspace(peak[1]-hw, peak[1]+hw, npx)
    img = np.full((npx, npx), np.nan)
    order = np.argsort(-dx_blk)   # coarse first
    for b in order:
        kz = np.where((zlo_b[b] <= peak[2]) & (zhi_b[b] > peak[2]))[0]
        if len(kz) == 0: continue
        if xc[b].max() < gx[0] or xc[b].min() > gx[-1]: continue
        if yc[b].max() < gy[0] or yc[b].min() > gy[-1]: continue
        plane = rho[b, kz[0], :, :]          # [j, i]
        d = dx_blk[b]
        i0 = np.clip(np.searchsorted(gx, xc[b]-d/2), 0, npx)
        i1 = np.clip(np.searchsorted(gx, xc[b]+d/2), 0, npx)
        j0 = np.clip(np.searchsorted(gy, yc[b]-d/2), 0, npx)
        j1 = np.clip(np.searchsorted(gy, yc[b]+d/2), 0, npx)
        for j in range(32):
            if j1[j] <= j0[j]: continue
            for i in range(32):
                if i1[i] > i0[i]:
                    img[j0[j]:j1[j], i0[i]:i1[i]] = plane[j, i]
        # ensure at least 1px for sub-pixel cells: paint block bbox center row/col handled above
    slices[tag] = dict(img=img, extent=[float(gx[0]), float(gx[-1]), float(gy[0]), float(gy[-1])])
    print(f'  slice {tag} done, filled {np.isfinite(img).mean()*100:.0f}%', flush=True)

np.savez_compressed(f'{OUT}/t4_data.npz',
    ts=ts, rhomax=rhomax, cyc=np.array(cyc), maxlev=np.array(maxlev), nblk=np.array(nblk),
    rmid=rmid, prof_rho=prof['rho'][1:len(rbins)], prof_T=prof['T'][1:len(rbins)],
    prof_B=prof['B'][1:len(rbins)], prof_vr=prof['vr'][1:len(rbins)], prof_xe=prof['xe'][1:len(rbins)],
    rbmid=rbmid, Bmean=Bmean[1:len(rb)],
    slice_full=slices['full']['img'], ext_full=slices['full']['extent'],
    slice_mid=slices['mid']['img'], ext_mid=slices['mid']['extent'],
    slice_core=slices['core']['img'], ext_core=slices['core']['extent'])

key = dict(t_final_code=tfin, cycle_final=cfin, n_snapshots=n, central=central,
           rho_max_cgs=central['rho_code']*RHO0, T_central_est_K=central['T_est_K'],
           dx_finest_cm=central['dx_finest_code']*L0, dx_finest_au=central['dx_finest_code']*L0/AU,
           tseries=[dict(t=float(a), rho_code=float(b), cyc=int(c), lev=int(d), nb=int(e))
                    for a, b, c, d, e in zip(ts, rhomax, cyc, maxlev, nblk)])
with open(f'{OUT}/t4_key_numbers.json', 'w') as fo:
    json.dump(key, fo, indent=1)
print('DONE', flush=True)
