#!/usr/bin/env python
"""Extract prod_t4_full data for the star-formation flow-diagram artifact.
Outputs: /home/bbg6470/.claude/jobs/873b4541/tmp/t4_data.npz + t4_key_numbers.json
"""
import h5py, numpy as np, glob, json, os, sys

RUN = '/beegfs/u/bbg6470/athenapk/runs/prod_t4_full'
OUT = '/home/bbg6470/.claude/jobs/873b4541/tmp'
RHO0 = 5.467e-19      # g/cm^3 per code unit
V0   = 1.9e4          # cm/s
L0   = 2.81e16        # cm
T0_S = 1.48e12        # s  (l0/v0 = 1.4789e12)
AU   = 1.495979e13    # cm
MU   = 2.33
MH   = 1.6726e-24
KB   = 1.380649e-16
TCOEF = V0*V0*MU*MH/KB   # T[K] = (P/rho)_code * TCOEF  (mu=2.33 estimate)

files = sorted(glob.glob(f'{RUN}/parthenon.out1.[0-9]*.phdf'))
n = len(files)
# sample: every 20th, plus last 40 densely (runaway phase)
idx = sorted(set(list(range(0, n, 20)) + list(range(max(0, n-40), n))))
ts, rhomax, cyc, maxlev, nblk = [], [], [], [], []
for i in idx:
    with h5py.File(files[i], 'r') as f:
        t = float(f['Info'].attrs['Time'])
        c = int(f['Info'].attrs['NCycle'])
        lev = int(f['Levels'][:].max())
        nb = int(f['Info'].attrs['NumMeshBlocks'])
        r = float(f['prim'][:, 0, :, :, :].max())
    ts.append(t); rhomax.append(r); cyc.append(c); maxlev.append(lev); nblk.append(nb)
    if len(ts) % 20 == 0:
        print(f'  {len(ts)}/{len(idx)} sampled', flush=True)
ts = np.array(ts); rhomax = np.array(rhomax)

# ---- final snapshot deep dive ----
fin = files[-1]
with h5py.File(fin, 'r') as f:
    prim = f['prim']            # [blk, comp, k, j, i]
    rho = prim[:, 0, :, :, :]
    ib, ik, ij, ii = np.unravel_index(np.argmax(rho), rho.shape)
    xg = f['Locations/x'][:]; yg = f['Locations/y'][:]; zg = f['Locations/z'][:]
    xc = 0.5*(xg[:, 1:]+xg[:, :-1]); yc = 0.5*(yg[:, 1:]+yg[:, :-1]); zc = 0.5*(zg[:, 1:]+zg[:, :-1])
    peak = (float(xc[ib, ii]), float(yc[ib, ij]), float(zc[ib, ik]))
    P    = prim[:, 4, :, :, :]
    vx   = prim[:, 1, :, :, :]; vy = prim[:, 2, :, :, :]; vz = prim[:, 3, :, :, :]
    Bx   = prim[:, 5, :, :, :]; By = prim[:, 6, :, :, :]; Bz = prim[:, 7, :, :, :]
    xe   = prim[:, 13, :, :, :]   # prim_scalar_4 = comp 9(psi)+... comps: 0rho 1-3v 4P 5-7B 8psi 9-13 scalars
    Er   = f['rad.Er'][:, 0, :, :, :] if f['rad.Er'].ndim == 5 else f['rad.Er'][:]
    lev_f = f['Levels'][:]
    tfin = float(f['Info'].attrs['Time']); cfin = int(f['Info'].attrs['NCycle'])

    nb_, nk, nj, ni = rho.shape
    XX = np.broadcast_to(xc[:, None, None, :], rho.shape)
    YY = np.broadcast_to(yc[:, None, :, None], rho.shape)
    ZZ = np.broadcast_to(zc[:, :, None, None], rho.shape)
    dx_blk = (xg[:, 1] - xg[:, 0])  # per-block cell size

    central = dict(
        rho_code=float(rho[ib, ik, ij, ii]), P_code=float(P[ib, ik, ij, ii]),
        T_est_K=float(P[ib, ik, ij, ii]/rho[ib, ik, ij, ii]*TCOEF),
        Bmag_code=float(np.sqrt(Bx[ib, ik, ij, ii]**2+By[ib, ik, ij, ii]**2+Bz[ib, ik, ij, ii]**2)),
        xe=float(xe[ib, ik, ij, ii]), Er_code=float(Er[ib, ik, ij, ii]),
        peak_xyz=peak, dx_finest_code=float(dx_blk.min()))

    # radial profiles about the peak (log bins 0.02 AU .. 50000 AU)
    R = np.sqrt((XX-peak[0])**2+(YY-peak[1])**2+(ZZ-peak[2])**2)
    Vol = np.broadcast_to((dx_blk**3)[:, None, None, None], rho.shape)
    r_au = (R*L0/AU).ravel()
    rho_f = rho[:].ravel(); P_f = P[:].ravel(); w = Vol.ravel()
    B_f = np.sqrt(Bx**2+By**2+Bz**2).ravel()
    # radial velocity (infall)
    eps = 1e-30
    vr_f = (((XX-peak[0])*vx + (YY-peak[1])*vy + (ZZ-peak[2])*vz)/(R+eps)).ravel()
    xe_f = xe[:].ravel()
    rbins = np.logspace(np.log10(0.02), np.log10(5e4), 60)
    ridx = np.digitize(r_au, rbins)
    prof = {k: np.full(len(rbins)+1, np.nan) for k in ('rho', 'T', 'B', 'vr', 'xe')}
    cnt = np.zeros(len(rbins)+1)
    for k, q in (('rho', rho_f), ('T', P_f/np.maximum(rho_f, eps)*TCOEF),
                 ('B', B_f), ('vr', vr_f), ('xe', xe_f)):
        s = np.bincount(ridx, weights=q*w, minlength=len(rbins)+1)
        sw = np.bincount(ridx, weights=w, minlength=len(rbins)+1)
        prof[k] = s/np.maximum(sw, 1e-300)
        cnt = np.bincount(ridx, minlength=len(rbins)+1)
    rmid = np.sqrt(rbins[:-1]*rbins[1:])

    # B-rho relation: mass-weighted mean |B| in log-rho bins
    lr = np.log10(np.maximum(rho_f, 1e-10))
    rb = np.linspace(-1, np.log10(rho_f.max()), 50)
    bidx = np.digitize(lr, rb)
    mw = rho_f*w
    Bs = np.bincount(bidx, weights=B_f*mw, minlength=len(rb)+1)
    Ws = np.bincount(bidx, weights=mw, minlength=len(rb)+1)
    Bmean = Bs/np.maximum(Ws, 1e-300)
    rbmid = 0.5*(rb[:-1]+rb[1:])

    # midplane slices (cells whose z-range straddles z_peak), three zooms (half-width, code)
    slices = {}
    zlo = np.broadcast_to(zg[:, :-1][:, :, None, None], rho.shape)
    zhi = np.broadcast_to(zg[:, 1:][:, :, None, None], rho.shape)
    mid = ((zlo <= peak[2]) & (zhi > peak[2])).ravel()
    xs = XX.ravel()[mid]; ys = YY.ravel()[mid]; rs = rho_f[mid]; dxs = np.broadcast_to(dx_blk[:, None, None, None], rho.shape).ravel()[mid]
    for tag, hw in (('full', 26.0), ('mid', 0.5), ('core', 0.004)):
        m = (np.abs(xs-peak[0]) < hw) & (np.abs(ys-peak[1]) < hw)
        npx = 400
        gx = np.linspace(peak[0]-hw, peak[0]+hw, npx)
        gy = np.linspace(peak[1]-hw, peak[1]+hw, npx)
        img = np.zeros((npx, npx)); filled = np.zeros((npx, npx), bool)
        # paint coarsest first, finest last (finer overwrites)
        order = np.argsort(-dxs[m])
        xi = xs[m][order]; yi = ys[m][order]; ri = rs[m][order]; di = dxs[m][order]
        for x0, y0, r0, d0 in zip(xi, yi, ri, di):
            i0 = np.searchsorted(gx, x0-d0/2); i1 = np.searchsorted(gx, x0+d0/2)
            j0 = np.searchsorted(gy, y0-d0/2); j1 = np.searchsorted(gy, y0+d0/2)
            i1 = max(i1, i0+1); j1 = max(j1, j0+1)
            img[j0:j1, i0:i1] = r0; filled[j0:j1, i0:i1] = True
        img[~filled] = np.nan
        slices[tag] = dict(img=img, extent=[gx[0], gx[-1], gy[0], gy[-1]], hw=hw)
        print(f'  slice {tag}: {m.sum()} cells painted', flush=True)

np.savez_compressed(f'{OUT}/t4_data.npz',
    ts=ts, rhomax=rhomax, cyc=np.array(cyc), maxlev=np.array(maxlev), nblk=np.array(nblk),
    rmid=rmid, prof_rho=prof['rho'][1:len(rbins)], prof_T=prof['T'][1:len(rbins)],
    prof_B=prof['B'][1:len(rbins)], prof_vr=prof['vr'][1:len(rbins)], prof_xe=prof['xe'][1:len(rbins)],
    rbmid=rbmid, Bmean=Bmean[1:len(rb)],
    slice_full=slices['full']['img'], ext_full=slices['full']['extent'],
    slice_mid=slices['mid']['img'], ext_mid=slices['mid']['extent'],
    slice_core=slices['core']['img'], ext_core=slices['core']['extent'])

key = dict(t_final_code=tfin, cycle_final=cfin, n_snapshots=n,
           central=central, rho_max_cgs=central['rho_code']*RHO0,
           T_central_est_K=central['T_est_K'],
           dx_finest_cm=central['dx_finest_code']*L0,
           dx_finest_au=central['dx_finest_code']*L0/AU)
with open(f'{OUT}/t4_key_numbers.json', 'w') as f:
    json.dump(key, f, indent=1)
print(json.dumps(key, indent=1))
print('DONE')
