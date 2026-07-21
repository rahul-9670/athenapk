#!/usr/bin/env python
"""Core rotation history of prod_t4_full: Omega(t) mass-weighted about the density peak,
for two core definitions (comoving rho > rho_max/100, and fixed rho > 1e-12 g/cm3).
N_rot = trapz(Omega dt)/2pi. Velocities read only for blocks containing core cells."""
import h5py, numpy as np, glob, json, time

RUN = '/beegfs/u/bbg6470/athenapk/runs/prod_t4_full'
OUT = '/home/bbg6470/.claude/jobs/873b4541/tmp'
RHO0, V0, L0 = 5.467e-19, 1.9e4, 2.81e16
T0S = L0 / V0                       # 1.4789e12 s
AU = 1.495979e13
FIX_CODE = 1e-12 / RHO0             # fixed threshold in code units (1.829e6)

snaps = [25, 55, 85, 120, 150, 300, 600, 900, 1200, 1350, 1440, 1480,
         1501, 1520, 1536, 1545]
files = sorted(glob.glob(f'{RUN}/parthenon.out1.[0-9]*.phdf'))

def core_omega(f, mask, rho, xc, yc, zc, dxb, pk):
    """mass-weighted Omega about peak pk for cells in mask (block-selective v read)."""
    bsel = np.unique(np.where(mask.any(axis=(1, 2, 3)))[0])
    n = 0; arrs = []
    for b in bsel:
        v = f['prim'][b, 1:4, :, :, :]          # [3,k,j,i]
        m = mask[b]
        if not m.any(): continue
        kk, jj, ii = np.where(m)
        x = xc[b][ii] - pk[0]; y = yc[b][jj] - pk[1]; z = zc[b][kk] - pk[2]
        mm = rho[b][m] * dxb[b]**3
        arrs.append((x, y, z, v[0][m], v[1][m], v[2][m], mm))
    X = np.concatenate([a[0] for a in arrs]); Y = np.concatenate([a[1] for a in arrs])
    Z = np.concatenate([a[2] for a in arrs])
    VX = np.concatenate([a[3] for a in arrs]); VY = np.concatenate([a[4] for a in arrs])
    VZ = np.concatenate([a[5] for a in arrs]); M = np.concatenate([a[6] for a in arrs])
    Mt = M.sum()
    vcm = np.array([(M*VX).sum(), (M*VY).sum(), (M*VZ).sum()]) / Mt
    ux, uy, uz = VX-vcm[0], VY-vcm[1], VZ-vcm[2]
    Lx = (M*(Y*uz - Z*uy)).sum(); Ly = (M*(Z*ux - X*uz)).sum()
    Lz = (M*(X*uy - Y*ux)).sum()
    L = np.array([Lx, Ly, Lz]); Lmag = np.linalg.norm(L)
    lh = L / (Lmag + 1e-300)
    # moment of inertia about the L axis
    rp2 = (X*X+Y*Y+Z*Z) - (X*lh[0]+Y*lh[1]+Z*lh[2])**2
    I = (M*rp2).sum()
    om = Lmag / (I + 1e-300)                     # rad per code time
    rrms = np.sqrt((M*(X*X+Y*Y+Z*Z)).sum()/Mt)   # code
    return om, Mt, rrms, lh, len(M)

rows = []
for i in snaps:
    t0 = time.time()
    with h5py.File(files[i], 'r') as f:
        t = float(f['Info'].attrs['Time']); c = int(f['Info'].attrs['NCycle'])
        rho = f['prim'][:, 0, :, :, :]
        xg = f['Locations/x'][:]; yg = f['Locations/y'][:]; zg = f['Locations/z'][:]
        xc = 0.5*(xg[:, 1:]+xg[:, :-1]); yc = 0.5*(yg[:, 1:]+yg[:, :-1])
        zc = 0.5*(zg[:, 1:]+zg[:, :-1])
        dxb = xg[:, 1] - xg[:, 0]
        ib, ik, ij, ii = np.unravel_index(np.argmax(rho), rho.shape)
        pk = (xc[ib, ii], yc[ib, ij], zc[ib, ik])
        rmax = rho.max()
        row = dict(snap=i, t=t, cyc=c, rho_max_cgs=float(rmax*RHO0))
        for tag, thr in (('com', rmax/100.0), ('fix', FIX_CODE)):
            mask = rho > thr
            if mask.sum() < 8:
                row[tag] = None; continue
            om, Mt, rrms, lh, ncell = core_omega(f, mask, rho, xc, yc, zc, dxb, pk)
            row[tag] = dict(omega_code=float(om), M_code=float(Mt),
                            r_rms_au=float(rrms*L0/AU), ncell=int(ncell),
                            Lhat=[float(v) for v in lh])
        rows.append(row)
    o1 = row['com']['omega_code'] if row['com'] else float('nan')
    o2 = row['fix']['omega_code'] if row['fix'] else float('nan')
    print(f"[{i:5d}] t={t:.6f} rho={rmax*RHO0:.2e} Om_com={o1:.3e} Om_fix={o2:.3e} "
          f"({time.time()-t0:.0f}s)", flush=True)

# integrate rotations
out = dict(rows=rows)
for tag in ('com', 'fix'):
    ts = np.array([r['t'] for r in rows if r[tag]])
    om = np.array([r[tag]['omega_code'] for r in rows if r[tag]])
    nrot = np.trapezoid(om, ts) / (2*np.pi)
    out[f'nrot_{tag}'] = float(nrot)
    last = [r for r in rows if r[tag]][-1]
    om_s = last[tag]['omega_code'] / T0S
    out[f'final_{tag}'] = dict(omega_rad_s=om_s, period_yr=2*np.pi/om_s/3.1557e7,
                               M_msun=last[tag]['M_code']*RHO0*L0**3/1.989e33,
                               r_rms_au=last[tag]['r_rms_au'])
    print(f"{tag}: N_rot = {nrot:.3f}; final period "
          f"{out[f'final_{tag}']['period_yr']:.3f} yr; "
          f"M = {out[f'final_{tag}']['M_msun']:.4f} Msun; "
          f"r_rms = {out[f'final_{tag}']['r_rms_au']:.3f} AU", flush=True)
with open(f'{OUT}/rotations.json', 'w') as fo:
    json.dump(out, fo, indent=1)
print('DONE')
