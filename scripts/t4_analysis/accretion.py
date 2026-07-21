#!/usr/bin/env python
"""Final-snapshot accretion diagnostics for the professor email:
Mdot(r) shell mass flux, magnetic flux Phi(r) through the equatorial disk about L_hat,
Phi-dot(r) advection estimate, core mass-to-flux ratio, B-rho power-law fit."""
import h5py, numpy as np, json

RUN = '/beegfs/u/bbg6470/athenapk/runs/prod_t4_full/parthenon.out1.01545.phdf'
OUT = '/home/bbg6470/.claude/jobs/873b4541/tmp'
RHO0, V0, L0 = 5.467e-19, 1.9e4, 2.81e16
T0S = L0/V0; AU = 1.495979e13
BUNIT = V0*np.sqrt(4*np.pi*RHO0)     # Gaussian-equivalent G per code B
MSUN, YR = 1.989e33, 3.1557e7
MUNIT = RHO0*L0**3                   # g per code mass

f = h5py.File(RUN, 'r')
rho = f['prim'][:, 0, :, :, :]
xg = f['Locations/x'][:]; yg = f['Locations/y'][:]; zg = f['Locations/z'][:]
xc = 0.5*(xg[:, 1:]+xg[:, :-1]); yc = 0.5*(yg[:, 1:]+yg[:, :-1]); zc = 0.5*(zg[:, 1:]+zg[:, :-1])
dxb = xg[:, 1]-xg[:, 0]
ib, ik, ij, ii = np.unravel_index(np.argmax(rho), rho.shape)
pk = (xc[ib, ii], yc[ib, ij], zc[ib, ik])
RMAX_AU = 50.0; RMAX = RMAX_AU*AU/L0

# blocks overlapping the analysis sphere
bx = np.abs(0.5*(xg[:, 0]+xg[:, -1]) - pk[0]) - 0.5*(xg[:, -1]-xg[:, 0])
by = np.abs(0.5*(yg[:, 0]+yg[:, -1]) - pk[1]) - 0.5*(yg[:, -1]-yg[:, 0])
bz = np.abs(0.5*(zg[:, 0]+zg[:, -1]) - pk[2]) - 0.5*(zg[:, -1]-zg[:, 0])
bsel = np.where((np.maximum(bx, 0)**2 + np.maximum(by, 0)**2 + np.maximum(bz, 0)**2) < RMAX**2)[0]
print(f'{len(bsel)} blocks in r<{RMAX_AU} AU', flush=True)

X = []; Y = []; Z = []; RH = []; VX = []; VY = []; VZ = []; BX = []; BY = []; BZ = []; VOL = []
for b in bsel:
    p = f['prim'][b]                    # [14,k,j,i]
    kk = np.broadcast_to(zc[b][:, None, None], p[0].shape)
    jj = np.broadcast_to(yc[b][None, :, None], p[0].shape)
    iix = np.broadcast_to(xc[b][None, None, :], p[0].shape)
    m = ((iix-pk[0])**2 + (jj-pk[1])**2 + (kk-pk[2])**2) < RMAX**2
    X.append((iix-pk[0])[m]); Y.append((jj-pk[1])[m]); Z.append((kk-pk[2])[m])
    RH.append(p[0][m]); VX.append(p[1][m]); VY.append(p[2][m]); VZ.append(p[3][m])
    BX.append(p[5][m]); BY.append(p[6][m]); BZ.append(p[7][m])
    VOL.append(np.full(m.sum(), dxb[b]**3))
X = np.concatenate(X); Y = np.concatenate(Y); Z = np.concatenate(Z)
RH = np.concatenate(RH).astype(np.float64); VOL = np.concatenate(VOL)
VX = np.concatenate(VX); VY = np.concatenate(VY); VZ = np.concatenate(VZ)
BX = np.concatenate(BX).astype(np.float64); BY = np.concatenate(BY); BZ = np.concatenate(BZ)
R = np.sqrt(X*X+Y*Y+Z*Z); M = RH*VOL
print(f'{len(R)} cells', flush=True)

# core bulk velocity (rho>rho_max/100) and L_hat
cm = RH > RH.max()/100
vcm = np.array([(M*VX)[cm].sum(), (M*VY)[cm].sum(), (M*VZ)[cm].sum()])/M[cm].sum()
ux, uy, uz = VX-vcm[0], VY-vcm[1], VZ-vcm[2]
Lv = np.array([ (M*(Y*uz-Z*uy))[cm].sum(), (M*(Z*ux-X*uz))[cm].sum(), (M*(X*uy-Y*ux))[cm].sum() ])
lh = Lv/np.linalg.norm(Lv)
vr = (X*ux + Y*uy + Z*uz)/np.maximum(R, 1e-12)

res = {'L_hat': [float(v) for v in lh], 'vcm_code': [float(v) for v in vcm]}
# Mdot(r): volume-weighted shell average of 4 pi r^2 rho vr
print('r[AU]  Mdot[Msun/yr]  Phi[G cm^2]  Phidot[G cm^2/yr]  lambda_core', flush=True)
rows = []
for rau in (0.5, 1.0, 3.0, 10.0, 30.0):
    rc = rau*AU/L0
    sh = np.abs(R-rc) < 0.15*rc
    if sh.sum() < 20: continue
    w = VOL[sh]
    rvr = (RH*vr)[sh]
    mdot_code = -4*np.pi*rc*rc*(w*rvr).sum()/w.sum()
    mdot = mdot_code * MUNIT/T0S * YR/MSUN            # Msun/yr
    # magnetic flux through equatorial disk of radius rc about L_hat (midplane cells)
    zpar = X*lh[0]+Y*lh[1]+Z*lh[2]
    dxl = VOL**(1./3.)
    mid = (np.abs(zpar) < dxl/2) & (np.sqrt(np.maximum(R*R-zpar*zpar, 0)) < rc)
    Bpar = BX*lh[0]+BY*lh[1]+BZ*lh[2]
    area = (VOL[mid]/dxl[mid])                         # dx^2 per midplane cell
    phi_code = (Bpar[mid]*area).sum()
    phi = phi_code * BUNIT * L0*L0                     # G cm^2
    # flux advection: ring average of 2 pi r |B_pol x ...| ~ 2 pi r <vr Bpar> (estimate)
    ring = mid & (np.abs(np.sqrt(np.maximum(R*R-zpar*zpar, 0))-rc) < 0.15*rc)
    if ring.sum() > 5:
        phidot_code = -2*np.pi*rc*np.average((vr*Bpar)[ring], weights=VOL[ring])
        phidot = phidot_code * BUNIT*L0*L0/T0S * YR    # G cm^2 / yr
    else:
        phidot = float('nan')
    # enclosed mass and mass-to-flux
    enc = R < rc
    Menc = M[enc].sum()*MUNIT
    lam = (Menc/max(abs(phi), 1e-30)) * (2*np.pi*np.sqrt(6.674e-8))  # / (1/2pi sqrt(G))
    rows.append(dict(r_au=rau, mdot_msun_yr=float(mdot), phi_Gcm2=float(phi),
                     phidot_Gcm2_yr=float(phidot), Menc_msun=float(Menc/MSUN),
                     lambda_norm=float(lam)))
    print(f'{rau:5.1f}  {mdot:.3e}  {phi:.3e}  {phidot:.3e}  {lam:.1f}', flush=True)
res['shells'] = rows

# B-rho power-law fit from the saved npz (flux-freezing regime 1e-16..1e-12)
D = np.load(f'{OUT}/t4_data.npz')
lr = D['rbmid']; lB = np.log10(np.maximum(D['Bmean']*BUNIT*1e6, 1e-30))
rho_c = 10**lr*RHO0
mfit = (rho_c > 1e-16) & (rho_c < 1e-12) & np.isfinite(lB)
k1, b1 = np.polyfit(np.log10(rho_c[mfit]), lB[mfit], 1)
mplat = (rho_c > 3e-12) & np.isfinite(lB)
k2, b2 = np.polyfit(np.log10(rho_c[mplat]), lB[mplat], 1)
res['kappa_fluxfreeze'] = float(k1); res['kappa_plateau'] = float(k2)
print(f'B~rho^k: flux-freezing regime k={k1:.3f}; plateau regime k={k2:.3f}', flush=True)

json.dump(res, open(f'{OUT}/accretion.json', 'w'), indent=1)
print('DONE')
