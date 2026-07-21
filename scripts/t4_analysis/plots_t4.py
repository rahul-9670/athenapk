#!/usr/bin/env python
"""Figures for the flow-diagram artifact, from t4_data.npz. Dark theme to match the page."""
import numpy as np, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

OUT = '/home/bbg6470/.claude/jobs/873b4541/tmp'
D = np.load(f'{OUT}/t4_data.npz')
RHO0, V0, L0 = 5.467e-19, 1.9e4, 2.81e16
AU = 1.495979e13
T0YR = 1.4789e12/3.1557e7      # code time unit in yr = 46,862 yr
BUNIT_G = V0*np.sqrt(4*np.pi*RHO0)   # Gaussian-equivalent field per code unit = 4.98e-5 G

# palette (dark ground #0C0F16)
INK, MUTED, GRID = '#E9ECF4', '#97A0B5', '#232B3E'
ACC = '#E0A33C'      # amber accent (single-series line)
CO2 = '#5B8DD9'      # cool blue secondary
plt.rcParams.update({
    'figure.facecolor': '#0C0F16', 'axes.facecolor': '#10141F',
    'axes.edgecolor': GRID, 'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': MUTED, 'ytick.color': MUTED, 'grid.color': GRID,
    'font.size': 11, 'axes.titlesize': 12, 'font.family': 'DejaVu Sans',
    'axes.grid': True, 'grid.linewidth': 0.5, 'grid.alpha': 0.6,
    'axes.linewidth': 0.8, 'savefig.facecolor': '#0C0F16', 'savefig.dpi': 160})

ts, rho = D['ts'], D['rhomax']*RHO0
# merge early-phase supplement if present
import glob as _g
for fn in sorted(_g.glob(f'{OUT}/t4_early*.json')):
    early = json.load(open(fn))
    te = np.array([r['t'] for r in early]); re_ = np.array([r['rho_code'] for r in early])*RHO0
    ts = np.concatenate([ts, te]); rho = np.concatenate([rho, re_])
o = np.argsort(ts); ts, rho = ts[o], rho[o]
tyr = ts*T0YR

# ---- F1: rho_max(t), full + runaway zoom ----
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.4, 3.9), constrained_layout=True)
a1.semilogy(tyr/1e3, rho, '-', color=ACC, lw=2, marker='o', ms=4, mfc=ACC, mec='none')
a1.set_xlabel('time  [kyr]'); a1.set_ylabel(r'max density  [g cm$^{-3}$]')
a1.set_title('Central density over the whole run', loc='left', color=INK)
for y, lab in ((1e-13, 'first core forms (ρ$_{crit}$)'), (1.65e-8, 'runaway resumes (L13)')):
    a1.axhline(y, color=MUTED, lw=0.7, ls='--', alpha=0.8)
    a1.text(0.4, y*1.8, lab, fontsize=8.5, color=MUTED)
# zoom: last samples vs time-before-end (log)
tend = ts[-1]
m = rho > 1e-12
dtb = np.maximum((tend - ts[m])*T0YR, 1e-3)  # yr before stop
a2.loglog(dtb, rho[m], '-', color=ACC, lw=2, marker='o', ms=4, mec='none')
a2.invert_xaxis()
a2.set_xlabel('time before stop  [yr]  (log, runs right→left)')
a2.set_title('The runaway: first core → second collapse', loc='left', color=INK)
a2.axhline(1.65e-8, color=MUTED, lw=0.7, ls='--', alpha=0.8)
fig.savefig(f'{OUT}/f1_rhomax.png'); plt.close(fig)

# ---- F2: density slices ----
fig, axs = plt.subplots(1, 3, figsize=(12.6, 4.3), constrained_layout=True)
meta = [('slice_full', 'ext_full', 'Full box · 0.47 pc'),
        ('slice_mid', 'ext_mid', 'Envelope · ±939 AU'),
        ('slice_core', 'ext_core', 'First core & runaway centre · ±7.5 AU')]
for ax, (k, ek, title) in zip(axs, meta):
    img = D[k]*RHO0; ext = D[ek]
    ext_au = [(e - (ext[0]+ext[1])/2)*L0/AU for e in ext[:2]] + \
             [(e - (ext[2]+ext[3])/2)*L0/AU for e in ext[2:]]
    vmin = np.nanpercentile(img, 2); vmax = np.nanmax(img)
    im = ax.imshow(img, origin='lower', extent=ext_au, cmap='inferno',
                   norm=LogNorm(vmin=max(vmin, 1e-22), vmax=vmax))
    ax.set_title(title, loc='left', fontsize=10.5, color=INK)
    ax.set_xlabel('AU'); ax.grid(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label(r'ρ [g cm$^{-3}$]', fontsize=9); cb.ax.tick_params(labelsize=8)
fig.savefig(f'{OUT}/f2_slices.png'); plt.close(fig)

# ---- F3: radial profiles 2x2 ----
r = D['rmid']
fig, axs = plt.subplots(2, 2, figsize=(10.4, 7.2), constrained_layout=True)
panels = [
    (axs[0, 0], D['prof_rho']*RHO0, r'ρ  [g cm$^{-3}$]', 'log'),
    (axs[0, 1], D['prof_T'], 'T (μ=2.33 estimate)  [K]', 'log'),
    (axs[1, 0], D['prof_B']*BUNIT_G*1e6, '|B|  [μG, Gaussian]', 'log'),
    (axs[1, 1], D['prof_vr']*V0/1e5, r'v$_r$  [km s$^{-1}$]', 'linear'),
]
for ax, q, lab, sc in panels:
    ax.plot(r, q, '-', color=ACC, lw=2)
    ax.set_xscale('log'); ax.set_yscale(sc)
    ax.set_xlabel('radius from density peak  [AU]'); ax.set_ylabel(lab)
    if sc == 'linear':
        ax.axhline(0, color=MUTED, lw=0.7)
axs[0, 0].plot(r[(r > 20) & (r < 2e4)], D['prof_rho'][(r > 20) & (r < 2e4)].max()*RHO0
               * ((r[(r > 20) & (r < 2e4)]/r[(r > 20) & (r < 2e4)][0])**-2), '--',
               color=CO2, lw=1.2)
axs[0, 0].text(0.05, 0.08, '– – ρ ∝ r⁻² (isothermal envelope)', transform=axs[0, 0].transAxes,
               fontsize=9, color=CO2)
fig.suptitle('Final snapshot (cycle 77250) — spherically averaged profiles about the density peak',
             color=INK, fontsize=12, x=0.02, ha='left')
fig.savefig(f'{OUT}/f3_profiles.png'); plt.close(fig)

# ---- F4: B–rho relation ----
fig, ax = plt.subplots(figsize=(6.4, 4.4), constrained_layout=True)
rr = 10**D['rbmid']*RHO0
bb = D['Bmean']*BUNIT_G*1e6
mm = bb > 0
ax.loglog(rr[mm], bb[mm], '-', color=ACC, lw=2)
x0 = 3e-17
for kap, style in ((0.5, '--'), (2/3., ':')):
    xx = np.logspace(-17, -8, 40)
    i0 = np.argmin(np.abs(rr[mm]-x0))
    yy = bb[mm][i0]*(xx/x0)**kap
    ax.loglog(xx, yy, style, color=CO2, lw=1.2)
ax.text(0.62, 0.28, 'B ∝ ρ$^{1/2}$', transform=ax.transAxes, color=CO2, fontsize=10)
ax.text(0.62, 0.16, 'B ∝ ρ$^{2/3}$', transform=ax.transAxes, color=CO2, fontsize=10)
ax.set_xlabel(r'ρ  [g cm$^{-3}$]'); ax.set_ylabel('mass-weighted ⟨|B|⟩  [μG]')
ax.set_title('Magnetic field–density relation (final snapshot)', loc='left', color=INK)
fig.savefig(f'{OUT}/f4_brho.png'); plt.close(fig)

print(json.dumps(dict(
    t_end_kyr=float(tyr[-1]/1e3),
    rho_end=float(rho[-1]),
    B_central_uG=float(D['prof_B'][0]*BUNIT_G*1e6),
    T_profile_max=float(np.nanmax(D['prof_T'])),
    vr_min_kms=float(np.nanmin(D['prof_vr'])*V0/1e5)), indent=1))
print('PLOTS DONE')
