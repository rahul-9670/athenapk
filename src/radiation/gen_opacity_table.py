#!/usr/bin/env python3
"""Generate a tabulated multigroup opacity table for the AthenaPK M1 radiation package.

STATE-OF-THE-ART opacity architecture: the frequency-resolved dust+gas physics lives HERE
(offline, documented, swappable), and the C++ side only interpolates 2D (log rho, log T) tables.
Drop in a real Semenov (2003) / DSHARP / Draine monochromatic dataset by replacing `kappa_nu`.

What this produces, per (log10 rho_cgs, log10 T_K) grid point, in CODE units:
  - gray  kappa_P(rho,T), kappa_R(rho,T)          (Planck- and Rosseland-mean absorption)
  - per-group kappa_P,g(rho,T), kappa_R,g(rho,T)  (g = 0 .. n_group-1)
  - scattering kappa_s(rho,T)
all obtained by integrating a monochromatic kappa_nu = kappa_BL(rho,T) * psihat(nu,T) over the
Planck weight B_nu (emission/absorption) and the Rosseland weight dB_nu/dT (flux), per group and
over the full band.

Physics / honesty:
  * MAGNITUDE is anchored to Bell & Lin (1994) so the GRAY ROSSELAND reproduces the trusted law
    exactly (psihat is normalized so <psihat>_Rosseland,full = 1 => kappa_R,gray == kappa_BL).
    Bell & Lin's cold kappa ~ 2e-4 T^2 IS the Rosseland mean of a beta=2 dust, so this is
    physically self-consistent, not a fudge.
  * The frequency SHAPE psi(nu,T) is a multi-component modified-blackbody dust emissivity
    (nu^beta far-IR rise saturating to the geometric/gray limit above a grain-size break) with a
    MULTI-SPECIES SUBLIMATION CASCADE (ice mantle ~150 K, refractory silicate/carbon ~1500 K):
    above each front the corresponding opacity is removed and the shape blends toward gray gas.
  * NEW vs the analytic model it replaces: kappa_P is the TRUE Planck mean (not planck_ross_ratio
    * kappa_R), and the per-group means come from real frequency integration -- both tabulated.

This is a physically-motivated multi-component model, NOT literal Mie theory on measured optical
constants; the ARCHITECTURE accepts such a dataset unchanged (swap `kappa_nu`).
"""
import numpy as np
import struct
import argparse

# --- CGS constants ---
h_P = 6.62607015e-27      # erg s
k_B = 1.380649e-16        # erg/K
h_over_k = h_P / k_B      # s K  (= 4.7992e-11)

# --- shared FHC code units (must match units/physical_units.hpp defaults) ---
RHO_UNIT = 5.467e-19      # g/cm^3
LEN_UNIT = 2.81e16        # cm
OPACITY_UNIT = RHO_UNIT * LEN_UNIT   # kappa_code = kappa_cgs * OPACITY_UNIT  (= 1.536e-2)


# ---------------------------------------------------------------------------------------
def bell_lin_kappa(rho, T):
    """Bell & Lin (1994) Rosseland-mean gray opacity [cm^2/g], rho,T in cgs. Regime-skip-robust
    walk (matches C++ BellLinKappaFixed). 8 piecewise power laws kappa = k0 rho^a T^b."""
    k0 = np.array([2.0e-4, 2.0e16, 0.1, 2.0e81, 1.0e-8, 1.0e-36, 1.5e20, 0.348])
    aa = np.array([0.0, 0.0, 0.0, 1.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0])
    bb = np.array([2.0, -7.0, 0.5, -24.0, 3.0, 10.0, -2.5, 0.0])

    def Tcross(p, q):
        return ((k0[p] / k0[q]) * rho**(aa[p] - aa[q]))**(1.0 / (bb[q] - bb[p]))
    i = 0
    while i < 7:
        j = i + 1
        while j < 7 and Tcross(j, j + 1) <= Tcross(i, j):
            j += 1
        if T < Tcross(i, j):
            return k0[i] * rho**aa[i] * T**bb[i]
        i = j
    return k0[7] * rho**aa[7] * T**bb[7]


# ---------------------------------------------------------------------------------------
# Dust frequency-shape parameters (documented, swappable).
BETA = 1.8               # far-IR dust emissivity index (silicate/carbon aggregates; 1.5-2)
NU_BREAK = 1.0e14        # grain-size turnover [Hz] (lambda ~ few micron); flat (gray) above
# multi-species sublimation cascade: (fraction of dust opacity, T_sub_lo, T_sub_hi) [K]
DUST_SPECIES = [
    (0.35, 140.0, 170.0),   # volatile ice mantles (H2O/CO2) -> far-IR emissivity carriers
    (0.65, 1400.0, 1600.0), # refractory silicate + carbon grains
]


def dust_survival(T):
    """Surviving fraction of the dust far-IR opacity vs T (multi-species sublimation cascade)."""
    S = 0.0
    for frac, Tlo, Thi in DUST_SPECIES:
        if T <= Tlo:
            S += frac
        elif T < Thi:
            S += frac * (Thi - T) / (Thi - Tlo)
        # else 0
    return S


def psi(x, T):
    """Monochromatic opacity SHAPE at x = h nu / (k T): dust nu^beta rise saturating to the
    geometric (gray) limit above the grain break, blended toward gray gas as dust sublimates.
    Returns an UN-normalized shape; the caller divides by its Rosseland-weighted mean."""
    S = dust_survival(T)
    xb = h_over_k * NU_BREAK / T
    dust_shape = np.minimum((x / xb)**BETA, 1.0)
    return S * dust_shape + (1.0 - S)   # (1-S) -> gray gas continuum


# --- Planck / Rosseland spectral weights, with the common factor e^{-xa} pulled OUT of the
# ratio: the returned integrand carries e^{xa-x} (<= 1 on [xa,xb]) instead of e^{-x}, so it never
# overflows/underflows even for deep-Wien bands (xa >> 1). Every use divides two of these with the
# same xa, so the e^{xa} cancels exactly. ross=1 => Rosseland weight (x^4 dB/dT), else Planck (x^3 B).
def spectral_weight(x, xa, ross):
    em = np.exp(xa - x)               # e^{xa-x}, always in (0,1]
    d = 1.0 - np.exp(-x)
    return (x**4 * em / (d * d)) if ross else (x**3 * em / d)


def band_mean(T, xa, xb, ross, harmonic=False, nq=400):
    """<psi(^-1)>_weight over [xa,xb]. Simpson in x; e^{-xa} factored out of the num/den ratio."""
    xa = max(xa, 1e-8)
    if xb <= xa:
        return 1.0
    x = np.linspace(xa, xb, nq + 1)
    w = np.ones(nq + 1); w[1:-1:2] = 4; w[2:-1:2] = 2
    base = spectral_weight(x, xa, ross)
    s = psi(x, T)
    if harmonic:
        s = 1.0 / s
    num = np.sum(w * s * base); den = np.sum(w * base)
    val = num / (den + 1e-300)
    return 1.0 / val if harmonic else val


def group_opacities(rho, T, edges):
    """Return (kP_gray, kR_gray, [kP_g...], [kR_g...]) [cm^2/g] for kappa_nu = kBL * psihat."""
    ng = len(edges) - 1
    kBL = bell_lin_kappa(rho, T)
    XFULL = 80.0
    # normalize psi so <psihat>_Rosseland,full = 1  => kappa_R,gray == kBL (anchor to Bell&Lin).
    # ross=1 Rosseland weight, ross=0 Planck weight.
    psiR_full = band_mean(T, 1e-8, XFULL, ross=1)              # <psi>_Ross,full
    psiP_full = band_mean(T, 1e-8, XFULL, ross=0)             # <psi>_Planck,full
    kR_gray = kBL                                             # by construction
    kP_gray = kBL * (psiP_full / psiR_full)                  # TRUE Planck mean
    kP_g, kR_g = [], []
    for g in range(ng):
        xa = h_over_k * edges[g] / T
        xb = (xa + 60.0) if g == ng - 1 else h_over_k * edges[g + 1] / T
        xb = min(xb, xa + 60.0)
        # per-group Planck mean = kBL * <psi>_{P,band}/psiR_full ; Rosseland = kBL * harmonic mean
        pP = band_mean(T, xa, xb, ross=0)
        pR = band_mean(T, xa, xb, ross=1, harmonic=True)
        kP_g.append(kBL * pP / psiR_full)
        kR_g.append(kBL * pR / psiR_full)
    return kP_gray, kR_gray, kP_g, kR_g


def scattering_kappa(rho, T):
    """Gray electron-scattering opacity [cm^2/g] (Bell&Lin regime 7 floor, ionized gas)."""
    return 0.348 if T > 1.0e4 else 0.0


# ---------------------------------------------------------------------------------------
def build(path, nr, nT, rho_min, rho_max, T_min, T_max, edges):
    ng = len(edges) - 1
    lr = np.linspace(np.log10(rho_min), np.log10(rho_max), nr)   # log10 rho_cgs
    lT = np.linspace(np.log10(T_min), np.log10(T_max), nT)       # log10 T_K
    # GRAY magnitude tables over (rho,T) [code units]: true Planck mean kP, Rosseland kR, scatter ks.
    kP = np.zeros((nr, nT)); kR = np.zeros((nr, nT)); ks = np.zeros((nr, nT))
    # PER-GROUP MULTIPLIERS over (group,T) [dimensionless]: kappa_{P,R},g / kappa_{P,R},gray.
    # rho-INDEPENDENT (the kBL magnitude cancels in the ratio -> shape <psi>_band/<psi>_full),
    # so tabulated in T only. Evaluated at a reference rho (the ratio is rho-free by construction).
    mP = np.zeros((ng, nT)); mR = np.zeros((ng, nT))
    rho_ref = 1.0e-12
    for j in range(nT):
        T = 10.0**lT[j]
        gP, gR, pg, rg = group_opacities(rho_ref, T, edges)
        for g in range(ng):
            mP[g, j] = pg[g] / (gP + 1e-300)   # kappa_P,g / kappa_P,gray
            mR[g, j] = rg[g] / (gR + 1e-300)   # kappa_R,g / kappa_R,gray
    for i in range(nr):
        rho = 10.0**lr[i]
        for j in range(nT):
            T = 10.0**lT[j]
            gP, gR, _, _ = group_opacities(rho, T, edges)
            kP[i, j] = gP * OPACITY_UNIT
            kR[i, j] = gR * OPACITY_UNIT
            ks[i, j] = scattering_kappa(rho, T) * OPACITY_UNIT
    # binary: int64 [ng, nr, nT]; double [lr0,dlr,lT0,dlT];
    #   gray kP[nr,nT], kR[nr,nT], ks[nr,nT]  (row-major float64, code units)
    #   mult mP[ng,nT], mR[ng,nT]             (row-major float64, dimensionless)
    with open(path, "wb") as f:
        np.array([ng, nr, nT], dtype=np.int64).tofile(f)
        np.array([lr[0], lr[1]-lr[0], lT[0], lT[1]-lT[0]], dtype=np.float64).tofile(f)
        kP.astype(np.float64).tofile(f)
        kR.astype(np.float64).tofile(f)
        ks.astype(np.float64).tofile(f)
        mP.astype(np.float64).tofile(f)
        mR.astype(np.float64).tofile(f)
    print(f"wrote {path}: ng={ng} grid {nr}x{nT}  logrho[{lr[0]:.1f},{lr[-1]:.1f}] "
          f"logT[{lT[0]:.2f},{lT[-1]:.2f}]  OPACITY_UNIT={OPACITY_UNIT:.4e}")
    # sanity anchors (cgs gray Rosseland vs Bell&Lin; per-group multipliers)
    for T in [10.0, 100.0, 1000.0, 1.0e4]:
        kBL = bell_lin_kappa(1e-12, T)
        gP, gR, pg, rg = group_opacities(1e-12, T, edges)
        mp = [pg[g]/gP for g in range(ng)]
        print(f"  T={T:>7.0f}K: kR(BL)={kBL:.3e} kR_tab={gR:.3e} kP_tab={gP:.3e} "
              f"kP/kR={gP/max(gR,1e-30):.3f}  multP={[f'{m:.3f}' for m in mp]}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/beegfs/u/bbg6470/athenapk/src/radiation/opacity_table.bin")
    ap.add_argument("--nr", type=int, default=60)
    ap.add_argument("--nT", type=int, default=120)
    ap.add_argument("--rho_min", type=float, default=1e-20)
    ap.add_argument("--rho_max", type=float, default=1e-2)
    ap.add_argument("--T_min", type=float, default=3.0)
    ap.add_argument("--T_max", type=float, default=1e5)
    # default 3-group FHC structure (edges in Hz; 0 and +inf on the ends)
    ap.add_argument("--edges", default="0,1e12,1e15,1e30")
    args = ap.parse_args()
    edges = [float(x) for x in args.edges.split(",")]
    build(args.out, args.nr, args.nT, args.rho_min, args.rho_max,
          args.T_min, args.T_max, edges)
