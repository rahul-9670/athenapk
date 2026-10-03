#!/usr/bin/env python3
# ========================================================================================
# AthenaPK - a performance portable block structured AMR astrophysical MHD code.
# Copyright (c) 2026, Athena-Parthenon Collaboration. All rights reserved.
# Licensed under the BSD 3-Clause License (the "LICENSE").
# ========================================================================================
"""
Generate the tabulated hydrogen/helium equation of state read by AthenaPK with
`<hydro> eos = hydrogen` (see src/eos/eos_table.hpp and docs/eos.md).

Physics (gas = hydrogen mass fraction X plus helium mass fraction Y = 1 - X):
  species H2, HI, HII, e, HeI, HeII, HeIII from coupled Saha equilibria at (rho, T):
    H2  <-> 2 HI      K_d = (pi m_H kT/h^2)^{3/2} 16 / (Zrot_ns Zvib) exp(-chi_d/kT)
    HI  <-> HII + e   K_H = (2 pi m_e kT/h^2)^{3/2} exp(-chi_H/kT)
    HeI <-> HeII + e  K_1 = 4 (2 pi m_e kT/h^2)^{3/2} exp(-chi_He1/kT)
    HeII<-> HeIII + e K_2 = (2 pi m_e kT/h^2)^{3/2} exp(-chi_He2/kT)
  with charge neutrality n_e = nHII + nHeII + 2 nHeIII and nucleus conservation.
  Internal energy (zero point: cold ground state H2 and neutral He):
    e = 3/2 kT n_tot + nH2 e_rotvib(T) + chi_d (nHI + nHII)/2
        + chi_H nHII + chi_He1 (nHeII + nHeIII) + chi_He2 nHeIII
  Pressure P = n_tot k T. The H2 rotational energy assumes a frozen 3:1 ortho:para mix.
  The adiabatic sound speed is cs^2 = (dP/drho)_esp + (P/rho^2)(dP/desp)_rho, evaluated
  by finite differences on the (rho, esp) grid.

Output: a flat binary file (format in eos_table.hpp) with, on log10 spaced grids,
  (log10 rho, log10 esp) -> P, cs^2, log10 T   and   (log10 rho, log10 T) -> esp.
Everything is written in cgs; AthenaPK converts to code units when it reads the file.

Usage:
  gen_eos_table.py --out eos_table.bin                 # 180 x 220 x 200 default table
  gen_eos_table.py --out <file> --nr 400 --ne 1000 --nT 920 --esp-pad 399 --T-pad 450
                                                       # the high resolution table, see
                                                       # build_hires_v5_cold.sh
  gen_eos_table.py --check                             # print gamma_eff(T), no table
"""

import argparse
import os

import numpy as np

# ---- cgs constants ----
m_H = 1.6726219e-24
m_e = 9.1093837e-28
k_B = 1.380649e-16
h_pl = 6.62607015e-27
eV = 1.602176634e-12
chi_d = 4.4781 * eV  # H2 dissociation
chi_H = 13.598 * eV  # H ionization
chi_He1 = 24.587 * eV  # He  -> He+
chi_He2 = 54.418 * eV  # He+ -> He++
theta_rot = 85.4  # H2 rotational temperature [K]
theta_vib = 5987.0  # H2 vibrational temperature [K]

# ---- composition (mu = 2.33 molecular: 1/mu = X/2 + Y/4) ----
X = 0.716
Y = 1.0 - X

# ---- file format, must match eos_table.hpp ----
EOS_MAGIC = 0x454F535441424C31  # ASCII "EOSTABL1"
EOS_VERSION = 2
EOS_FLAG_CGS = 1  # bit 0: axes and arrays are stored in cgs

# Internal working scales. The tables are computed in these units and converted back to
# cgs when written, so they do not affect what AthenaPK reads beyond round-off. They are
# kept fixed so that a rebuild reproduces existing tables bit for bit.
rho0 = 5.467e-19  # g/cm^3
v0 = 1.9e4  # cm/s
e_unit = rho0 * v0 * v0  # erg/cm^3


def lam3(T, m):  # (2 pi m k T / h^2)^{3/2}
    return (2.0 * np.pi * m * k_B * T / (h_pl * h_pl)) ** 1.5


def Zrot_ns(T):
    # Nuclear spin weighted H2 rotational partition function for the dissociation Saha
    # equation: para (even J, weight 1) plus ortho (odd J, weight 3), energies measured
    # from the para J=0 level (consistent with chi_d referenced to v=0, J=0).
    J = np.arange(0, 60)
    w = np.where(J % 2 == 0, 1.0, 3.0)
    x = J * (J + 1) * theta_rot / np.maximum(T, 1.0)
    return np.sum(w * (2 * J + 1) * np.exp(-x))


def _erot_species(T, Jvals):
    # Rotational energy per molecule of one spin species [erg], measured from the species'
    # own ground state (para J=0, ortho J=1) so that e_rot -> 0 for T -> 0.
    E = (Jvals * (Jvals + 1) - Jvals[0] * (Jvals[0] + 1)) * k_B * theta_rot
    x = Jvals * (Jvals + 1) * theta_rot / np.maximum(T, 1.0)
    w = (2 * Jvals + 1) * np.exp(-x)
    Z = np.sum(w)
    return np.sum(w * E) / Z if Z > 0 else 0.0


def erot_per_H2(T):
    # Rotational energy per molecule [erg] for a frozen 3:1 ortho:para mixture (ortho/para
    # conversion is slow, so the gas is not in rotational spin equilibrium).
    Jeven = np.arange(0, 60, 2)
    Jodd = np.arange(1, 60, 2)
    return 0.25 * _erot_species(T, Jeven) + 0.75 * _erot_species(T, Jodd)


def Zvib(T):
    return 1.0 / (1.0 - np.exp(-theta_vib / np.maximum(T, 1.0)))


def evib_per_H2(T):
    # Harmonic vibrational energy per molecule [erg]. Below ~4 K the exponential overflows
    # to inf and the result is exactly 0, which is the correct limit, so the overflow
    # warning is silenced. (The equivalent exp(-x)/(1-exp(-x)) form differs in the last
    # bit and would change existing tables.)
    with np.errstate(over="ignore"):
        return k_B * theta_vib / (np.exp(theta_vib / np.maximum(T, 1.0)) - 1.0)


def solve_saha(rho, T):
    """Return the number densities [cm^-3] of all species and n_e at (rho [g/cm^3], T [K])."""
    nH_tot = rho * X / m_H  # H nuclei per cm^3
    nHe_tot = rho * Y / (4.0 * m_H)
    le = lam3(T, m_e)
    # H2 <-> 2H: n_H^2/n_H2 = (pi m_H kT/h^2)^{3/2} g_H^2 / (g_el,H2 Zrot_ns Zvib)
    # exp(-chi_d/kT) with g_H = 4 (electron spin x nuclear spin) and g_el,H2 = 1. This
    # reproduces the ~8% dissociation of H2 at 3000 K and 1 atm (JANAF).
    Kd = (
        (np.pi * m_H * k_B * T / (h_pl * h_pl)) ** 1.5
        * 16.0
        / (Zrot_ns(T) * Zvib(T))
        * np.exp(-chi_d / (k_B * T))
    )
    KH = le * np.exp(-chi_H / (k_B * T))
    K1 = 4.0 * le * np.exp(-chi_He1 / (k_B * T))
    K2 = 1.0 * le * np.exp(-chi_He2 / (k_B * T))

    def neutrality(ne):
        ne = max(ne, 1e-300)
        # H: 2 nHI^2/Kd + nHI (1 + KH/ne) - nH_tot = 0, solved with the stable root
        # nHI = 2 nH_tot / (b + sqrt(b^2 + 8 nH_tot/Kd)), so that Kd -> 0 (cold) gives
        # nHI -> 0 (all H2) instead of an overflow.
        rH = KH / ne
        b = 1.0 + rH
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            disc = np.sqrt(b * b + 8.0 * nH_tot / Kd) if Kd > 0 else np.inf
            nHI = 2.0 * nH_tot / (b + disc)
        if not np.isfinite(nHI):
            nHI = 0.0
        nHII = rH * nHI
        # He
        r1 = K1 / ne
        r2 = K2 / ne
        nHeI = nHe_tot / (1.0 + r1 + r1 * r2)
        nHeII = r1 * nHeI
        nHeIII = r1 * r2 * nHeI
        ne_new = nHII + nHeII + 2.0 * nHeIII
        return ne_new, (nHI, nHII, nHeI, nHeII, nHeIII)

    # bisection in log10(n_e) on the neutrality residual
    lo, hi = -40.0, np.log10(max(nH_tot + 2 * nHe_tot, 1e-30)) + 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        ne = 10.0**mid
        ne_new, _ = neutrality(ne)
        if ne_new > ne:
            lo = mid
        else:
            hi = mid
    ne = 10.0 ** (0.5 * (lo + hi))
    _, (nHI, nHII, nHeI, nHeII, nHeIII) = neutrality(ne)
    nH2 = max((nH_tot - nHI - nHII) / 2.0, 0.0)  # H nucleus conservation
    return dict(H2=nH2, HI=nHI, HII=nHII, e=ne, HeI=nHeI, HeII=nHeII, HeIII=nHeIII)


def P_e_of_rho_T(rho, T):
    """Pressure and internal energy density [erg/cm^3] at (rho [g/cm^3], T [K])."""
    s = solve_saha(rho, T)
    n_tot = s["H2"] + s["HI"] + s["HII"] + s["e"] + s["HeI"] + s["HeII"] + s["HeIII"]
    P = n_tot * k_B * T
    e = 1.5 * k_B * T * n_tot
    e += s["H2"] * (erot_per_H2(T) + evib_per_H2(T))
    e += chi_d * (s["HI"] + s["HII"]) / 2.0
    e += chi_H * s["HII"]
    e += chi_He1 * (s["HeII"] + s["HeIII"]) + chi_He2 * s["HeIII"]
    return P, e


def build_table(
    path,
    nr=180,
    ne=220,
    nT=200,
    rho_phys_min=1e-20,
    rho_phys_max=1e0,
    T_min=8.0,
    T_max=3.0e5,
    esp_pad_pts=0,
    T_pad_pts=0,
    T_pad_hi_pts=0,
):
    """Compute the table and write it to `path` (binary format of eos_table.hpp).

    The esp axis spans the internal energies of the (rho, T) grid between T_min and T_max.
    The pad arguments extend the T axis (T_pad_pts below T_min, T_pad_hi_pts above T_max)
    and the esp axis (esp_pad_pts below) at the SAME spacing, so the nodes of an unpadded
    table keep their exact values. Changing T_min, T_max or nT instead would move every
    node.
    """
    v0sq = v0 * v0
    rho_code = np.logspace(
        np.log10(rho_phys_min / rho0), np.log10(rho_phys_max / rho0), nr
    )
    rho_cgs = rho_code * rho0

    # T axis
    lT_base = np.log10(T_min)
    dlT = (np.log10(T_max) - lT_base) / (nT - 1)
    lT = np.linspace(lT_base, np.log10(T_max), nT)
    if T_pad_pts > 0:
        lT = np.concatenate([lT_base - dlT * np.arange(T_pad_pts, 0, -1), lT])
    if T_pad_hi_pts > 0:
        lT = np.concatenate(
            [lT, np.log10(T_max) + dlT * np.arange(1, T_pad_hi_pts + 1)]
        )
    Tgrid = 10.0**lT
    nT_tot = lT.size

    # (rho, T) grid: esp(rho, T) and P(rho, T), in the internal units
    esp_rhoT = np.zeros((nr, nT_tot))
    P_rhoT = np.zeros((nr, nT_tot))
    for i, rc in enumerate(rho_cgs):
        for j, T in enumerate(Tgrid):
            P, ed = P_e_of_rho_T(rc, T)
            esp_rhoT[i, j] = (ed / rc) / v0sq
            P_rhoT[i, j] = P / e_unit
    log10T = np.log10(Tgrid)

    # esp axis. Its range and spacing come from the unpadded T columns only, extra nodes
    # are added at the same spacing.
    hi0 = esp_rhoT.shape[1] - T_pad_hi_pts if T_pad_hi_pts > 0 else esp_rhoT.shape[1]
    base = esp_rhoT[:, T_pad_pts:hi0]
    esp_min = base.min() * 0.999
    esp_max = base.max() * 1.001
    le_base = np.linspace(np.log10(esp_min), np.log10(esp_max), ne)
    dle_base = le_base[1] - le_base[0]
    le_all = le_base
    if esp_pad_pts > 0:
        le_all = np.concatenate(
            [le_base[0] - dle_base * np.arange(esp_pad_pts, 0, -1), le_all]
        )
    if T_pad_hi_pts > 0:
        esp_hi = esp_rhoT.max() * 1.001
        n_up = int(np.ceil((np.log10(esp_hi) - le_base[-1]) / dle_base))
        if n_up > 0:
            le_all = np.concatenate(
                [le_all, le_base[-1] + dle_base * np.arange(1, n_up + 1)]
            )
    esp_code = 10.0**le_all
    ne_tot = le_all.size
    # The padded esp nodes must lie inside the (rho, T) coverage, otherwise they would just
    # repeat the T_min edge value.
    if esp_pad_pts > 0 and esp_code[0] < esp_rhoT.min():
        raise SystemExit(
            "T_pad_pts=%d is not enough for esp_pad_pts=%d: the lowest esp node %.4e erg/g "
            "is below the coldest tabulated esp %.4e erg/g. Increase T_pad_pts."
            % (T_pad_pts, esp_pad_pts, esp_code[0] * v0sq, esp_rhoT.min() * v0sq)
        )
    # The (rho, esp) grid is rectangular, but the esp reached at T_max depends on rho (it
    # is largest at the lowest density). Rows at higher density therefore have esp nodes
    # above their own esp(T_max), where T and P are flat in esp. Report how wide that is.
    nsat = np.zeros(nr, dtype=int)
    for i in range(nr):
        nsat[i] = int(np.count_nonzero(esp_code > esp_rhoT[i].max()))
    bad = np.where(nsat >= 2)[0]
    if len(bad):
        T_top = 10.0 ** lT[-1]
        print(
            "NOTE: %d of %d density rows (rho = %.3e..%.3e g/cm^3) have >= 2 esp nodes above "
            "their own esp(T = %.4g K), widest flat run %d of %d esp nodes. P is flat in esp "
            "there. Raise T_max (--T-pad-hi) if a run approaches %.4g K."
            % (
                len(bad),
                nr,
                rho_cgs[bad[0]],
                rho_cgs[bad[-1]],
                T_top,
                int(nsat.max()),
                ne_tot,
                T_top,
            )
        )

    # (rho, esp) grid: invert esp(T) at each rho (monotone), then P(T)
    P_re = np.zeros((nr, ne_tot))
    logT_re = np.zeros((nr, ne_tot))
    for i in range(nr):
        Tof = np.interp(
            np.log10(esp_code),
            np.log10(esp_rhoT[i]),
            log10T,
            left=log10T[0],
            right=log10T[-1],
        )
        logT_re[i] = Tof
        P_re[i] = np.interp(Tof, log10T, P_rhoT[i])
    # cs^2 = (dP/drho)_esp + (P/rho^2)(dP/desp)_rho, esp being the SPECIFIC internal energy.
    # Ideal gas check: P = (g-1) rho esp gives cs^2 = g P / rho.
    lr = np.log10(rho_code)
    le = np.log10(esp_code)
    dP_dlr = np.gradient(P_re, lr, axis=0)
    dP_dle = np.gradient(P_re, le, axis=1)
    RHO = rho_code[:, None]
    ESP = esp_code[None, :]
    dP_drho_e = dP_dlr / (np.log(10) * RHO)
    dP_desp = dP_dle / (np.log(10) * ESP)
    P_over_rho2 = P_re / (RHO * RHO)
    cs2_re = dP_drho_e + P_over_rho2 * dP_desp
    cs2_re = np.maximum(cs2_re, 1e-12)

    # Write in cgs, header first (see eos_table.hpp). The specific energy unit is computed
    # as v0^2 here (not e_unit / rho0, which differs in the last bit).
    esp_unit = v0 * v0
    with open(path, "wb") as fb:
        np.array(
            [EOS_MAGIC, EOS_VERSION, EOS_FLAG_CGS, nr, ne_tot, nT_tot], dtype=np.int64
        ).tofile(fb)
        np.array(
            [
                lr[0] + np.log10(rho0),
                lr[1] - lr[0],
                le[0] + np.log10(esp_unit),
                le[1] - le[0],
                log10T[0],
                log10T[1] - log10T[0],
                0.0,
                0.0,
            ],
            dtype=np.float64,
        ).tofile(fb)
        (P_re * e_unit).astype(np.float64).tofile(fb)
        (cs2_re * esp_unit).astype(np.float64).tofile(fb)
        logT_re.astype(np.float64).tofile(fb)
        (esp_rhoT * esp_unit).astype(np.float64).tofile(fb)
    print("wrote %s (%d bytes)" % (path, os.path.getsize(path)))
    print(
        "  (rho, esp) grid %d x %d, (rho, T) grid %d x %d (esp_pad=%d, T_pad=%d, T_pad_hi=%d)"
        % (nr, ne_tot, nr, nT_tot, esp_pad_pts, T_pad_pts, T_pad_hi_pts)
    )
    print(
        "  log10 rho in [%.4f, %.4f] g/cm^3, log10 esp in [%.4f, %.4f] erg/g, "
        "log10 T in [%.4f, %.4f] K"
        % (
            lr[0] + np.log10(rho0),
            lr[-1] + np.log10(rho0),
            le[0] + np.log10(esp_unit),
            le[-1] + np.log10(esp_unit),
            log10T[0],
            log10T[-1],
        )
    )
    print(
        "  cs2 range [%.3e, %.3e] cm^2/s^2"
        % (cs2_re.min() * esp_unit, cs2_re.max() * esp_unit)
    )
    return path


def check():
    """Print gamma_eff = 1 + P/e at a few densities: ~7/5 for warm H2, dips at H2
    dissociation (~2000 K) and H ionization (~1e4 K), ~5/3 for atomic or ionized gas."""
    for rho in [1e-10, 1e-6, 1e-2]:
        print("\n# rho=%.0e g/cm^3 : T[K]  x_H2  x_HII  gamma_eff" % rho)
        for T in [50, 100, 300, 1000, 2000, 3000, 6000, 1e4, 2e4, 5e4, 1e5]:
            P, e = P_e_of_rho_T(rho, T)
            s = solve_saha(rho, T)
            nHt = rho * X / m_H
            print(
                "  %8.0f  %.3f  %.3e  %.4f"
                % (T, 2 * s["H2"] / nHt, s["HII"] / nHt, 1.0 + P / e)
            )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", help="output table file")
    ap.add_argument("--nr", type=int, default=180, help="number of density nodes")
    ap.add_argument(
        "--ne", type=int, default=220, help="number of esp nodes (before padding)"
    )
    ap.add_argument(
        "--nT", type=int, default=200, help="number of T nodes (before padding)"
    )
    ap.add_argument(
        "--esp-pad", type=int, default=0, help="esp nodes added below the axis"
    )
    ap.add_argument(
        "--T-pad", type=int, default=0, help="T nodes added below T_min = 8 K"
    )
    ap.add_argument(
        "--T-max", type=float, default=3.0e5, help="upper end of the T axis [K]"
    )
    ap.add_argument("--T-pad-hi", type=int, default=0, help="T nodes added above T_max")
    ap.add_argument("--check", action="store_true", help="print gamma_eff(T) and exit")
    args = ap.parse_args()
    if args.check:
        check()
    elif args.out:
        build_table(
            args.out,
            nr=args.nr,
            ne=args.ne,
            nT=args.nT,
            T_max=args.T_max,
            esp_pad_pts=args.esp_pad,
            T_pad_pts=args.T_pad,
            T_pad_hi_pts=args.T_pad_hi,
        )
    else:
        ap.error("either --out or --check is required")
