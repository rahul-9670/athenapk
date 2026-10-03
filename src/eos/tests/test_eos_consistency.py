#!/usr/bin/env python3
# ========================================================================================
# AthenaPK - a performance portable block structured AMR astrophysical MHD code.
# Copyright (c) 2026, Athena-Parthenon Collaboration. All rights reserved.
# Licensed under the BSD 3-Clause License (the "LICENSE").
# ========================================================================================
"""Consistency checks for the tabulated hydrogen EOS (src/eos/gen_eos_table.py).

  1. Charge neutrality of the Saha solution (to round-off).
  2. H and He nucleus conservation (to round-off).
  3. Thermodynamic consistency of P and e from the same Saha solution:
     -rho^2 (de/drho)_T = T (dP/dT)_rho - P   (e = specific internal energy).
  4. Table fidelity: P(rho, esp) from the table reproduces the Saha P(rho, T).
  5. Sound speed: cs^2 from the table matches the isentropic derivative of the Saha EOS.

Checks 1 to 3 only need the generator. Checks 4 and 5 need a table file, given as the
first argument or via $EOS_TABLE_BIN, and are skipped without one. 4 and 5 gate the
median error; the worst case sits at the H2 dissociation and H ionization kinks, where
bilinear interpolation is least accurate, and is printed for information.

usage: test_eos_consistency.py [table.bin]
Exit code 0 if all checks pass.
"""

import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.normpath(os.path.join(HERE, "..", "gen_eos_table.py"))
TABLE = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("EOS_TABLE_BIN", "")

sys.dont_write_bytecode = True  # no __pycache__ next to the generator
spec = importlib.util.spec_from_file_location("gen_eos_table", GEN)
G = importlib.util.module_from_spec(spec)
spec.loader.exec_module(G)


def load_table(path):
    """Read a table file (format in eos_table.hpp). Everything is in cgs."""
    with open(path, "rb") as f:
        iv = np.fromfile(f, dtype=np.int64, count=6)
        assert (
            iv[0] == G.EOS_MAGIC and iv[1] == G.EOS_VERSION and (iv[2] & G.EOS_FLAG_CGS)
        )
        nr, ne, nT = (int(x) for x in iv[3:6])
        lr0, dlr, le0, dle, lT0, dlT, _, _ = np.fromfile(f, dtype=np.float64, count=8)
        P = np.fromfile(f, dtype=np.float64, count=nr * ne).reshape(nr, ne)
        cs2 = np.fromfile(f, dtype=np.float64, count=nr * ne).reshape(nr, ne)
        logT = np.fromfile(f, dtype=np.float64, count=nr * ne).reshape(nr, ne)
        espT = np.fromfile(f, dtype=np.float64, count=nr * nT).reshape(nr, nT)
        assert f.read() == b"", "trailing bytes in table file"
    return dict(
        nr=nr,
        ne=ne,
        nT=nT,
        lr0=lr0,
        dlr=dlr,
        le0=le0,
        dle=dle,
        lT0=lT0,
        dlT=dlT,
        P=P,
        cs2=cs2,
        logT=logT,
        espT=espT,
    )


def bilin(A, x0, dx, y0, dy, x, y):
    """Same clamped bilinear interpolation as EosTable::Bilin."""
    n1, n2 = A.shape
    fi = (x - x0) / dx
    i = min(max(int(fi), 0), n1 - 2)
    ti = min(max(fi - i, 0.0), 1.0)
    fj = (y - y0) / dy
    j = min(max(int(fj), 0), n2 - 2)
    tj = min(max(fj - j, 0.0), 1.0)
    return (
        (1 - ti) * (1 - tj) * A[i, j]
        + ti * (1 - tj) * A[i + 1, j]
        + (1 - ti) * tj * A[i, j + 1]
        + ti * tj * A[i + 1, j + 1]
    )


def in_table(t, lr, le):
    return (
        t["lr0"] <= lr <= t["lr0"] + (t["nr"] - 1) * t["dlr"]
        and t["le0"] <= le <= t["le0"] + (t["ne"] - 1) * t["dle"]
    )


def main():
    fails = 0
    # protostellar collapse regime, from the envelope to the first core, crossing H2
    # dissociation and H ionization
    rhos = np.logspace(-16, -4, 13)  # g/cm^3
    Ts = np.logspace(1.2, 4.3, 16)  # ~16 K .. 20000 K

    print("--- 1/2: charge neutrality and nucleus conservation ---")
    worst_neut = worst_H = worst_He = 0.0
    for rho in rhos:
        nH_tot = rho * G.X / G.m_H
        nHe_tot = rho * G.Y / (4.0 * G.m_H)
        for T in Ts:
            s = G.solve_saha(rho, T)
            net = s["HII"] + s["HeII"] + 2.0 * s["HeIII"] - s["e"]
            # normalized by n_H (n_e itself goes to zero in cold gas)
            worst_neut = max(worst_neut, abs(net) / nH_tot)
            H = 2 * s["H2"] + s["HI"] + s["HII"]
            He = s["HeI"] + s["HeII"] + s["HeIII"]
            worst_H = max(worst_H, abs(H - nH_tot) / nH_tot)
            worst_He = max(worst_He, abs(He - nHe_tot) / nHe_tot)
    ok = worst_neut < 1e-10 and worst_H < 1e-10 and worst_He < 1e-10
    print(
        f"  |sum Z n|/n_H = {worst_neut:.2e}, H {worst_H:.2e}, He {worst_He:.2e}  "
        f"{'PASS' if ok else 'FAIL'}"
    )
    fails += not ok

    print("--- 3: thermodynamic identity -rho^2 (de/drho)_T = T (dP/dT)_rho - P ---")
    worst = 0.0
    h = 1e-3  # relative finite difference step
    n = 0
    for rho in rhos[1:-1]:
        for T in Ts[1:-1]:
            _, ep = G.P_e_of_rho_T(rho * (1 + h), T)
            _, em = G.P_e_of_rho_T(rho * (1 - h), T)
            dedrho = (ep / (rho * (1 + h)) - em / (rho * (1 - h))) / (2 * h * rho)
            lhs = -rho * rho * dedrho
            PTp, _ = G.P_e_of_rho_T(rho, T * (1 + h))
            PTm, _ = G.P_e_of_rho_T(rho, T * (1 - h))
            P0, _ = G.P_e_of_rho_T(rho, T)
            rhs = T * (PTp - PTm) / (2 * h * T) - P0
            worst = max(worst, abs(lhs - rhs) / max(abs(lhs), abs(rhs), P0))
            n += 1
    ok = worst < 3e-3  # limited by the finite differences
    print(
        f"  worst relative residual over {n} states = {worst:.2e} (tol 3e-3)  "
        f"{'PASS' if ok else 'FAIL'}"
    )
    fails += not ok

    if not TABLE:
        print("--- 4/5: SKIPPED (no table file given) ---")
    else:
        t = load_table(TABLE)
        print(f"--- 4: table P(rho, esp) vs Saha P(rho, T)   [{TABLE}] ---")
        errs = []
        wloc = (0.0, 0.0, 0.0)
        for rho in np.logspace(-16, -4, 25):
            lr = np.log10(rho)
            for T in np.logspace(1.2, 4.3, 30):
                P_saha, e_saha = G.P_e_of_rho_T(rho, T)
                le = np.log10(e_saha / rho)
                if not in_table(t, lr, le):
                    continue
                P_tab = bilin(t["P"], t["lr0"], t["dlr"], t["le0"], t["dle"], lr, le)
                err = abs(P_tab - P_saha) / P_saha
                errs.append(err)
                if err > wloc[0]:
                    wloc = (err, T, rho)
        e = np.array(errs)
        ok = len(e) > 0 and np.median(e) < 1e-2
        print(
            f"  n={len(e)}  median={np.median(e):.2e}  p99={np.percentile(e, 99):.2e}  "
            f"max={e.max():.2e} at T={wloc[1]:.0f} K, rho={wloc[2]:.0e} g/cm^3  "
            f"(gate: median < 1e-2)  {'PASS' if ok else 'FAIL'}"
        )
        fails += not ok

        print("--- 5: table cs^2 vs isentropic cs^2 of the Saha EOS ---")
        errs = []
        wloc = (0.0, 0.0, 0.0)

        def PU(r, T):
            P, e = G.P_e_of_rho_T(r, T)
            return P, e / r

        for rho in np.logspace(-15, -5, 21):
            lr = np.log10(rho)
            for T in np.logspace(1.4, 4.2, 24):
                P0, u0 = PU(rho, T)
                le = np.log10(u0)
                if not in_table(t, lr, le):
                    continue
                Pp, up = PU(rho * (1 + h), T)
                Pm, um = PU(rho * (1 - h), T)
                dP_drho_T = (Pp - Pm) / (2 * h * rho)
                du_drho_T = (up - um) / (2 * h * rho)
                PTp, uTp = PU(rho, T * (1 + h))
                PTm, uTm = PU(rho, T * (1 - h))
                dP_dT_rho = (PTp - PTm) / (2 * h * T)
                du_dT_rho = (uTp - uTm) / (2 * h * T)
                if du_dT_rho == 0.0:
                    continue
                # change variables from (rho, T) to (rho, u)
                dP_drho_u = dP_drho_T - dP_dT_rho * du_drho_T / du_dT_rho
                dP_du_rho = dP_dT_rho / du_dT_rho
                cs2_saha = dP_drho_u + (P0 / (rho * rho)) * dP_du_rho
                cs2_tab = bilin(
                    t["cs2"], t["lr0"], t["dlr"], t["le0"], t["dle"], lr, le
                )
                if cs2_saha <= 0 or cs2_tab <= 0:
                    continue
                err = abs(cs2_tab - cs2_saha) / cs2_saha
                errs.append(err)
                if err > wloc[0]:
                    wloc = (err, T, rho)
        e = np.array(errs)
        ok = len(e) > 0 and np.median(e) < 3e-2
        print(
            f"  n={len(e)}  median={np.median(e):.2e}  max={e.max():.2e} at "
            f"T={wloc[1]:.0f} K, rho={wloc[2]:.0e} g/cm^3  (gate: median < 3e-2)  "
            f"{'PASS' if ok else 'FAIL'}"
        )
        fails += not ok

    print(f"\n{'ALL CHECKS PASS' if fails == 0 else f'{fails} CHECK(S) FAILED'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
