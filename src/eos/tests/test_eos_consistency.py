#!/usr/bin/env python
"""Flagship Phase 3 EOS consistency gate for the tabulated multi-Saha protostellar EOS.

Read-only characterization of the EXISTING EOS (src/eos/gen_eos_table.py + the shipped
src/eos/eos_table.bin). No numerics are changed; this quantifies whether the EOS meets the
Phase-3 acceptance criteria before any Helmholtz-free-energy rework is considered.

Gates:
  1. Charge neutrality of the Saha solution to round-off (Phase-3 gate).
  2. Element (H, He nucleus) conservation to round-off (Phase-3 gate).
  3. Thermodynamic consistency: the fundamental identity (du/dv)|_T = T (dP/dT)|_rho - P
     (u = specific internal energy, v = 1/rho). ANY EOS derived from one free energy must
     satisfy it; here P and e are produced by the same Saha solver, so this checks that they
     are mutually thermodynamically consistent (a Maxwell-relation test, entropy-free).
  4. Shipped-binary fidelity: the eos_table.bin P(rho,esp) reproduces the Saha P(rho,T) to
     the table's interpolation accuracy (catches a stale/corrupt table or a unit-scale bug).

Run: /beegfs/u/bbg6470/venvs/analysis_env/bin/python src/eos/tests/test_eos_consistency.py
"""
import os
import sys
import importlib.util

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.normpath(os.path.join(HERE, "..", "gen_eos_table.py"))
# Table under test: argv[1] or $EOS_TABLE_BIN overrides the shipped table, so a candidate
# hi-res table can be validated against the same gates before any production swap.
TABLE = (sys.argv[1] if len(sys.argv) > 1 else
         os.environ.get("EOS_TABLE_BIN",
                        os.path.normpath(os.path.join(HERE, "..", "eos_table.bin"))))

# import gen_eos_table as a module (its build action is guarded by __main__)
spec = importlib.util.spec_from_file_location("gen_eos_table", GEN)
G = importlib.util.module_from_spec(spec)
spec.loader.exec_module(G)


def load_table(path):
    """Read the flat binary written by EosTable::Load (see eos_table.hpp)."""
    with open(path, "rb") as f:
        hdr = np.fromfile(f, dtype=np.int64, count=3)
        nr, ne, nT = (int(x) for x in hdr)
        g = np.fromfile(f, dtype=np.float64, count=6)
        lr0, dlr, le0, dle, lT0, dlT = g
        P = np.fromfile(f, dtype=np.float64, count=nr * ne).reshape(nr, ne)
        cs2 = np.fromfile(f, dtype=np.float64, count=nr * ne).reshape(nr, ne)
        logT = np.fromfile(f, dtype=np.float64, count=nr * ne).reshape(nr, ne)
        espT = np.fromfile(f, dtype=np.float64, count=nr * nT).reshape(nr, nT)
    return dict(nr=nr, ne=ne, nT=nT, lr0=lr0, dlr=dlr, le0=le0, dle=dle, lT0=lT0, dlT=dlT,
                P=P, cs2=cs2, logT=logT, espT=espT)


def bilin(A, x0, dx, y0, dy, x, y):
    n1, n2 = A.shape
    fi = (x - x0) / dx
    i = min(max(int(np.floor(fi)), 0), n1 - 2)
    ti = min(max(fi - i, 0.0), 1.0)
    fj = (y - y0) / dy
    j = min(max(int(np.floor(fj)), 0), n2 - 2)
    tj = min(max(fj - j, 0.0), 1.0)
    return ((1 - ti) * (1 - tj) * A[i, j] + ti * (1 - tj) * A[i + 1, j] +
            (1 - ti) * tj * A[i, j + 1] + ti * tj * A[i + 1, j + 1])


def main():
    fails = 0
    # sample grid over the physically relevant FHC collapse regime.
    rhos = np.logspace(-16, -4, 13)      # g/cm^3 (envelope -> first core)
    Ts = np.logspace(1.2, 4.3, 16)       # ~16 K .. 20000 K (crosses H2 diss + H ioniz)

    # ---- Gate 1 & 2: charge neutrality + element conservation ----
    print("--- GATE 1/2: charge neutrality + element conservation (Saha) ---")
    worst_neut = 0.0
    worst_H = 0.0
    worst_He = 0.0
    for rho in rhos:
        nH_tot = rho * G.X / G.m_H
        nHe_tot = rho * G.Y / (4.0 * G.m_H)
        for T in Ts:
            s = G.solve_saha(rho, T)
            net = s['HII'] + s['HeII'] + 2.0 * s['HeIII'] - s['e']  # sum Z_j n_j (should be 0)
            # normalize by n_H_tot, NOT n_e: at cold T ionization -> 0 and n_e is pinned at its
            # 1e-40 bisection floor, so |net|/n_e is ill-defined there. n_H_tot is the physical
            # charge scale and stays finite.
            worst_neut = max(worst_neut, abs(net) / nH_tot)
            H = 2 * s['H2'] + s['HI'] + s['HII']
            He = s['HeI'] + s['HeII'] + s['HeIII']
            worst_H = max(worst_H, abs(H - nH_tot) / nH_tot)
            worst_He = max(worst_He, abs(He - nHe_tot) / nHe_tot)
    g12 = worst_neut < 1e-10 and worst_H < 1e-10 and worst_He < 1e-10
    print(f"  worst |sum Z n|/n_H = {worst_neut:.2e} ; H-cons = {worst_H:.2e} ; "
          f"He-cons = {worst_He:.2e}  {'PASS' if g12 else 'FAIL'}")
    fails += not g12

    # ---- Gate 3: thermodynamic consistency (du/dv)|_T = T (dP/dT)|_rho - P ----
    # u = specific internal energy = e_density/rho ; v = 1/rho.  (du/dv)|_T = -rho^2 (du/drho)|_T.
    print("--- GATE 3: thermodynamic identity  -rho^2 (du/drho)_T == T (dP/dT)_rho - P ---")
    worst_thermo = 0.0
    hln = 1e-3  # log-space finite-difference step
    checked = 0
    for rho in rhos[1:-1]:
        for T in Ts[1:-1]:
            # central differences in log rho (fixed T) and log T (fixed rho)
            Pp, ep = G.P_e_of_rho_T(rho * (1 + hln), T)
            Pm, em = G.P_e_of_rho_T(rho * (1 - hln), T)
            up = ep / (rho * (1 + hln)); um = em / (rho * (1 - hln))
            dudrho_T = (up - um) / (2 * hln * rho)
            lhs = -rho * rho * dudrho_T
            PTp, _ = G.P_e_of_rho_T(rho, T * (1 + hln))
            PTm, _ = G.P_e_of_rho_T(rho, T * (1 - hln))
            dPdT_rho = (PTp - PTm) / (2 * hln * T)
            P0, _ = G.P_e_of_rho_T(rho, T)
            rhs = T * dPdT_rho - P0
            scale = max(abs(lhs), abs(rhs), P0)
            worst_thermo = max(worst_thermo, abs(lhs - rhs) / scale)
            checked += 1
    # tolerance: dominated by O(hln^2) truncation + composition-derivative curvature.
    g3 = worst_thermo < 3e-3
    print(f"  worst rel residual over {checked} states = {worst_thermo:.2e}  "
          f"(tol 3e-3, finite-diff limited)  {'PASS' if g3 else 'FAIL'}")
    fails += not g3

    # ---- Gate 4: shipped-binary fidelity vs the Saha EOS ----
    # Two parts: (a) TYPICAL fidelity (median) must be sub-percent; (b) the worst-case is a
    # CHARACTERIZATION -- bilinear interp on the finite log grid degrades to ~5-7% right at the
    # H2-dissociation (~2200 K) and H-ionization (~7000-10000 K) gamma-softening kinks. That
    # degradation is the actionable Phase-3 EOS item (finer/adaptive grid at the transitions),
    # NOT a bug: the underlying Saha EOS is thermodynamically consistent (gates 1-3).
    print("--- GATE 4: eos_table.bin P(rho,esp) reproduces Saha P(rho,T) ---")
    if not os.path.exists(TABLE):
        print(f"  SKIP: {TABLE} not found")
    else:
        t = load_table(TABLE)
        rho0 = G.rho0; v0 = G.v0
        e_unit = rho0 * v0 * v0
        errs = []
        worst_loc = (0.0, 0.0, 0.0)
        for rho in np.logspace(-16, -4, 25):
            lr = np.log10(rho / rho0)
            if lr < t['lr0'] or lr > t['lr0'] + (t['nr'] - 1) * t['dlr']:
                continue
            for T in np.logspace(1.2, 4.3, 30):
                Psaha, esaha = G.P_e_of_rho_T(rho, T)     # erg/cm^3
                le = np.log10((esaha / rho) / (v0 * v0))    # specific -> code v0^2
                if le < t['le0'] or le > t['le0'] + (t['ne'] - 1) * t['dle']:
                    continue
                P_tab_cgs = bilin(t['P'], t['lr0'], t['dlr'], t['le0'], t['dle'], lr, le) * e_unit
                re = abs(P_tab_cgs - Psaha) / Psaha
                errs.append(re)
                if re > worst_loc[0]:
                    worst_loc = (re, T, rho)
        e = np.array(errs)
        med = float(np.median(e)); p99 = float(np.percentile(e, 99)); mx = float(e.max())
        g4 = med < 1e-2  # gate on TYPICAL fidelity; worst-at-kinks is characterized, not gated
        print(f"  n={len(e)}  median={med:.2e}  p99={p99:.2e}  max={mx:.2e}  "
              f"(gate: median<1e-2)  {'PASS' if g4 else 'FAIL'}")
        print(f"  CHARACTERIZATION: worst {worst_loc[0]:.1%} at T={worst_loc[1]:.0f} K "
              f"rho={worst_loc[2]:.0e} (a dissociation/ionization kink) -> Phase-3 EOS item.")
        fails += not g4

        # ---- Gate 5: sound-speed (entropy-along-adiabat) consistency ----
        # The tabulated cs2 IS the adiabatic sound speed = the isentrope slope (dP/drho)|_s,
        # equivalently cs2 = (dP/drho)|_u + (P/rho^2)(dP/du)|_rho with u the SPECIFIC internal
        # energy. Recompute it INDEPENDENTLY from the Saha EOS by pointwise finite differences
        # of P_e_of_rho_T in (rho,T) + the chain rule to fixed-u partials, and compare to the
        # table's cs2 (which the generator built by grid-gradient of P over the (rho,esp) grid,
        # then the C++ reads bilinearly). Agreement validates the shipped sound speed -- the
        # entropy-along-adiabat consistency the Riemann solver relies on -- independently of P.
        print("--- GATE 5: sound-speed consistency  cs2_table == cs2_Saha(isentrope) ---")
        hln = 1e-3
        cs_errs = []
        cs_worst = (0.0, 0.0, 0.0)

        def _PU(r, tK):
            P, e = G.P_e_of_rho_T(r, tK)
            return P, e / r  # P[cgs], specific u[cgs]

        for rho in np.logspace(-15, -5, 21):
            lr = np.log10(rho / rho0)
            if lr < t['lr0'] or lr > t['lr0'] + (t['nr'] - 1) * t['dlr']:
                continue
            for T in np.logspace(1.4, 4.2, 24):
                P0, u0 = _PU(rho, T)
                Pp, up = _PU(rho * (1 + hln), T); Pm, um = _PU(rho * (1 - hln), T)
                dP_drho_T = (Pp - Pm) / (2 * hln * rho)
                du_drho_T = (up - um) / (2 * hln * rho)
                PTp, uTp = _PU(rho, T * (1 + hln)); PTm, uTm = _PU(rho, T * (1 - hln))
                dP_dT_rho = (PTp - PTm) / (2 * hln * T)
                du_dT_rho = (uTp - uTm) / (2 * hln * T)
                if du_dT_rho == 0.0:
                    continue
                dT_drho_u = -du_drho_T / du_dT_rho
                dP_drho_u = dP_drho_T + dP_dT_rho * dT_drho_u
                dP_du_rho = dP_dT_rho / du_dT_rho
                cs2_saha_code = (dP_drho_u + (P0 / (rho * rho)) * dP_du_rho) / (v0 * v0)
                le = np.log10((u0) / (v0 * v0))
                if le < t['le0'] or le > t['le0'] + (t['ne'] - 1) * t['dle']:
                    continue
                cs2_tab = bilin(t['cs2'], t['lr0'], t['dlr'], t['le0'], t['dle'], lr, le)
                if cs2_saha_code <= 0 or cs2_tab <= 0:
                    continue
                re = abs(cs2_tab - cs2_saha_code) / cs2_saha_code
                cs_errs.append(re)
                if re > cs_worst[0]:
                    cs_worst = (re, T, rho)
        ce = np.array(cs_errs)
        cmed = float(np.median(ce)); cmx = float(ce.max())
        # cs2 carries a derivative of P so it is noisier (FD + grid-gradient); gate the median.
        g5 = cmed < 3e-2
        print(f"  n={len(ce)}  median={cmed:.2e}  max={cmx:.2e}  (gate: median<3e-2)  "
              f"{'PASS' if g5 else 'FAIL'}")
        print(f"  CHARACTERIZATION: worst {cs_worst[0]:.1%} at T={cs_worst[1]:.0f} K "
              f"rho={cs_worst[2]:.0e} (cs2 is derivative-sensitive at the kinks).")
        fails += not g5

    print(f"\n{'ALL GATES PASS' if fails == 0 else 'GATE FAILURES'} ({fails} failures)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
