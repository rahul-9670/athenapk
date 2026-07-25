#!/usr/bin/env python
"""Flagship Phase 3 conductivity gate: independent check of the Wardle tensor.

Read-only. An INDEPENDENT pure-Python reimplementation of the Pandey & Wardle (2008)
conductivity tensor and the Ohmic/Hall/ambipolar diffusivities that
src/hydro/diffusion/ionization.hpp::Diffusivities computes on device. Gates 1-5 cover the
gas-phase balance (electrons + one ion) and the tensor ASSEMBLY; gates 6-8 cover the
GRAIN-INCLUSIVE model (MRN bins + grain charging), closing the Phase-3 conductivity item.

The grain reimplementation is ALGORITHMICALLY INDEPENDENT of the C++, not a transliteration:
observing that the per-bin charge balance n_e v_e e^psi = n_i v_i (1-psi) has NO bin dependence
(the pi a^2 capture cross-sections cancel between electrons and ions), psi is a single global
unknown and Z_k = psi tau_k. The whole coupled system then collapses to a 1-D bracketed root
find in r = n_e/n_i (brentq), versus the C++ relaxed fixed point on n_i with an inner Newton
on psi. Same physics, different algorithm -- so agreement is a real cross-check.

Gates:
  1. sigma_O identity: the tensor sigma_O = (ec/B) sum_j n_j Z_j beta_j must equal the scalar
     conductivity sum_j n_j Z_j^2 e^2/(m_j nu_j) EXACTLY -- the B cancels. Validates the
     beta/ecB assembly.
  2. eta_O identity: eta_O = c^2/(4 pi sigma_O).
  3. parity in B: eta_O, eta_A are EVEN in B; eta_H is ODD (sign follows B). This is the Hall
     sign structure the whistler test exercises.
  4. tensor algebra: sigma_perp^2 = sigma_H^2 + sigma_P^2 (definitional, guards the assembly).
  5. canonical curve: over a collapse rho-track, eta_A dominates at low density and Ohmic rises
     to dominate at high density (the standard non-ideal-MHD crossover; Wardle 2007 / Marchand+
     2016). Qualitative ordering gate.
  6. GRAINS -- charge neutrality n_i = n_e + sum_k (-Z_k) ng_k and the ionization balance
     P = alpha n_e n_i + Ki n_i hold on the independent solution (self-consistency of the
     Python solve, to roundoff).
  7. GRAINS -- physics structure: grains are negatively charged (psi<0), |Z_k| grows linearly
     with bin radius a_k (Z = psi tau, tau ~ a), grains SUPPRESS the free-electron density
     vs the gas-phase-only balance, and they vanish above T_subl.
  8. GRAINS -- direct C++-vs-Python numerical cross-check of the charge state (n_e, n_i, Z_k)
     AND eta_{O,H,A}, reading the xcheck_conductivity dump (f_dg=0.01). This is the
     grain-inclusive counterpart of the gas-phase cross-check, and the item Phase 3 listed
     as remaining.

Run (gates 1-7):
  /beegfs/u/bbg6470/venvs/analysis_env/bin/python src/hydro/diffusion/tests/test_conductivity.py
Run incl. gate 8 (needs the C++ dump):
  g++ -O2 -std=c++17 -Ihost_shims xcheck_conductivity.cpp -o <out>/xcheck
  <out>/xcheck 0.01 > <out>/xcheck_grains.csv
  ... test_conductivity.py --xcheck-csv <out>/xcheck_grains.csv
"""
import sys
import numpy as np
from scipy.optimize import brentq

# CGS constants -- MUST match ionization.hpp namespace cgs.
e_chg = 4.80320425e-10
m_e = 9.10938370e-28
m_H = 1.67262192e-24
c_light = 2.99792458e10
k_B = 1.38064900e-16
PI = 3.14159265358979324

# IonizationModel gas-phase defaults (ionization.hpp).
ZETA = 1.0e-16       # CR rate (production value)
MU_N = 2.33
M_ION = 24.3 * m_H
ALPHA0 = 2.4e-7
ALPHA_TEXP = -0.69
SIGMA_EN = 1.0e-15
SIGV_IN = 1.9e-9


def vbar(T, mass):
    return np.sqrt(8.0 * k_B * T / (PI * mass))


def gas_charges(rho, T):
    """Gas-phase CR<->recombination balance: n_e = n_i = sqrt(zeta n_n / alpha(T))."""
    n_n = rho / (MU_N * m_H)
    alpha = ALPHA0 * (T / 300.0) ** ALPHA_TEXP
    n_i = np.sqrt(ZETA * n_n / alpha)
    return n_n, n_i, n_i  # n_n, n_e, n_i


def tensor(rho, T, B):
    """Independent Pandey&Wardle sigma_{O,H,P} and eta_{O,H,A} (electrons+ions)."""
    n_n, n_e, n_i = gas_charges(rho, T)
    m_n = MU_N * m_H
    nu_e = (m_n / (m_e + m_n)) * n_n * (SIGMA_EN * vbar(T, m_e))
    nu_i = (m_n / (M_ION + m_n)) * n_n * SIGV_IN
    ecB = e_chg * c_light / B

    def beta(Zsigned, mass, nu):
        return (Zsigned * e_chg * B) / (mass * c_light) / nu

    sO = sH = sP = 0.0
    for (nj, Zj, mj, nuj) in [(n_e, -1.0, m_e, nu_e), (n_i, 1.0, M_ION, nu_i)]:
        b = beta(Zj, mj, nuj)
        sO += nj * Zj * b
        sH += nj * Zj / (1.0 + b * b)
        sP += nj * Zj * b / (1.0 + b * b)
    sO *= ecB; sH *= ecB; sP *= ecB
    pref = c_light * c_light / (4.0 * PI)
    sperp2 = sH * sH + sP * sP
    eta_O = pref / sO
    eta_H = pref * sH / sperp2
    eta_A = pref * sP / sperp2 - eta_O
    return dict(sO=sO, sH=sH, sP=sP, sperp2=sperp2, eta_O=eta_O, eta_H=eta_H, eta_A=eta_A,
                n_e=n_e, nu_e=nu_e, nu_i=nu_i, n_i=n_i, M_ION=M_ION)


# ---------------------------------------------------------------------------------------
# GRAIN-INCLUSIVE MODEL (gates 6-8). MRN dust + grain charging, independent formulation.
# ---------------------------------------------------------------------------------------
RHO_GRAIN = 3.0      # grain material density [g/cm^3]
T_SUBL = 1.5e3       # grain sublimation temperature [K]
STICK_E = 1.0
STICK_I = 1.0
X_K = 1.0e-7         # potassium abundance rel. H nuclei
CHI_K = 4.34         # K first ionization potential [eV]
X_H = 0.716          # hydrogen MASS fraction
EV = 1.602176634e-12
H_PL = 6.62607015e-27
XE_FLOOR = 1.0e-20


def mrn_bins(a_min=5.0e-7, a_max=2.5e-5, p=3.5, n_bin=5):
    """MRN bin geometric-mean radii a_k and dust-MASS fractions mw_k (dn/da ~ a^-p).

    The bin discretization is a modelling choice shared with the C++ (not the physics under
    test), so gate 8 uses the a_k/mw_k the C++ dump reports; this reproduces them.
    """
    edges = np.logspace(np.log10(a_min), np.log10(a_max), n_bin + 1)

    def integ(lo, hi, q):
        return np.log(hi / lo) if abs(q + 1.0) < 1e-9 else \
            (hi ** (q + 1.0) - lo ** (q + 1.0)) / (q + 1.0)

    a_k = np.sqrt(edges[:-1] * edges[1:])
    mw = np.array([integ(edges[k], edges[k + 1], 3.0 - p) for k in range(n_bin)])
    return a_k, mw / mw.sum()


def saha_thermal(n_n, T):
    """Thermal free electrons from K (4.34 eV) + H (13.598 eV), sharing one n_e.

    Independent solve: brentq on the neutrality residual (C++ uses fixed 64-step bisection).
    """
    le3 = (2.0 * PI * m_e * k_B * T / (H_PL * H_PL)) ** 1.5
    fK = le3 * np.exp(-CHI_K * EV / (k_B * T))
    fH = le3 * np.exp(-13.598 * EV / (k_B * T))
    n_K = X_K * n_n
    n_H = n_n * MU_N * X_H
    hi = n_K + n_H
    if hi <= 0.0:
        return 0.0

    def res(ne):
        return n_K * fK / (fK + ne) + n_H * fH / (fH + ne) - ne

    if res(hi) > 0.0:
        return hi
    lo = hi * 1e-300
    # Cold gas: both Boltzmann factors underflow to exactly 0 -> no thermal electrons at all
    # (residual is -ne everywhere, so there is no sign change to bracket).
    if res(lo) <= 0.0:
        return 0.0
    return brentq(res, lo, hi, xtol=1e-300, rtol=1e-15, maxiter=500)


def grain_charges(rho, T, f_dg=0.01, a_k=None, mw_k=None):
    """Self-consistent (n_e, n_i, Z_k, ng_k) with MRN grains -- INDEPENDENT formulation.

    Physics (identical to ionization.hpp::SolveCharges):
      per-bin charge balance : n_e v_e exp(psi) = n_i v_i (1 - psi),  Z_k = psi tau_k
      ionization balance     : P = alpha n_e n_i + Ki n_i,  Ki = (1-psi) sum_k ng_k pi a_k^2 v_i
      neutrality             : n_e = n_i - sum_k (-Z_k) ng_k
      thermal Saha added on top; grains sublimate above T_subl.

    Algorithm (DIFFERENT from the C++): psi is bin-independent, so with r = n_e/n_i the system
    reduces to one scalar equation solved by brentq -- no relaxation, no outer fixed point.
    """
    if a_k is None:
        a_k, mw_k = mrn_bins()
    a_k = np.asarray(a_k, float)
    mw_k = np.asarray(mw_k, float)

    n_n = rho / (MU_N * m_H)
    alpha = ALPHA0 * (T / 300.0) ** ALPHA_TEXP
    P = ZETA * n_n
    ve = vbar(T, m_e)
    vi = vbar(T, M_ION)

    rho_d = f_dg * (n_n * MU_N * m_H)
    m_g = (4.0 / 3.0) * PI * a_k ** 3 * RHO_GRAIN
    ng = np.zeros_like(a_k) if T >= T_SUBL else mw_k * rho_d / m_g
    tau = a_k * k_B * T / (e_chg * e_chg)

    if ng.sum() <= 0.0:                      # no grains: pure gas-phase balance
        n_i = np.sqrt(P / alpha)
        n_e = n_i
        Zk = np.zeros_like(a_k)
    else:
        S = float((tau * ng).sum())                       # sum tau_k ng_k
        K0 = float((ng * PI * a_k ** 2 * vi * STICK_I).sum())

        def psi_of_r(r):
            """Solve r (ve/vi) e^psi = 1 - psi for psi (<0); monotonic -> bracketed."""
            c = r * ve / vi

            def f(psi):
                return c * np.exp(np.clip(psi, -200.0, 50.0)) - (1.0 - psi)

            lo, hi = -200.0, 0.0
            if f(hi) <= 0.0:                 # psi >= 0 branch (very electron-poor)
                return brentq(f, 0.0, 50.0, rtol=1e-15, maxiter=500)
            return brentq(f, lo, hi, rtol=1e-15, maxiter=500)

        def residual(r):
            psi = psi_of_r(r)
            n_i = -psi * S / (1.0 - r)       # from neutrality  n_i(1-r) = Q = -psi S
            return alpha * r * n_i * n_i + (1.0 - psi) * K0 * n_i - P

        # physical bracket: psi < 0 (grains negative) requires r > v_i/v_e
        r_lo = (vi / ve) * (1.0 + 1e-12)
        r = brentq(residual, r_lo, 1.0 - 1e-14, xtol=1e-18, rtol=8.9e-16, maxiter=500)
        psi = psi_of_r(r)
        n_i = -psi * S / (1.0 - r)
        n_e = r * n_i
        Zk = psi * tau

    n_th = saha_thermal(n_n, T)
    n_e += n_th
    n_i += n_th
    n_e = max(n_e, XE_FLOOR * n_n)
    return n_e, n_i, Zk, ng


def tensor_grains(rho, T, B, f_dg=0.01, a_k=None, mw_k=None):
    """Pandey&Wardle sigma/eta including charged MRN grain bins (independent reference)."""
    if a_k is None:
        a_k, mw_k = mrn_bins()
    a_k = np.asarray(a_k, float)
    n_e, n_i, Zk, ng = grain_charges(rho, T, f_dg, a_k, mw_k)
    n_n = rho / (MU_N * m_H)
    m_n = MU_N * m_H
    nu_e = (m_n / (m_e + m_n)) * n_n * (SIGMA_EN * vbar(T, m_e))
    nu_i = (m_n / (M_ION + m_n)) * n_n * SIGV_IN
    vn = vbar(T, m_n)
    ecB = e_chg * c_light / B

    def beta(Zsigned, mass, nu):
        return (Zsigned * e_chg * B) / (mass * c_light) / nu

    sO = sH = sP = 0.0
    for (nj, Zj, mj, nuj) in [(n_e, -1.0, m_e, nu_e), (n_i, 1.0, M_ION, nu_i)]:
        b = beta(Zj, mj, nuj)
        sO += nj * Zj * b
        sH += nj * Zj / (1.0 + b * b)
        sP += nj * Zj * b / (1.0 + b * b)
    for k in range(len(a_k)):
        if ng[k] <= 0.0 or abs(Zk[k]) < 1e-30:
            continue
        m_gk = (4.0 / 3.0) * PI * a_k[k] ** 3 * RHO_GRAIN
        nu_g = (m_n / (m_gk + m_n)) * n_n * (PI * a_k[k] ** 2 * vn)
        bg = beta(Zk[k], m_gk, nu_g)
        sO += ng[k] * Zk[k] * bg
        sH += ng[k] * Zk[k] / (1.0 + bg * bg)
        sP += ng[k] * Zk[k] * bg / (1.0 + bg * bg)
    sO *= ecB; sH *= ecB; sP *= ecB
    pref = c_light * c_light / (4.0 * PI)
    sperp2 = sH * sH + sP * sP
    eta_O = pref / sO
    eta_H = pref * sH / sperp2
    eta_A = pref * sP / sperp2 - eta_O
    return dict(n_e=n_e, n_i=n_i, Zk=Zk, ng=ng, sO=sO, sH=sH, sP=sP,
                eta_O=eta_O, eta_H=eta_H, eta_A=eta_A)


def scalar_sigma(rho, T):
    """Scalar (unmagnetized) conductivity sum_j n_j Z_j^2 e^2/(m_j nu_j)."""
    n_n, n_e, n_i = gas_charges(rho, T)
    m_n = MU_N * m_H
    nu_e = (m_n / (m_e + m_n)) * n_n * (SIGMA_EN * vbar(T, m_e))
    nu_i = (m_n / (M_ION + m_n)) * n_n * SIGV_IN
    return (n_e * e_chg**2 / (m_e * nu_e)) + (n_i * e_chg**2 / (M_ION * nu_i))


def main():
    fails = 0
    rhos = np.logspace(-18, -6, 25)
    Ts = np.array([15.0, 30.0, 100.0])
    Bs = np.array([1e-6, 1e-4, 1e-2, 1.0])  # G

    # Gate 1 & 2: sigma_O == scalar conductivity (B cancels); eta_O = c^2/(4pi sigma_O).
    print("--- GATE 1/2: sigma_O == scalar conductivity (B-independent); eta_O identity ---")
    w_sO = w_etaO = 0.0
    for rho in rhos:
        for T in Ts:
            ss = scalar_sigma(rho, T)
            for B in Bs:
                d = tensor(rho, T, B)
                w_sO = max(w_sO, abs(d['sO'] - ss) / ss)
                w_etaO = max(w_etaO, abs(d['eta_O'] - c_light**2 / (4 * PI * ss)) /
                             (c_light**2 / (4 * PI * ss)))
    g12 = w_sO < 1e-12 and w_etaO < 1e-12
    print(f"  worst |sigma_O - sigma_scalar|/sigma = {w_sO:.2e} ; eta_O identity = {w_etaO:.2e}"
          f"  {'PASS' if g12 else 'FAIL'}")
    fails += not g12

    # Gate 3: parity in B (eta_O,eta_A even; eta_H odd).
    print("--- GATE 3: parity in B (eta_O,eta_A even; eta_H odd -> Hall sign) ---")
    w_par = 0.0
    for rho in [1e-16, 1e-13, 1e-10]:
        for T in Ts:
            for B in Bs:
                p = tensor(rho, T, +B); m = tensor(rho, T, -B)
                w_par = max(w_par, abs(p['eta_O'] - m['eta_O']) / abs(p['eta_O']))
                w_par = max(w_par, abs(p['eta_A'] - m['eta_A']) / (abs(p['eta_A']) + 1e-300))
                w_par = max(w_par, abs(p['eta_H'] + m['eta_H']) / (abs(p['eta_H']) + 1e-300))
    g3 = w_par < 1e-12
    print(f"  worst parity residual = {w_par:.2e}  {'PASS' if g3 else 'FAIL'}")
    fails += not g3

    # Gate 4: tensor algebra sigma_perp^2 = sigma_H^2 + sigma_P^2 (guards assembly).
    print("--- GATE 4: sigma_perp^2 == sigma_H^2 + sigma_P^2 ---")
    w_alg = 0.0
    for rho in rhos:
        d = tensor(rho, 15.0, 1e-4)
        w_alg = max(w_alg, abs(d['sperp2'] - (d['sH']**2 + d['sP']**2)) / d['sperp2'])
    g4 = w_alg < 1e-12
    print(f"  worst residual = {w_alg:.2e}  {'PASS' if g4 else 'FAIL'}")
    fails += not g4

    # Gate 5: canonical crossover -- eta_A dominates at low rho, Ohmic rises at high rho.
    print("--- GATE 5: canonical eta(rho) crossover (ambipolar -> Ohmic) ---")
    # flux-frozen field B = B0 (rho/rho0)^(1/2), collapse track, T=15 K.
    rho0 = 5.467e-19; B0 = 4.98e-5 * 0.15  # ~ B0z=0.15 code -> cgs
    lo = tensor(1e-17, 15.0, B0 * (1e-17 / rho0) ** 0.5)
    hi = tensor(1e-8, 15.0, B0 * (1e-8 / rho0) ** 0.5)
    aOlo = abs(lo['eta_A']) / abs(lo['eta_O'])   # ambipolar/Ohmic at low rho (expect >> 1)
    aOhi = abs(hi['eta_A']) / abs(hi['eta_O'])   # at high rho (expect << low)
    g5 = aOlo > 1.0 and aOhi < aOlo
    print(f"  eta_A/eta_O: low-rho={aOlo:.2e}  high-rho={aOhi:.2e}  "
          f"(ambipolar dominant low, drops toward Ohmic)  {'PASS' if g5 else 'FAIL'}")
    fails += not g5

    # ---------------- GRAIN-INCLUSIVE GATES (Phase-3 remaining item) ----------------
    # Gate 6: self-consistency of the independent grain solve (neutrality + ionization
    # balance must hold on the returned state, to roundoff).
    print("--- GATE 6: grain solve self-consistency (neutrality + ionization balance) ---")
    a_k, mw_k = mrn_bins()
    w_neu = w_ion = 0.0
    for rho in [1e-18, 1e-16, 1e-14, 1e-12, 1e-10, 1e-8]:
        for T in [10.0, 15.0, 50.0, 300.0, 1000.0]:
            n_e, n_i, Zk, ng = grain_charges(rho, T, 0.01, a_k, mw_k)
            n_n = rho / (MU_N * m_H)
            n_th = saha_thermal(n_n, T)
            # subtract the thermal channel (added on top of the CR/grain balance)
            ne_c, ni_c = n_e - n_th, n_i - n_th
            Q = float((-Zk * ng).sum())
            w_neu = max(w_neu, abs(ni_c - ne_c - Q) / max(abs(ni_c), 1e-300))
            alpha = ALPHA0 * (T / 300.0) ** ALPHA_TEXP
            P = ZETA * n_n
            vi = vbar(T, M_ION)
            psi = Zk[0] / (a_k[0] * k_B * T / (e_chg * e_chg)) if ng.sum() > 0 else 0.0
            Ki = float((ng * PI * a_k ** 2 * vi * STICK_I).sum()) * (1.0 - psi)
            w_ion = max(w_ion, abs(alpha * ne_c * ni_c + Ki * ni_c - P) / P)
    g6 = w_neu < 1e-10 and w_ion < 1e-8
    print(f"  worst neutrality residual = {w_neu:.2e} ; ionization-balance residual = "
          f"{w_ion:.2e}  {'PASS' if g6 else 'FAIL'}")
    fails += not g6

    # Gate 7: grain physics structure -- negative charging, Z ~ a, electron depletion,
    # sublimation above T_subl.
    print("--- GATE 7: grain physics structure (sign, Z~a, e- depletion, sublimation) ---")
    n_e, n_i, Zk, ng = grain_charges(1e-14, 15.0, 0.01, a_k, mw_k)
    neg = bool(np.all(Zk < 0.0))
    # Z_k = psi tau_k with tau ~ a  =>  Z_k/a_k constant across bins
    ratio = Zk / a_k
    lin = float(np.max(np.abs(ratio - ratio[0])) / abs(ratio[0]))
    # grains suppress free electrons vs the gas-phase-only balance
    _, ne_gas, _ = gas_charges(1e-14, 15.0)
    depl = n_e < ne_gas
    # above T_subl grains are gone
    _, _, Zk_hot, ng_hot = grain_charges(1e-14, 2000.0, 0.01, a_k, mw_k)
    subl = bool(np.all(ng_hot == 0.0))
    g7 = neg and lin < 1e-10 and depl and subl
    print(f"  Z_k<0 all bins={neg} ; |Z_k/a_k| spread={lin:.2e} (expect ~0, Z=psi*tau~a) ; "
          f"n_e {n_e:.3e} < gas-only {ne_gas:.3e} = {depl} ; sublimated>T_subl={subl}  "
          f"{'PASS' if g7 else 'FAIL'}")
    fails += not g7

    # Gate 8: direct C++-vs-Python grain-inclusive numerical cross-check.
    csv = None
    if "--xcheck-csv" in sys.argv:
        csv = sys.argv[sys.argv.index("--xcheck-csv") + 1]
    print("--- GATE 8: C++-vs-Python GRAIN-INCLUSIVE cross-check (charge state + eta) ---")
    if csv is None:
        print("  SKIPPED (pass --xcheck-csv <dump>; build xcheck_conductivity.cpp first)")
    else:
        import csv as _csv
        a_dump, mw_dump, rows = [], [], []
        with open(csv) as fh:
            for line in fh:
                if line.startswith("# bin"):
                    p = line.split()
                    a_dump.append(float(p[3].split("=")[1]))
                    mw_dump.append(float(p[4].split("=")[1]))
                elif not line.startswith("#"):
                    rows.append(line)
        rd = list(_csv.DictReader(rows))
        a_dump = np.array(a_dump); mw_dump = np.array(mw_dump)
        nb = len(a_dump)
        w = {"n_e": 0.0, "n_i": 0.0, "Z": 0.0, "eta_O": 0.0, "eta_H": 0.0, "eta_A": 0.0}
        worst_at = {}
        for r in rd:
            rho, T, B = float(r["rho"]), float(r["T"]), float(r["B"])
            d = tensor_grains(rho, T, B, 0.01, a_dump, mw_dump)
            for key, val in [("n_e", d["n_e"]), ("n_i", d["n_i"]),
                             ("eta_O", d["eta_O"]), ("eta_H", d["eta_H"]),
                             ("eta_A", d["eta_A"])]:
                ref = float(r[key])
                rel = abs(val - ref) / max(abs(ref), 1e-300)
                if rel > w[key]:
                    w[key] = rel; worst_at[key] = (rho, T, B)
            for k in range(nb):
                ref = float(r[f"Z{k}"])
                if abs(ref) > 1e-30:
                    rel = abs(d["Zk"][k] - ref) / abs(ref)
                    if rel > w["Z"]:
                        w["Z"] = rel; worst_at["Z"] = (rho, T, B, k)
        for key in ("n_e", "n_i", "Z", "eta_O", "eta_H", "eta_A"):
            print(f"    worst rel diff {key:6s} = {w[key]:.3e}   at {worst_at.get(key)}")
        # The charge state and eta_O are exact to solver precision. eta_H and eta_A are
        # cancellation-limited: sigma_H is a signed sum over e/i/grains whose net is ~1e-13
        # of the summed |term| magnitude near the grain-induced Hall sign reversal (a ~1e13x
        # relative-error amplifier), and eta_A = c^2/4pi sigma_P/sigma_perp^2 - eta_O
        # subtracts two nearly equal numbers. Both are held to the Phase-3 <1% criterion;
        # where they are least accurate they are also physically subdominant (at the worst
        # point eta_H is 5 decades below eta_O).
        g8 = (w["n_e"] < 1e-8 and w["n_i"] < 1e-8 and w["Z"] < 1e-8 and
              w["eta_O"] < 1e-8 and w["eta_H"] < 1e-2 and w["eta_A"] < 1e-2)
        print(f"  {len(rd)} sample points, {nb} grain bins  {'PASS' if g8 else 'FAIL'}")
        fails += not g8

    print(f"\n{'ALL GATES PASS' if fails == 0 else 'GATE FAILURES'} ({fails} failures)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
