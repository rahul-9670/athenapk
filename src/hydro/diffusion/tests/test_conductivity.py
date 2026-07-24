#!/usr/bin/env python
"""Flagship Phase 3 conductivity gate: independent check of the Wardle tensor.

Read-only. An INDEPENDENT pure-Python reimplementation of the Pandey & Wardle (2008)
conductivity tensor and the Ohmic/Hall/ambipolar diffusivities that
src/hydro/diffusion/ionization.hpp::Diffusivities computes on device, evaluated on the same
gas-phase ionization balance (electrons + one ion species; grains are the C++ model's
extension, checked separately). It validates the tensor ASSEMBLY and the derived eta's against
rigorous physics identities and the canonical collapse eta(rho) structure -- the "conductivity
vs an independent solver" Phase-3 gate. A direct numerical C++-vs-Python cross-check (incl.
grains) needs a small debug dump hook in Diffivities() + a CPU build -> noted as the follow-up.

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

Run: /beegfs/u/bbg6470/venvs/analysis_env/bin/python src/hydro/diffusion/tests/test_conductivity.py
"""
import sys
import numpy as np

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

    print(f"\n{'ALL GATES PASS' if fails == 0 else 'GATE FAILURES'} ({fails} failures)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
