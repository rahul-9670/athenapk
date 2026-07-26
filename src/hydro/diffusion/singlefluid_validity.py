#!/usr/bin/env python3
"""Phase 6 gate — single-fluid (non-ideal MHD) VALIDITY MAP.

Answers the Phase-6 question directly: is the single-fluid strong-coupling approximation
(the Wardle-tensor non-ideal MHD the code evolves) valid throughout the FHC collapse, or is
there a regime where the ions/electrons/grains decouple from the neutrals and a true multifluid
treatment (separate ion inertia, generalized Ohm law) is required?

Physics is grounded in the code's own ionization model (`ionization.hpp`): same constants
(mu_n=2.33, m_ion=24.3 m_H, sigv_in=1.9e-9 Langevin), same signed Hall parameter
beta_j = Z_j e B / (m_j c) / nu_j, same collision frequency nu_jn = (m_n/(m_j+m_n)) n_n <sigv>.

Key simplification (why this is robust): the ion-neutral collision rate is Langevin (nu_in
propto n_n, T-INDEPENDENT), so the ion Hall parameter beta_i and the neutral-ion coupling
nu_ni/omega_ff depend ONLY on (rho, B, x_e) -- no temperature needed. Electron/grain terms are
reported where T is available but do not gate the single-fluid verdict.

Validity criterion (per cell) -- CORRECT strong-coupling condition:
  * The two-fluid (ion+neutral) system reduces to single-fluid non-ideal MHD (the Wardle-tensor
    AD/Hall/Ohm the code evolves) when ION INERTIA is negligible: the ions reach terminal drift
    within a dynamical time.  chi_ion = nu_in / omega_ff >> 1, with nu_in = (m_n/(m_i+m_n)) n_n
    <sigv> the ION-neutral collision rate (propto n_n, LARGE in dense gas -- THIS is why
    single-fluid AD holds throughout star formation) and omega_ff = sqrt(4 pi G rho).
    chi_ion >> 1 => single-fluid VALID; chi_ion <~ 1 => ion inertia => TRUE MULTIFLUID.
  * PITFALL (guarded against): nu_ni/omega_ff (NEUTRAL-ion, propto rho_i) being < 1 does NOT mean
    multifluid -- it only means AD is ACTIVE (neutrals not perfectly flux-frozen), the intended
    single-fluid non-ideal regime. Using nu_ni as the criterion falsely flags ~all dense gas.
  * Hall parameter |beta_i| labels the diffusion regime (resistive/Hall/ambipolar); all are
    single-fluid where chi_ion >> 1.
"""
import sys, glob, argparse
import numpy as np
import h5py

# --- CGS constants (match ionization.hpp) ---
m_H = 1.67262192e-24; e_chg = 4.80320471e-10; c_light = 2.99792458e10
G = 6.674e-8; k_B = 1.380649e-16; m_e = 9.1093837e-28
# --- ionization-model params (ionization.hpp defaults) ---
MU_N = 2.33; M_ION = 24.3 * m_H; SIGV_IN = 1.9e-9  # cm^3/s Langevin
# --- FHC code units ---
RHO0 = 5.467e-19; V0 = 1.9e4; L0 = 2.81e16
B_UNIT = 4.98e-5  # G per code B (Heaviside-Lorentz)


def load(fn):
    with h5py.File(fn, "r") as h:
        t = float(h["Info"].attrs["Time"])
        p = np.array(h["prim"])  # [blk, comp, k,j,i]; 0=rho 1-3=v 4=P 5-7=B 8=psi 9-13=scalars
    rho = p[:, 0].reshape(-1) * RHO0                         # g/cm^3
    B = np.sqrt(p[:, 5]**2 + p[:, 6]**2 + p[:, 7]**2).reshape(-1) * B_UNIT  # G
    xe = p[:, 13].reshape(-1)                                # scalar_4 = x_e (chemistry)
    return t, rho, B, np.clip(xe, 1e-20, None)


def validity(rho, B, xe):
    n_n = rho / (MU_N * m_H)                 # cm^-3
    m_n = MU_N * m_H
    nu_in = (m_n / (M_ION + m_n)) * n_n * SIGV_IN     # ION-neutral collision freq [s^-1] (propto n_n)
    rho_i = xe * n_n * M_ION                          # ion mass density
    nu_ni = (rho_i / rho) * nu_in                     # neutral-ion collision freq (propto rho_i)
    omega_ff = np.sqrt(4.0 * np.pi * G * rho)         # dynamical rate [s^-1]
    # SINGLE-FLUID (strong-coupling) validity: ion inertia negligible <=> ions equilibrate within
    # a dynamical time <=> nu_in >> omega_ff. THIS is the correct criterion (nu_in propto n_n, so
    # it is easily satisfied in dense gas -- why AD single-fluid MHD holds in star formation).
    chi_ion = nu_in / omega_ff                        # >> 1 => single-fluid VALID
    chi_ni = nu_ni / omega_ff                         # neutral flux-freezing (< 1 => AD active; NOT a failure)
    beta_i = e_chg * np.maximum(B, 1e-20) / (M_ION * c_light) / nu_in  # ion Hall parameter
    return dict(n_n=n_n, nu_in=nu_in, nu_ni=nu_ni, omega_ff=omega_ff,
                chi_ion=chi_ion, chi_ni=chi_ni, beta_i=beta_i)


def report(fn):
    t, rho, B, xe = load(fn)
    d = validity(rho, B, xe)
    ci, cn, bi = d["chi_ion"], d["chi_ni"], d["beta_i"]
    print(f"=== Phase 6 single-fluid validity map ===")
    print(f"file t={t:.4f}  cells={rho.size}  rho[{rho.min():.2e},{rho.max():.2e}] g/cm^3")
    print(f"\nregime by ION strong-coupling chi_ion = nu_in/omega_ff (>>1 => single-fluid valid):")
    frac_ok = np.mean(ci > 100.0); frac_marg = np.mean((ci <= 100.0) & (ci > 10.0))
    frac_bad = np.mean(ci <= 10.0)
    print(f"  chi_ion>100 (single-fluid SOLID)      : {100*frac_ok:6.2f}% of cells")
    print(f"  10<chi<100  (single-fluid ok)         : {100*frac_marg:6.2f}%")
    print(f"  chi_ion<10  (ion inertia -> MULTIFLUID): {100*frac_bad:6.2f}%")
    print(f"\nvalidity vs density (mass-representative bins):")
    print(f"  {'rho[cgs]':>10} {'x_e(med)':>9} {'beta_i(med)':>11} {'chi_ion(med)':>12} {'chi_ni':>9} {'regime':>16}")
    edges = np.logspace(np.log10(rho.min()), np.log10(rho.max()), 9)
    for a, b in zip(edges[:-1], edges[1:]):
        m = (rho >= a) & (rho < b)
        if m.sum() < 5: continue
        bmed = np.median(bi[m]); cimed = np.median(ci[m]); cnmed = np.median(cn[m]); xmed = np.median(xe[m])
        reg = "single-fluid" if cimed > 10 else ("marginal" if cimed > 1 else "MULTIFLUID?")
        hall = "Hall" if 0.3 < abs(bmed) < 3 else ("ambip" if abs(bmed) >= 3 else "resist")
        print(f"  {np.sqrt(a*b):>10.2e} {xmed:>9.2e} {bmed:>11.2e} {cimed:>12.2e} {cnmed:>9.1e} {reg:>10}[{hall}]")
    frac_fail = np.mean(ci <= 1.0)
    verdict = "SINGLE-FLUID NON-IDEAL MHD ADEQUATE (chi_ion>1 everywhere)" if frac_fail < 1e-3 else \
              f"MULTIFLUID needed in {100*frac_fail:.3f}% of cells (chi_ion<1, ion inertia)"
    print(f"\nVERDICT: {verdict}")
    print(f"  min chi_ion over all cells = {ci.min():.2e} (single-fluid needs chi_ion >~ 1)")
    print(f"  (chi_ni<1 at high rho is EXPECTED: AD active, neutrals not flux-frozen -- still single-fluid)")
    return d


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("phdf", help="a science phdf (prim with B + x_e scalar)")
    args = ap.parse_args()
    report(args.phdf)
