#!/usr/bin/env python3
"""
Jeans analysis for AMR runs, both stable (oscillation) and unstable (growth).
Projects density onto the seeded sin(kx) mode via volume-weighted integration
over all meshblocks at any refinement level.

Usage: python3 scripts/jeans_analysis_amr_v2.py <stable|unstable> <prefix>
"""

import sys, os, glob
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..",
    "external/parthenon/scripts/python/packages/parthenon_tools/parthenon_tools"))
import phdf


def project_sin_mode(fn, k_wave, Lx, Ly, Lz):
    """
    Volume-weighted projection of rho onto sin(k*x) and cos(k*x).
    For rho = rho_bar + A sin(kx + phi), returns (A*cos(phi), A*sin(phi)) in
    the (a_sin, a_cos) components, and amp = hypot(a_sin, a_cos).
    Also returns the volume-averaged rho_mean for conservation diagnostics.
    """
    f = phdf.phdf(fn)
    prim = f.Get("prim")
    locs = np.array(f.BlockBounds)
    nb = f.NumBlocks

    try:
        sizes = np.array(f.MeshBlockSize)
    except Exception:
        sizes = np.array([f.MeshBlockSize] * 3)
    nx_cell = int(sizes[0])
    ny_cell = int(sizes[1])
    nz_cell = int(sizes[2])

    V_total = Lx * Ly * Lz

    numer_sin = 0.0
    numer_cos = 0.0
    mean_rho_num = 0.0
    vol_total = 0.0

    for b in range(nb):
        x1lo, x1hi = locs[b, 0], locs[b, 1]
        x2lo, x2hi = locs[b, 2], locs[b, 3]
        x3lo, x3hi = locs[b, 4], locs[b, 5]
        dx = (x1hi - x1lo) / nx_cell
        rho_flat = prim[b, 0, :]
        rho_3d = rho_flat.reshape((nz_cell, ny_cell, nx_cell))
        rho_1d = rho_3d.mean(axis=(0, 1))
        xc = x1lo + (np.arange(nx_cell) + 0.5) * dx
        A_transverse = (x2hi - x2lo) * (x3hi - x3lo)
        numer_sin += np.sum(rho_1d * np.sin(k_wave * xc)) * dx * A_transverse
        numer_cos += np.sum(rho_1d * np.cos(k_wave * xc)) * dx * A_transverse
        mean_rho_num += np.sum(rho_1d) * dx * A_transverse
        vol_total += (x1hi - x1lo) * A_transverse

    a_sin = 2.0 / V_total * numer_sin
    a_cos = 2.0 / V_total * numer_cos
    rho_mean = mean_rho_num / vol_total
    return f.Time, a_sin, a_cos, rho_mean


def count_blocks_by_level(fn):
    """Return (nblocks_total, nblocks_at_each_level) for diagnostic."""
    f = phdf.phdf(fn)
    try:
        levels = np.array(f.Levels)
        unique, counts = np.unique(levels, return_counts=True)
        return f.NumBlocks, dict(zip(unique.tolist(), counts.tolist()))
    except Exception:
        return f.NumBlocks, {}


def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    regime = sys.argv[1]
    prefix = sys.argv[2]

    files = sorted(glob.glob(f"{prefix}.*.phdf"))
    files = [f for f in files if ".final." not in f]
    if not files:
        print(f"No files matching {prefix}.*.phdf"); sys.exit(1)
    print(f"Found {len(files)} snapshots")

    rho0 = 1.0
    four_pi_G = 1.0
    cs = 0.1 if regime == "unstable" else 1.0
    Lx = 1.0
    Ly = 0.03125
    Lz = 0.03125
    nwave = 1
    k_wave = 2.0 * np.pi * nwave / Lx

    omega2 = k_wave**2 * cs**2 - four_pi_G * rho0
    if omega2 > 0:
        omega_theory = np.sqrt(omega2)
        print(f"Analytic: k={k_wave:.4f}, omega^2={omega2:.4f}, mode=oscillating")
        print(f"  omega_theory = {omega_theory:.6f}  (period {2*np.pi/omega_theory:.6f})")
    else:
        gamma_theory = np.sqrt(-omega2)
        print(f"Analytic: k={k_wave:.4f}, omega^2={omega2:.4f}, mode=growing")
        print(f"  gamma_theory = {gamma_theory:.6f}  (e-fold {1/gamma_theory:.6f})")

    # AMR diagnostic: block count over time
    nb_first, lvls_first = count_blocks_by_level(files[0])
    nb_mid, lvls_mid = count_blocks_by_level(files[len(files)//2])
    nb_last, lvls_last = count_blocks_by_level(files[-1])
    print(f"\nBlock counts: t=0 -> {nb_first} blocks {lvls_first}, "
          f"t=mid -> {nb_mid} blocks {lvls_mid}, "
          f"t=end -> {nb_last} blocks {lvls_last}")
    # If block counts change, adaptive refinement is actually working
    if nb_first != nb_mid or nb_mid != nb_last:
        print("  (Block count varies — adaptive refinement is active.)")
    else:
        print("  (Block count constant — adaptive refinement did not trigger.)")

    times, a_sins, a_coss, rho_means = [], [], [], []
    for fn in files:
        t, a_sin, a_cos, rho_mean = project_sin_mode(fn, k_wave, Lx, Ly, Lz)
        times.append(t); a_sins.append(a_sin); a_coss.append(a_cos)
        rho_means.append(rho_mean)

    times = np.array(times)
    A_sin = np.array(a_sins)
    A_cos = np.array(a_coss)
    amps = np.hypot(A_sin, A_cos)
    rho_means = np.array(rho_means)

    print(f"\nTime range: [{times[0]:.3f}, {times[-1]:.3f}]")
    print(f"Amp(rho mode) range: [{amps.min():.3e}, {amps.max():.3e}]")
    print(f"rho mean over time: min={rho_means.min():.6e}, max={rho_means.max():.6e}")
    print(f"  (should be very close to {rho0:.4f})")

    for idx in [0, len(times)//4, len(times)//2, 3*len(times)//4, len(times)-1]:
        print(f"  t = {times[idx]:6.3f},  amp = {amps[idx]:.4e}")

    if omega2 < 0:
        # Growing mode: fit log-linear in the late window
        i0 = int(0.70 * len(times))
        i1 = int(0.95 * len(times))
        t_fit = times[i0:i1]
        logA = np.log(amps[i0:i1])
        p = np.polyfit(t_fit, logA, 1)
        gamma_measured = p[0]
        err = 100.0 * abs(gamma_measured - gamma_theory) / gamma_theory
        print(f"\nFit over t in [{t_fit[0]:.3f}, {t_fit[-1]:.3f}]:")
        print(f"  Measured growth rate: {gamma_measured:.6f}")
        print(f"  Theory:               {gamma_theory:.6f}")
        print(f"  Error:                {err:.3f}%")
        print(">>> " + ("PASS (error < 2%)" if err < 2.0 else "FAIL"))
    else:
        # Oscillating mode: FFT of the sin-projection signal.
        # Use A_sin (signed) as the time signal; subtract its mean.
        S = A_sin - np.mean(A_sin)
        dt = np.mean(np.diff(times))
        freqs = np.fft.rfftfreq(len(S), d=dt)
        spectrum = np.abs(np.fft.rfft(S))
        # Skip DC bin
        if len(spectrum) < 2:
            print("Too few snapshots for FFT"); sys.exit(1)
        peak_idx = np.argmax(spectrum[1:]) + 1
        f_peak = freqs[peak_idx]
        omega_measured = 2.0 * np.pi * f_peak

        err = 100.0 * abs(omega_measured - omega_theory) / omega_theory
        print(f"\nFFT of sin-mode projection:")
        print(f"  Peak frequency bin:   {peak_idx} / {len(spectrum)}")
        print(f"  Measured omega:       {omega_measured:.6f}")
        print(f"  Theory:               {omega_theory:.6f}")
        print(f"  Error:                {err:.3f}%")
        print(">>> " + ("PASS (error < 5%)" if err < 5.0 else "FAIL"))
        # Note: stable AMR pass threshold is 5%, not 2%, because
        # refinement events add noise beyond pure truncation error.


if __name__ == "__main__":
    main()
