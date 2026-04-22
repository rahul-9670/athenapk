#!/usr/bin/env python3
"""
Jeans dispersion analysis, v2: handles flattened (nblocks, nvar, nflat) phdf shape.

Usage: python3 scripts/jeans_analysis_v2.py <stable|unstable> <prefix>
"""

import sys, os, glob
import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), "..", "/Users/user/athenapk/external/parthenon/scripts/python/packages/parthenon_tools/parthenon_tools/"))
import phdf


def load_rho_1d(fn, nx_block, ny_block, nz_block):
    """Load rho and return a 1D array along x after averaging over y,z.

    Assumes single meshblock spanning full domain in x.
    phdf returns prim with shape (nblocks, nvar, ncells_flat) where
    ncells_flat = nz * ny * nx with x being the fastest-moving index.
    """
    f = phdf.phdf(fn)
    prim = f.Get("prim")
    # single block
    rho_flat = prim[0, 0, :]   # (nvar index 0 = density)
    assert rho_flat.size == nx_block * ny_block * nz_block, \
        f"expected {nx_block*ny_block*nz_block}, got {rho_flat.size}"
    # Reshape: Parthenon stores cells in (k, j, i) order, x fastest
    rho_3d = rho_flat.reshape((nz_block, ny_block, nx_block))
    # Average over y (axis=1) and z (axis=0)
    rho_1d = rho_3d.mean(axis=(0, 1))
    return f.Time, rho_1d


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    regime = sys.argv[1]
    prefix = sys.argv[2]

    files = sorted(glob.glob(f"{prefix}.*.phdf"))
    # Filter out .final. which has no clean timestamp ordering with numeric ones
    files = [f for f in files if ".final." not in f]
    if not files:
        print(f"No files matching {prefix}.*.phdf")
        sys.exit(1)
    print(f"Found {len(files)} snapshots")

    # Problem parameters — hardcoded to match jeans_*.in
    rho0 = 1.0
    four_pi_G = 1.0
    if regime == "stable":
        cs = 1.0
    elif regime == "unstable":
        cs = 0.1
    else:
        print("regime must be 'stable' or 'unstable'"); sys.exit(1)
    Lx = 1.0
    nwave = 1
    k = 2.0 * np.pi * nwave / Lx

    # Mesh from input file
    nx_block = 128
    ny_block = 4
    nz_block = 4

    # Analytic dispersion
    omega2 = k*k*cs*cs - four_pi_G * rho0
    if omega2 > 0:
        omega_theory = np.sqrt(omega2)
        mode = "oscillating"
        print(f"Analytic: k={k:.4f}, omega^2 = {omega2:.4f}, mode = {mode}")
        print(f"  Analytic omega = {omega_theory:.6f} (period {2*np.pi/omega_theory:.6f})")
    else:
        gamma_theory = np.sqrt(-omega2)
        mode = "growing"
        print(f"Analytic: k={k:.4f}, omega^2 = {omega2:.4f}, mode = {mode}")
        print(f"  Analytic growth rate = {gamma_theory:.6f} (e-fold {1/gamma_theory:.6f})")

    # Extract
    times, amps_sine = [], []
    x = (np.arange(nx_block) + 0.5) * (Lx / nx_block)
    sin_kx = np.sin(k * x)
    cos_kx = np.cos(k * x)

    for fn in files:
        t, rho_1d = load_rho_1d(fn, nx_block, ny_block, nz_block)
        # Project onto sin and cos
        drho = rho_1d - rho_1d.mean()
        a_sin = (2.0 / nx_block) * np.sum(drho * sin_kx)
        a_cos = (2.0 / nx_block) * np.sum(drho * cos_kx)
        amp = np.hypot(a_sin, a_cos)
        times.append(t)
        amps_sine.append(amp)

    t = np.array(times)
    A = np.array(amps_sine)

    print(f"Time range: [{t.min():.3f}, {t.max():.3f}]")
    print(f"Amp(rho mode) range: [{A.min():.3e}, {A.max():.3e}]")

    # Print the growth trajectory at a few points
    for i in [0, len(t)//4, len(t)//2, 3*len(t)//4, len(t)-1]:
        print(f"  t = {t[i]:6.3f},  amp = {A[i]:.4e}")

    if mode == "growing":
        # Skip initial transient (first 20%) and possible saturation (last 10%)
        i0, i1 = int(0.70 * len(t)), int(0.95 * len(t))
        if i1 - i0 < 5:
            print("Not enough points for fit")
            sys.exit(1)
        logA = np.log(A[i0:i1])
        p = np.polyfit(t[i0:i1], logA, 1)
        gamma_measured = p[0]
        err = 100.0 * abs(gamma_measured - gamma_theory) / gamma_theory
        print(f"\nFit over t in [{t[i0]:.3f}, {t[i1-1]:.3f}]:")
        print(f"  Measured growth rate: {gamma_measured:.6f}")
        print(f"  Theory:               {gamma_theory:.6f}")
        print(f"  Error:                {err:.3f}%")
        print(">>> " + ("PASS (error < 2%)" if err < 2.0 else "FAIL"))
    else:
        if len(t) < 8:
            print("Too few snapshots for FFT"); sys.exit(1)
        dt = np.mean(np.diff(t))
        freqs = np.fft.rfftfreq(len(A), d=dt)
        # Use signed amplitude (sin projection) for cleaner oscillation signal
        times2, signed = [], []
        for fn in files:
            t_, rho_1d = load_rho_1d(fn, nx_block, ny_block, nz_block)
            drho = rho_1d - rho_1d.mean()
            a_sin = (2.0 / nx_block) * np.sum(drho * sin_kx)
            times2.append(t_); signed.append(a_sin)
        S = np.array(signed) - np.mean(signed)
        spectrum = np.abs(np.fft.rfft(S))
        peak_idx = np.argmax(spectrum[1:]) + 1
        f_peak = freqs[peak_idx]
        omega_measured = 2.0 * np.pi * f_peak
        err = 100.0 * abs(omega_measured - omega_theory) / omega_theory
        print(f"\n  Measured omega:  {omega_measured:.6f}")
        print(f"  Theory:          {omega_theory:.6f}")
        print(f"  Error:           {err:.3f}%")
        print(">>> " + ("PASS (error < 2%)" if err < 2.0 else "FAIL"))


if __name__ == "__main__":
    main()
