#!/usr/bin/env python3
"""
Overlay rho_max(t) and B_max(t) between AthenaPK and Athena++ BE collapse runs.

Usage: python3 compare_codes.py <athenapk_dir> <athena++_dir>
  AthenaPK  dir: contains parthenon.out0.*.phdf
  Athena++  dir: contains *.athdf (typically Collapse.out2.*.athdf)
"""
import sys
import glob
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_athenapk(run_dir):
    """AthenaPK: prim[nblocks, 9, nz, ny, nx] with [rho, v1, v2, v3, P, B1, B2, B3, psi]"""
    files = sorted(glob.glob(f"{run_dir}/parthenon.out0.?????.phdf"))
    times, rho_max, B_max = [], [], []
    for fn in files:
        with h5py.File(fn, "r") as f:
            t = float(f["Info"].attrs["Time"])
            prim = f["prim"][:]          # (nblocks, nvars, nz, ny, nx)
            rho = prim[:, 0, ...]
            B1 = prim[:, 5, ...]
            B2 = prim[:, 6, ...]
            B3 = prim[:, 7, ...]
            Bm = np.sqrt(B1**2 + B2**2 + B3**2)
            times.append(t)
            rho_max.append(float(rho.max()))
            B_max.append(float(Bm.max()))
    return np.array(times), np.array(rho_max), np.array(B_max), files


def read_athenapp(run_dir):
    """Athena++: prim[7, nblocks, nz, ny, nx] with [rho, press, v1, v2, v3, phi, defect]
                 B[3, nblocks, nz, ny, nx] with [Bcc1, Bcc2, Bcc3]"""
    files = sorted(glob.glob(f"{run_dir}/*.athdf"))
    times, rho_max, B_max = [], [], []
    for fn in files:
        with h5py.File(fn, "r") as f:
            t = float(f.attrs["Time"])
            prim = f["prim"][:]          # (nvars=7, nblocks, nz, ny, nx)
            B = f["B"][:]                # (3, nblocks, nz, ny, nx)
            rho = prim[0, ...]
            Bm = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
            times.append(t)
            rho_max.append(float(rho.max()))
            B_max.append(float(Bm.max()))
    return np.array(times), np.array(rho_max), np.array(B_max), files


def main(apk_dir, ath_dir):
    print("Reading AthenaPK...")
    t1, rhomax1, Bmax1, files1 = read_athenapk(apk_dir)
    print(f"  {len(t1)} snapshots, t=[{t1[0]:.3f}, {t1[-1]:.3f}]")

    print("Reading Athena++...")
    t2, rhomax2, Bmax2, files2 = read_athenapp(ath_dir)
    print(f"  {len(t2)} snapshots, t=[{t2[0]:.3f}, {t2[-1]:.3f}]")

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    ax.semilogy(t1, rhomax1, 'b-',  label='AthenaPK', linewidth=2)
    ax.semilogy(t2, rhomax2, 'r--', label='Athena++', linewidth=2)
    ax.set_xlabel('time (code units)')
    ax.set_ylabel(r'$\rho_{\max}$')
    ax.legend()
    ax.set_title('Maximum density evolution')
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(t1, Bmax1, 'b-',  label='AthenaPK', linewidth=2)
    ax.plot(t2, Bmax2, 'r--', label='Athena++', linewidth=2)
    ax.set_xlabel('time (code units)')
    ax.set_ylabel(r'$|B|_{\max}$')
    ax.legend()
    ax.set_title('Maximum field strength evolution')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('comparison.png', dpi=120, bbox_inches='tight')
    print("\nSaved comparison.png")

    # --- Table: side-by-side rho_max ---
    print("\n--- rho_max(t) comparison ---")
    print(f"{'t':>8} {'AthenaPK':>12} {'Athena++':>12} {'ratio':>8}")
    rhomax2_i = np.interp(t1, t2, rhomax2)
    for t_, r1, r2 in zip(t1[::4], rhomax1[::4], rhomax2_i[::4]):
        ratio = r1 / r2 if r2 > 0 else np.nan
        print(f"{t_:8.3f} {r1:12.4e} {r2:12.4e} {ratio:8.3f}")

    # --- Table: side-by-side B_max ---
    print("\n--- B_max(t) comparison ---")
    print(f"{'t':>8} {'AthenaPK':>12} {'Athena++':>12} {'ratio':>8}")
    Bmax2_i = np.interp(t1, t2, Bmax2)
    for t_, b1, b2 in zip(t1[::4], Bmax1[::4], Bmax2_i[::4]):
        ratio = b1 / b2 if b2 > 0 else np.nan
        print(f"{t_:8.3f} {b1:12.4e} {b2:12.4e} {ratio:8.3f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 compare_codes.py <athenapk_dir> <athena++_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
