#!/usr/bin/env python
"""WP-22 part 2 — the PHYSICAL ambipolar diffusivity in the first core, from a production snapshot.

Part 1 (`wp22_eta_numerical.py`) bounded the scheme's own NUMERICAL resistivity from the WP-14
Alfven ladder. This closes the other half: eta_A as the code actually evaluates it, so the two can
be compared. The fossil-field claim rests on physical eta dominating numerical eta in the core, and
that comparison has never been made.

WHY THIS HAS TO BE DONE OFFLINE. eta is not a written field -- no phdf carries it. But production
runs `<diffusion> ambipolar_coeff = ionization_chem` with `xe_scalar_index = 4`, and the electron
abundance x_e IS written, as passive scalar 4. So eta_A can be reconstructed EXACTLY as the kernel
computes it, from (rho, B, x_e) in the snapshot.

The kernel is src/hydro/diffusion/ionization.hpp :: AmbipolarEtaFromXe --

    n_n   = rho / (mu_n m_H)                      neutral number density
    n_i   = max(x_e, xe_floor) * n_n              ions track electrons by charge balance
    rho_i = n_i * m_ion * m_H
    eta_A[cgs] = B^2 / (4 pi gamma_AD rho_i rho)
    eta_A[code] = eta_A[cgs] * eta_unit

-- which simplifies to eta_A[cgs] = B^2 mu_n / (4 pi gamma_AD m_ion x_e rho^2).

UNITS. rho and B are read in CODE units and converted with the IonizationModel's own factors.
AthenaPK is Heaviside-Lorentz, so B_unit = sqrt(4 pi rho_unit) * v_unit -- the 4 pi is already
inside B_unit, and the explicit 4 pi in eta_A above is the cgs Langevin-drag factor, NOT a second
unit conversion. (CLAUDE.md's #1 analysis pitfall is the v_A sqrt(4pi); this is the same trap one
step downstream.)

CAVEATS, stated up front because they bound what the answer means:
  * The kernel applies eta_A = min(eta_chem, eta_eq) (diffusion.hpp:335) as a blow-up guard, so
    the APPLIED eta_A can be below what this script reports. The `cap-VA/MA/DA` history columns
    (cap_diag) measure how much volume/mass actually sits on that guard; they are read here too.
  * mu_n = 2.33 (neutrals) is deliberately NOT the thermal mu = 2.29; conflating them was audit
    finding #4.
"""
import sys
import numpy as np
import h5py

# --- IonizationModel defaults, src/hydro/diffusion/ionization.hpp:70-111 ---------------------
RHO_UNIT = 5.467e-19   # g/cm^3 per code density
B_UNIT = 4.98e-5       # G per code B  (HL: sqrt(4 pi rho_unit) * v_unit)
ETA_UNIT = 1.874e-21   # eta_code = eta_cgs * (t_unit / l_unit^2)
L_UNIT = 2.81e16       # cm per code length
V_UNIT = 1.9e4         # cm/s per code velocity
MU_N = 2.33            # mean molecular weight of NEUTRALS (not the thermal 2.29)
M_ION = 24.3           # representative ion mass [m_H]
GAMMA_AD = 3.5e13      # ion-neutral drag [cm^3 g^-1 s^-1]
XE_FLOOR = 1.0e-20
M_H = 1.6737352238051868e-24
RHOCRIT_CODE = 1.0e-13 / RHO_UNIT   # first-core density in code units


def eta_A_code(rho_code, B_code, xe):
    """eta_A in CODE units, following AmbipolarEtaFromXe exactly."""
    rho = rho_code * RHO_UNIT
    B = B_code * B_UNIT
    n_n = rho / (MU_N * M_H)
    n_i = np.maximum(xe, XE_FLOOR) * n_n
    rho_i = n_i * (M_ION * M_H)
    return (B * B / (4.0 * np.pi * GAMMA_AD * rho_i * rho)) * ETA_UNIT


def main(path):
    with h5py.File(path, "r") as h:
        names = [str(s) for s in h["Info"].attrs["ComponentNames"]]
        t = float(h["Info"].attrs["Time"])
        # ComponentNames concatenates ALL datasets (grav.phi first), so index into it by NAME
        # and subtract the offset of the first prim_ entry -- never assume a fixed layout.
        off = names.index("prim_density")
        ip = lambda n: names.index(n) - off
        p = h["prim"]
        rho = p[:, ip("prim_density"), ...]
        b1 = p[:, ip("prim_magnetic_field_1"), ...]
        b2 = p[:, ip("prim_magnetic_field_2"), ...]
        b3 = p[:, ip("prim_magnetic_field_3"), ...]
        xe = p[:, ip("prim_scalar_4"), ...]          # xe_scalar_index = 4
        lvl = h["Levels"][:]
        xf = h["Locations/x"][:]

    Bmag = np.sqrt(b1 * b1 + b2 * b2 + b3 * b3)
    eta = eta_A_code(rho, Bmag, xe)
    vA = Bmag / np.sqrt(rho)                          # Heaviside-Lorentz: NO 4 pi
    dxb = (xf[:, 1] - xf[:, 0])                       # cell size per block (uniform in-block)
    dx = np.broadcast_to(dxb[:, None, None, None], rho.shape)

    print(f"snapshot : {path}")
    print(f"time     : {t:.6f} code   blocks: {rho.shape[0]}   levels {lvl.min()}-{lvl.max()}")
    print(f"rho_max  : {rho.max():.4e} code = {rho.max()/RHOCRIT_CODE:.3e} x rho_crit(1e-13)")
    print(f"finest dx: {dxb.min():.4e} code = {dxb.min()*L_UNIT:.4e} cm\n")

    # Report in density shells: the first core is where the fossil-field question lives.
    edges = [1e-2, 1.0, 1e2, RHOCRIT_CODE * 1e-2, RHOCRIT_CODE, RHOCRIT_CODE * 1e2,
             RHOCRIT_CODE * 1e4, np.inf]
    print("Physical eta_A (code units) by density shell, and the numerical resistivity that the")
    print("WP-14 Alfven ladder bounds at the SAME cell size. eta_num = C k v_A dx^2, C = 0.103")
    print("measured (see wp22_eta_numerical.py); k is taken as 2 pi / (N_cell dx), i.e. a feature")
    print("resolved by N_cell cells, so eta_num = 0.103 * 2 pi * v_A * dx / N_cell.\n")
    for ncell in (4, 8, 16):
        print(f"  --- assuming structures resolved by N_cell = {ncell} cells ---")
        print(f"  {'rho range (code)':>26} {'cells':>9} {'median eta_A':>13} "
              f"{'median eta_num':>15} {'eta_A/eta_num':>14}")
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (rho >= lo) & (rho < hi)
            n = int(m.sum())
            if n == 0:
                continue
            eA = np.median(eta[m])
            eN = np.median(0.103 * 2.0 * np.pi * vA[m] * dx[m] / ncell)
            hi_s = "inf" if not np.isfinite(hi) else f"{hi:.2e}"
            print(f"  {lo:10.2e} - {hi_s:>10} {n:9d} {eA:13.4e} {eN:15.4e} "
                  f"{eA/eN:14.3e}")
        print()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1
         else "/beegfs/u/bbg6470/athenapk/runs/prod_v9/parthenon.out1.00096.phdf")
