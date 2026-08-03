#!/usr/bin/env python
"""WP-22 part 3 — the APPLIED ambipolar diffusivity, read straight out of the code.

Part 2 (`wp22_eta_physical.py`) reconstructed eta_A offline from (rho, B, x_e). That is only the
CHEMISTRY branch; what the flux kernel applies is

    eta_A = min(eta_chem, eta_eq, eta_ad_cap)            [diffusion.hpp:310-341]

with eta_ad_cap inactive in production (defaults to numeric_limits<Real>::max(), hydro.cpp:1085)
but eta_eq -- the equilibrium NICIL/Wardle ceiling -- active, and NOT reproducible offline without
porting the whole grain + Saha charge solve. So part 2's number is an UPPER BOUND, which is the
wrong direction for a claim that physical eta DOMINATES numerical eta.

This script removes that caveat. The code already caches the applied (eta_O, eta_H, eta_A) per cell
in the Parthenon field `nonideal_eta` (hydro.cpp:1391-1395) -- it is simply never written. Job
2448571 restarted prod_v9 with that field added to the output list, so the applied value is now a
phdf field and can be compared directly. The chem-branch value is recomputed alongside it so the
suppression factor eta_chem/eta_applied is measured rather than assumed.
"""
import sys
import numpy as np
import h5py

sys.path.insert(0, "/beegfs/u/bbg6470/athenapk/docs/validation/scripts")
from wp22_eta_physical import eta_A_code, RHOCRIT_CODE, L_UNIT  # noqa: E402

C_NUM = 0.103  # eta_num = C k v_A dx^2, measured on the WP-14 Alfven ladder (+-3% over 8x in dx)


def main(path):
    with h5py.File(path, "r") as h:
        names = [str(s) for s in h["Info"].attrs["ComponentNames"]]
        t = float(h["Info"].attrs["Time"])
        off = names.index("prim_density")
        ip = lambda n: names.index(n) - off
        p = h["prim"]
        rho = p[:, ip("prim_density"), ...]
        b1 = p[:, ip("prim_magnetic_field_1"), ...]
        b2 = p[:, ip("prim_magnetic_field_2"), ...]
        b3 = p[:, ip("prim_magnetic_field_3"), ...]
        xe = p[:, ip("prim_scalar_4"), ...]
        # nonideal_eta is its own dataset; components are (eta_ohmic, eta_hall, eta_ambipolar).
        eta_app = h["nonideal_eta"][:, 2, ...]
        eta_ohm = h["nonideal_eta"][:, 0, ...]
        xf = h["Locations/x"][:]

    Bmag = np.sqrt(b1 * b1 + b2 * b2 + b3 * b3)
    eta_chem = eta_A_code(rho, Bmag, xe)
    vA = Bmag / np.sqrt(rho)                       # Heaviside-Lorentz, no 4 pi
    dx = np.broadcast_to((xf[:, 1] - xf[:, 0])[:, None, None, None], rho.shape)

    print(f"snapshot : {path}")
    print(f"time     : {t:.6f}   blocks {rho.shape[0]}   finest dx {dx.min():.4e} code "
          f"= {dx.min()*L_UNIT:.3e} cm")
    print(f"rho_max  : {rho.max():.4e} code = {rho.max()/RHOCRIT_CODE:.3e} x rho_crit\n")

    edges = [1e-2, 1.0, 1e2, RHOCRIT_CODE * 1e-2, RHOCRIT_CODE, RHOCRIT_CODE * 1e2,
             RHOCRIT_CODE * 1e4, np.inf]
    for ncell in (4, 8):
        print(f"  === structures resolved by N_cell = {ncell} cells "
              f"(eta_num = {C_NUM} * 2pi * v_A * dx / N_cell) ===")
        print(f"  {'rho range (code)':>26} {'cells':>9} {'eta_APPLIED':>12} {'eta_chem':>12} "
              f"{'suppress':>9} {'eta_num':>12} {'APPLIED/num':>12}")
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (rho >= lo) & (rho < hi)
            n = int(m.sum())
            if n == 0:
                continue
            eA = np.median(eta_app[m])
            eC = np.median(eta_chem[m])
            eN = np.median(C_NUM * 2.0 * np.pi * vA[m] * dx[m] / ncell)
            hi_s = "inf" if not np.isfinite(hi) else f"{hi:.2e}"
            sup = eC / eA if eA > 0 else np.inf
            print(f"  {lo:10.2e} - {hi_s:>10} {n:9d} {eA:12.4e} {eC:12.4e} {sup:9.2f}x "
                  f"{eN:12.4e} {eA/eN if eN > 0 else np.inf:12.3e}")
        print()

    # Ohmic is separately capped at eta_ohm_cap_code = 0.1 in the production deck; report how much
    # of the domain sits exactly on that ceiling, since a capped cell returns the ceiling exactly.
    on_cap = np.isclose(eta_ohm, 0.1, rtol=1e-12)
    print(f"Ohmic cap (eta_ohm_cap_code = 0.1): {on_cap.sum()} of {eta_ohm.size} cells "
          f"({100.0*on_cap.mean():.3f} %) sit exactly on the ceiling.")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1
         else "/beegfs/u/bbg6470/athenapk/runs/wp22_eta/dump/parthenon.out1.00096.phdf")
