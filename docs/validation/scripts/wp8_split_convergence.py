#!/usr/bin/env python
"""WP-8 completion — does the DENSITY SPLIT actually restore convergence?

WHAT WAS ALREADY ESTABLISHED (WP08_dissipation_nonconvergence.md, 2026-07-31).
`mag-dissO` / `mag-dissA` / `mag-Jsq` do not converge on the njeans ladder while `mag-ME`
converges to 0.2 %. Root cause: they are not volume integrals in any useful sense. The
integrand spans ~7 decades and 90 % of it is carried by ~1e-7 of the box, so each refinement
resolves a previously-unresolved region and changes the answer. Two competing explanations
(diagnostic stencil; stale eta cache) were falsified. The proposed remedy was to split the
budget at a density threshold -- core and envelope are then regions defined by PHYSICS, not by
the grid, and each can converge (or visibly fail to) on its own. Columns `mag-dissO-hi/lo`,
`mag-dissA-hi/lo`, `mag-Vhi` and the inverse-participation-ratio f_eff were implemented behind
`hydro/mag_diag_rho_split` (default 0 = OFF, so the OFF state stays byte-identical).

WHAT WAS NEVER DONE, AND IS WHAT THIS SCRIPT DOES.
**The split was never applied to the ladder.** The non-convergence table was measured with the
old GLOBAL columns, and the split was verified only for self-consistency on a 12-cycle L=52
smoke deck (hi + lo == global to the history file's text precision). So the remedy's PREMISE --
that splitting restores convergence -- is untested. A remedy that is implemented, gated and
documented but never shown to work is not a closed finding.

WHY Jsq AND NOT dissO/dissA. Re-running the ladder with the split enabled costs GPU time this
project does not have spare. Everything needed is already on disk: 135 snapshots, 986 GB. But
the ladder snapshots carry only `prim` -- there is no `nonideal_eta` field (that was added for
WP-22 part 3, job 2448571, on a different run), and reconstructing eta_A offline is only an
UPPER bound because the equilibrium NICIL/Wardle ceiling needs the full grain + Saha charge
solve (wp22_eta_physical.py docstring). Rather than test the remedy through a quantity we can
only bound, test it on **Jsq = integral |J|^2 dV**, which needs only B and the grid and is
therefore EXACT offline -- and which is itself one of the non-converging quantities
(-50.0 %, -39.6 %). If Jsq-hi and Jsq-lo converge while global Jsq does not, the split's premise
is demonstrated on a quantity nobody has to trust me about. If they do not, the remedy does not
work and WP-8 is not closed by it.

FALSIFIER, and it matters. A split can look convergent for a trivial reason: if essentially all
of the integral lands in ONE bin, that bin simply reproduces the global number and the other is
noise. So this script reports, for every leg and bin, the fraction of the global integral the
bin carries AND its f_eff. A "converged" bin holding 99.9 % of the integral with f_eff ~ 1e-7 is
the original pathology wearing a different label, not a fix.

METHOD. J = curl B by the 1-cell face difference -- the stencil WP-8 verified against the 2*dx
centred one (they agreed to 1.5 % and both gave identical non-convergence, so this choice is not
load-bearing). Differences are taken INSIDE each block only: phdf carries no ghost zones, so the
outermost cell layer of every block is dropped rather than differenced across a block boundary
with wrong data. That is a real approximation and is reported as the fraction of cells dropped.

Legs are compared at MATCHED DENSITY, not matched time (the three legs reach a given rho_max at
different t). The snapshot nearest the target in log rho is used; the residual epoch mismatch is
printed, because it -- not the split -- is the leading systematic if it is large.

USAGE
    wp8_split_convergence.py [rho_target_cgs]      default 1e-12 (the epoch of the WP-8 table)
"""
import sys
import glob
import numpy as np
import h5py

RHO0 = 5.467e-19                      # code density unit [g/cm^3]
LEGS = ["nj4", "nj8", "nj16"]
LADDER = "/beegfs/u/bbg6470/athenapk/runs/convergence_ladder"
# Split at rho_crit = 1e-13 g/cm^3, the first-core threshold: above it is the optically thick
# core the fossil-field question lives in, below it the envelope. This is the same physical
# boundary the code's mag_diag_rho_split is meant to be set to, not a tuned number.
RHO_SPLIT_CODE = 1.0e-13 / RHO0


def curl_sq_interior(b1, b2, b3, dx):
    """|curl B|^2 on block interiors. Arrays are [block, k, j, i]; dx is per-block, uniform."""
    d = dx[:, None, None, None]
    # central differences, interior only (drop one layer each side on every axis)
    dB3_dy = (b3[:, 1:-1, 2:, 1:-1] - b3[:, 1:-1, :-2, 1:-1]) / (2 * d)
    dB2_dz = (b2[:, 2:, 1:-1, 1:-1] - b2[:, :-2, 1:-1, 1:-1]) / (2 * d)
    dB1_dz = (b1[:, 2:, 1:-1, 1:-1] - b1[:, :-2, 1:-1, 1:-1]) / (2 * d)
    dB3_dx = (b3[:, 1:-1, 1:-1, 2:] - b3[:, 1:-1, 1:-1, :-2]) / (2 * d)
    dB2_dx = (b2[:, 1:-1, 1:-1, 2:] - b2[:, 1:-1, 1:-1, :-2]) / (2 * d)
    dB1_dy = (b1[:, 1:-1, 2:, 1:-1] - b1[:, 1:-1, :-2, 1:-1]) / (2 * d)
    jx = dB3_dy - dB2_dz
    jy = dB1_dz - dB3_dx
    jz = dB2_dx - dB1_dy
    return jx * jx + jy * jy + jz * jz


def f_eff(q, dV, Vtot):
    """Inverse participation ratio: fraction of the box that would carry the whole integral
    if q were uniform on it. ~1 volume-filling; ~1e-7 a point sample."""
    I1 = float((q * dV).sum())
    I2 = float((q * q * dV).sum())
    if I2 <= 0.0:
        return float("nan"), I1
    return I1 * I1 / (Vtot * I2), I1


def scan(leg, rho_target_cgs):
    """Nearest snapshot in log rho_max, read from the PRECOMPUTED epoch_scan.txt.

    Do NOT reopen every snapshot to find rho_max. The first version of this function did, and it
    reads the full density array of all 135 phdf files -- ~900 GB -- to recover three numbers that
    `epoch_scan.txt` already tabulates (it exists for exactly this, written by the WP-8 epoch
    sweep). Killed after 13 minutes without finishing a single leg.
    """
    best = None
    for L in open(f"{LADDER}/epoch_scan.txt"):
        if L.startswith("#") or "ERR" in L:
            continue
        p = L.split()
        if len(p) < 6 or p[0] != leg:
            continue
        f, t, rmax = f"{LADDER}/{leg}/{p[1]}", float(p[2]), float(p[5])
        if rmax <= 0:
            continue
        d = abs(np.log10(rmax) - np.log10(rho_target_cgs))
        if best is None or d < best[0]:
            best = (d, f, rmax / RHO0, t)   # rho_max returned in CODE units
    return best


def measure(path):
    with h5py.File(path, "r") as h:
        names = [str(s) for s in h["Info"].attrs["ComponentNames"]]
        off = names.index("prim_density")
        ip = lambda n: names.index(n) - off
        p = h["prim"]
        rho = p[:, ip("prim_density"), ...]
        b1 = p[:, ip("prim_magnetic_field_1"), ...]
        b2 = p[:, ip("prim_magnetic_field_2"), ...]
        b3 = p[:, ip("prim_magnetic_field_3"), ...]
        xf = h["Locations/x"][:]
    dxb = xf[:, 1] - xf[:, 0]
    jsq = curl_sq_interior(b1, b2, b3, dxb)
    rho_i = rho[:, 1:-1, 1:-1, 1:-1]
    dV = np.broadcast_to((dxb ** 3)[:, None, None, None], jsq.shape)
    Vtot = float(dV.sum())
    dropped = 1.0 - jsq.size / rho.size

    hi = rho_i > RHO_SPLIT_CODE
    out = {}
    out["V"] = Vtot
    out["dropped"] = dropped
    out["Vhi_frac"] = float(dV[hi].sum()) / Vtot
    fe, I = f_eff(jsq, dV, Vtot)
    out["glob"] = (I, fe)
    for lab, m in (("hi", hi), ("lo", ~hi)):
        if m.sum() == 0:
            out[lab] = (0.0, float("nan"))
            continue
        fe_b, I_b = f_eff(jsq[m], dV[m], float(dV[m].sum()))
        out[lab] = (I_b, fe_b)
    return out


def main():
    rho_t_cgs = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-12
    rho_t = rho_t_cgs / RHO0
    print(f"WP-8 split convergence test — matched epoch rho_max = {rho_t_cgs:.3e} g/cm^3 "
          f"({rho_t:.4e} code)")
    print(f"density split = {1e-13:.0e} g/cm^3 ({RHO_SPLIT_CODE:.4e} code) = rho_crit\n")

    res = {}
    for leg in LEGS:
        b = scan(leg, rho_t_cgs)
        if b is None:
            print(f"{leg}: no usable snapshot"); continue
        d, f, rmax, t = b
        m = measure(f)
        res[leg] = m
        print(f"{leg:>5}  {f.split('/')[-1]:<28} t={t:.6f}  rho_max={rmax*RHO0:.4e} g/cm^3 "
              f"(epoch miss {10**d:.2f}x)")
        print(f"        interior cells kept {100*(1-m['dropped']):.1f}%   "
              f"V(rho>split)/V = {m['Vhi_frac']:.3e}")
        for lab in ("glob", "hi", "lo"):
            I, fe = m[lab]
            frac = I / m["glob"][0] if m["glob"][0] else float("nan")
            print(f"        Jsq-{lab:<4} = {I: .6e}   f_eff = {fe:.3e}   "
                  f"carries {100*frac:6.2f}% of global")
        print()

    print("=" * 78)
    print("CONVERGENCE (successive relative change; the WP-8 table's own metric)")
    print("=" * 78)
    print(f"{'quantity':<12} {'nj4->nj8':>12} {'nj8->nj16':>12}   verdict")
    for lab in ("glob", "hi", "lo"):
        try:
            a, b_, c = (res[L][lab][0] for L in LEGS)
        except KeyError:
            continue
        r1 = (b_ - a) / a if a else float("nan")
        r2 = (c - b_) / b_ if b_ else float("nan")
        # ORDER MATTERS. Check monotonicity FIRST. An earlier version tested "did the magnitude
        # shrink" first, so a SIGN FLIP with a smaller second step (-80.9% -> +25.6%) was reported
        # as CONVERGING -- which is precisely the non-monotone pathology WP-8 exists to name. With
        # three points and a sign change you cannot claim convergence at all, however small the
        # second step is.
        monot = np.sign(r1) == np.sign(r2)
        shrink = abs(r2) < abs(r1) / 2
        if not monot:
            v = "NOT MONOTONE - no convergence claim possible"
        elif shrink:
            v = "CONVERGING"
        else:
            v = "monotone, not converged"
        print(f"Jsq-{lab:<8} {100*r1:11.1f}% {100*r2:11.1f}%   {v}")
    print()
    print("Read the f_eff column before believing any 'CONVERGING' verdict: a bin that carries")
    print("~100% of the integral with f_eff ~ 1e-7 is the original pathology relabelled, not a")
    print("fix. The split only means something if BOTH bins carry a real share of the budget.")


if __name__ == "__main__":
    main()
