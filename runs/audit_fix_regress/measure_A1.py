#!/usr/bin/env python
"""
AUDIT A1 -- measure, before fixing.

refinement/jeans.cpp:54 and refinement/nonideal.cpp:59 compute the sound speed for the
Truelove criterion as sqrt(gamma*p/rho) with the IDEAL gamma, while the production decks
run <hydro> eos = hydrogen (tabulated multi-Saha). The dynamics uses the tabulated
sound speed (AdiabaticGLMMHDEOS::SoundSpeed -> eos_tab.AsqFromRhoPres); the refinement
criterion does not.

This reproduces BOTH criteria cell by cell on a real production snapshot and reports:
  (1) the c_s ratio distribution,
  (2) how many BLOCKS change their AMR tag -- the only number that actually matters,
      because the criterion is a per-block min-reduction, not a per-cell test.

Run:  python measure_A1.py [snapshot.phdf]
"""
import sys, os
import numpy as np
import h5py

SNAP = sys.argv[1] if len(sys.argv) > 1 else \
    "/beegfs/u/bbg6470/athenapk/runs/prod_v9/parthenon.out1.final.phdf"
EOS_BIN = "/beegfs/u/bbg6470/athenapk/src/eos/eos_table_hires.bin"
GAMMA = 1.4          # <hydro> gamma in the production decks
NJEANS = 8.0         # <refinement> njeans


# --------------------------------------------------------------------------- EOS table
def read_eos(path):
    with open(path, "rb") as f:
        nr, ne, nT = np.fromfile(f, dtype=np.int64, count=3)
        lr0, dlr, le0, dle, lT0, dlT = np.fromfile(f, dtype=np.float64, count=6)
        n = int(nr) * int(ne)
        P = np.fromfile(f, dtype=np.float64, count=n).reshape(nr, ne)
        cs2 = np.fromfile(f, dtype=np.float64, count=n).reshape(nr, ne)
        logT = np.fromfile(f, dtype=np.float64, count=n).reshape(nr, ne)
    return dict(nr=int(nr), ne=int(ne), lr0=lr0, dlr=dlr, le0=le0, dle=dle,
                P=P, cs2=cs2, logT=logT)


def bilin(A, lr0, dlr, le0, dle, x, y):
    """Same edge-clamped bilinear interpolation the device does (eos_table.hpp::bilin)."""
    n1, n2 = A.shape
    fi = (x - lr0) / dlr
    i = np.clip(fi.astype(np.int64), 0, n1 - 2)
    ti = np.clip(fi - i, 0.0, 1.0)
    fj = (y - le0) / dle
    j = np.clip(fj.astype(np.int64), 0, n2 - 2)
    tj = np.clip(fj - j, 0.0, 1.0)
    return ((1 - ti) * (1 - tj) * A[i, j] + ti * (1 - tj) * A[i + 1, j] +
            (1 - ti) * tj * A[i, j + 1] + ti * tj * A[i + 1, j + 1])


def asq_from_rho_pres(E, rho, pres):
    """EosTable::AsqFromRhoPres = Cs2of(rho, EspFromP(rho, pres)).

    EspFromP is a bisection on the monotone esp axis of P(rho, esp); vectorised here as
    a bisection over the same log-esp grid, which is what the device does."""
    lr = np.log10(np.maximum(rho, 1e-300))
    lo = np.full(rho.shape, E["le0"])
    hi = np.full(rho.shape, E["le0"] + E["dle"] * (E["ne"] - 1))
    for _ in range(60):                       # >> device iteration count; converged
        mid = 0.5 * (lo + hi)
        Pm = bilin(E["P"], E["lr0"], E["dlr"], E["le0"], E["dle"], lr, mid)
        too_low = Pm < pres
        lo = np.where(too_low, mid, lo)
        hi = np.where(too_low, hi, mid)
    le = 0.5 * (lo + hi)
    return bilin(E["cs2"], E["lr0"], E["dlr"], E["le0"], E["dle"], lr, le), le


def temperature_K(E, rho, le):
    lr = np.log10(np.maximum(rho, 1e-300))
    return 10.0 ** bilin(E["logT"], E["lr0"], E["dlr"], E["le0"], E["dle"], lr, le)


# --------------------------------------------------------------------------- main
def main():
    E = read_eos(EOS_BIN)
    print("snapshot :", SNAP)
    with h5py.File(SNAP, "r") as g:
        names = [n.decode() if isinstance(n, bytes) else n
                 for n in g["Info"].attrs["ComponentNames"]]
        i0 = names.index("prim_density")           # prim block starts here
        prim = g["prim"]
        nb = prim.shape[0]
        levels = g["Levels"][:]
        t = float(g["Info"].attrs["Time"])
        # cell size per block: from the mesh extent and the level
        try:
            xmin = float(g["Info"].attrs["RootGridDomain"][0])
            xmax = float(g["Info"].attrs["RootGridDomain"][1])
            nx = int(g["Info"].attrs["RootGridSize"][0])
        except Exception:
            xmin, xmax, nx = -8.0, 8.0, 256
        dx_root = (xmax - xmin) / nx
        print("t = %.6f   blocks = %d   levels %d..%d   dx_root = %.4g"
              % (t, nb, levels.min(), levels.max(), dx_root))

        nj_ideal = np.empty(nb)
        nj_table = np.empty(nb)
        cs_ratio_all, T_all, w_all = [], [], []

        # component offsets inside `prim` (density is component index 0 of the prim block)
        c_rho, c_pre = 0, 4
        c_b1, c_b2, c_b3 = 5, 6, 7

        for b in range(nb):
            blk = prim[b]
            rho = blk[c_rho].astype(np.float64)
            pre = blk[c_pre].astype(np.float64)
            bsq = (blk[c_b1].astype(np.float64) ** 2 +
                   blk[c_b2].astype(np.float64) ** 2 +
                   blk[c_b3].astype(np.float64) ** 2)
            good = (rho > 0) & (pre > 0)
            if not good.any():
                nj_ideal[b] = nj_table[b] = np.inf
                continue
            r, p, bb = rho[good], pre[good], bsq[good]

            cs_i = np.sqrt(GAMMA * p / r)                       # what the code does now
            a2, le = asq_from_rho_pres(E, r, p)                 # what the dynamics uses
            cs_t = np.sqrt(np.maximum(a2, 1e-300))
            va = np.sqrt(bb / r)

            dx = dx_root / (2.0 ** int(levels[b]))
            fac = 2.0 * np.pi / dx
            nj_ideal[b] = fac * np.min((cs_i + va) / np.sqrt(r))
            nj_table[b] = fac * np.min((cs_t + va) / np.sqrt(r))

            cs_ratio_all.append(cs_i / cs_t)
            T_all.append(temperature_K(E, r, le))
            w_all.append(r)

        ratio = np.concatenate(cs_ratio_all)
        Tk = np.concatenate(T_all)
        rhow = np.concatenate(w_all)

    print("\n--- (1) sound-speed ratio  c_s(ideal gamma) / c_s(EOS table), per cell ---")
    for q in (50, 90, 99, 99.9):
        print("   p%-5s %.4f" % (q, np.percentile(ratio, q)))
    print("   max   %.4f   min %.4f   cells %d" % (ratio.max(), ratio.min(), ratio.size))
    hot = Tk > 1500.0
    print("   cells above 1500 K (dissociation onset): %d (%.3g %%)"
          % (hot.sum(), 100.0 * hot.mean()))
    if hot.any():
        print("   ratio there: median %.4f  p99 %.4f  max %.4f"
              % (np.median(ratio[hot]), np.percentile(ratio[hot], 99), ratio[hot].max()))
    print("   max gas temperature in the snapshot: %.4g K" % Tk.max())

    def tag(nj):
        t = np.zeros(nj.shape, dtype=np.int8)      # 0 = same
        t[nj < NJEANS] = 1                          # refine
        t[nj > 2.5 * NJEANS] = -1                   # derefine
        return t

    ti, tt = tag(nj_ideal), tag(nj_table)
    flip = ti != tt
    print("\n--- (2) AMR tag, per block (njeans = %.0f) ---" % NJEANS)
    print("   blocks: %d    tags CHANGED by the fix: %d (%.3g %%)"
          % (len(ti), flip.sum(), 100.0 * flip.mean()))
    for a, bl in ((1, "refine"), (0, "same"), (-1, "derefine")):
        for c, blc in ((1, "refine"), (0, "same"), (-1, "derefine")):
            n = int(((ti == a) & (tt == c)).sum())
            if n and a != c:
                print("       %-8s -> %-8s : %d blocks" % (bl, blc, n))
    print("   totals  ideal: refine %d same %d derefine %d"
          % ((ti == 1).sum(), (ti == 0).sum(), (ti == -1).sum()))
    print("   totals  table: refine %d same %d derefine %d"
          % ((tt == 1).sum(), (tt == 0).sum(), (tt == -1).sum()))

    print("\n--- (3) falsification check ---")
    print("   The fix must only ever refine MORE, never less.")
    worse = ((ti == 1) & (tt != 1)).sum() + ((ti == 0) & (tt == -1)).sum()
    print("   blocks where the fix refines LESS: %d  %s"
          % (worse, "<-- VIOLATION" if worse else "(none, as required)"))
    print("   nj_table <= nj_ideal in every block: %s"
          % bool(np.all(nj_table <= nj_ideal * (1 + 1e-12))))


if __name__ == "__main__":
    main()
