"""AthenaPK unit-handling for analysis scripts.

Two things this module does:

1) Read a run's cgs unit conversion factors. Source priority:
     a) units.json  (written automatically by patched collapse_be.cpp)
     b) run log     (fallback for old runs without the patch)
     c) Units.from_values()  (manual, when you know the numbers)

2) Open a Parthenon .phdf snapshot and return arrays already in cgs units,
   so analysis scripts never need to remember conversion factors.

Quick usage
-----------
    from athenapk_units import Units, read_phdf

    # Just the conversion factors:
    u = Units("/path/to/run/dir")
    print(u)                        # human-readable summary
    rho_cgs = rho_code * u.density  # g/cm^3

    # Full snapshot read, everything in cgs:
    snap = read_phdf("/path/to/run/parthenon.out0.00050.phdf")
    print(snap['t_yr'])             # time in years
    print(snap['rho'].shape)        # density, g/cm^3
    print(snap['B_microG'].shape)   # |B| in microGauss (MHD only)

    # Code units instead:
    snap = read_phdf("...", cgs=False)

    # CLI self-test / verifier:
    python3 athenapk_units.py <run_dir>
    python3 athenapk_units.py read <phdf_file>
"""
import json
import re
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Units class
# ---------------------------------------------------------------------------

class Units:
    """CGS unit conversion factors for one AthenaPK run."""

    # Physical constants in cgs
    YR       = 3.15569e7       # s
    AU       = 1.4959787e13    # cm
    PC       = 3.0857e18       # cm
    KPC      = 3.0857e21       # cm
    MSUN     = 1.9891e33       # g
    G_NEWTON = 6.67259e-8      # cm^3 g^-1 s^-2

    def __init__(self, run_dir):
        """Load units from *run_dir*.

        Looks for units.json first; falls back to parsing a run log in the
        same directory for old runs that predate the units patch.
        """
        run_dir = Path(run_dir)
        json_path = run_dir / "units.json"

        if json_path.exists():
            self._load_json(json_path)
            self.source = "units.json"
        else:
            log = self._find_log(run_dir)
            if log is None:
                raise FileNotFoundError(
                    f"No units.json found in {run_dir} and no run log found.\n"
                    "Options:\n"
                    "  1. Rebuild with the units patch and re-run (writes units.json).\n"
                    "  2. Use Units.from_values(...) to supply values manually.\n"
                    "  3. Make sure run.out / *.log is in the same directory."
                )
            self._load_log(log)
            self.source = f"{log.name} (parsed from log)"

    # ------------------------------------------------------------------
    # Private loaders
    # ------------------------------------------------------------------

    def _load_json(self, path):
        with open(path) as f:
            d = json.load(f)
        self.length   = float(d["code_length_cgs"])
        self.time     = float(d["code_time_cgs"])
        self.mass     = float(d["code_mass_cgs"])
        self.density  = float(d["code_density_cgs"])
        self.velocity = float(d.get("code_velocity_cgs",
                                    self.length / self.time))
        self.energy   = float(d.get("code_energy_cgs",
                                    self.mass * self.velocity ** 2))
        self.bfield   = float(d.get("code_bfield_cgs",
                                    (self.density ** 0.5) * self.velocity))

    def _find_log(self, run_dir):
        """Return the first log-like file that contains the normalization block."""
        for pattern in ("run.out", "*.log", "*.out", "slurm-*.out"):
            for h in sorted(run_dir.glob(pattern)):
                if h.suffix in {".phdf", ".rhst", ".athdf", ".xdmf"}:
                    continue
                try:
                    txt = h.read_text(errors="ignore")
                    if "Normalization Units" in txt:
                        return h
                except Exception:
                    continue
        return None

    def _load_log(self, log_path):
        """Parse normalization values from the run log printout."""
        text = Path(log_path).read_text(errors="ignore")

        def find(label, unit):
            # Pass the plain unit string (e.g. "g/cm^3"); re.escape handles
            # the special characters so we never double-escape.
            pattern = (
                rf"(?m)^\s*{label}\s*:\s*([0-9.eE+\-]+)\s*"
                rf"\[\s*{re.escape(unit)}\s*\]"
            )
            m = re.search(pattern, text)
            return float(m.group(1)) if m else None

        mass_msun = find("Mass",    "Msun")
        length_au = find("Length",  "au")
        time_yr   = find("Time",    "yr")
        rho0_cgs  = find("Density", "g/cm^3")   # plain string — NOT r"g/cm\^3"

        missing = [k for k, v in
                   dict(Mass=mass_msun, Length=length_au,
                        Time=time_yr, Density=rho0_cgs).items()
                   if v is None]
        if missing:
            raise ValueError(
                f"Could not parse {missing} from {log_path}.\n"
                "Check that the file contains the '--- Normalization Units ---' block."
            )

        self.length   = length_au * self.AU
        self.time     = time_yr   * self.YR
        self.mass     = mass_msun * self.MSUN
        self.density  = rho0_cgs
        self.velocity = self.length / self.time
        self.energy   = self.mass * self.velocity ** 2
        self.bfield   = (self.density ** 0.5) * self.velocity

    # ------------------------------------------------------------------
    # Alternate constructor
    # ------------------------------------------------------------------

    @classmethod
    def from_values(cls, length, time, mass, density,
                    velocity=None, energy=None, bfield=None):
        """Construct without a units.json or log file.

        Useful for old runs where you have the normalization values from
        the printed log output but don't want to point at the file.

        Example
        -------
            u = Units.from_values(
                length  = 2.806e16,   # cm
                time    = 1.477e12,   # s
                mass    = 1.208e31,   # g
                density = 5.467e-19,  # g/cm^3
            )
        """
        obj = cls.__new__(cls)
        obj.length   = float(length)
        obj.time     = float(time)
        obj.mass     = float(mass)
        obj.density  = float(density)
        obj.velocity = float(velocity) if velocity is not None else obj.length / obj.time
        obj.energy   = float(energy)   if energy   is not None else obj.mass * obj.velocity ** 2
        obj.bfield   = float(bfield)   if bfield   is not None else (obj.density ** 0.5) * obj.velocity
        obj.source   = "user-supplied"
        return obj

    # ------------------------------------------------------------------
    # Derived quantities
    # ------------------------------------------------------------------

    @property
    def pressure(self):
        """Pressure unit in erg/cm^3 = dyn/cm^2."""
        return self.density * self.velocity ** 2

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self):
        return (
            f"Units(source={self.source})\n"
            f"  length   = {self.length:.4e} cm     = {self.length/self.AU:.3f} au\n"
            f"  time     = {self.time:.4e} s      = {self.time/self.YR:.3f} yr\n"
            f"  mass     = {self.mass:.4e} g      = {self.mass/self.MSUN:.4f} Msun\n"
            f"  density  = {self.density:.4e} g/cm^3\n"
            f"  velocity = {self.velocity:.4e} cm/s\n"
            f"  pressure = {self.pressure:.4e} erg/cm^3\n"
            f"  bfield   = {self.bfield:.4e} G    = {self.bfield*1e6:.3f} uG"
        )


# ---------------------------------------------------------------------------
# Snapshot reader
# ---------------------------------------------------------------------------

def _resolve_components(component_names, prim_size):
    """Map AthenaPK prim ComponentNames → index in the prim array.

    ComponentNames may include entries for separate datasets (e.g. grav.phi)
    that are NOT columns in the prim array. We only keep names that start with
    'prim_' and stop once we've assigned prim_size indices.
    """
    idx = {}
    prim_counter = 0
    for raw in component_names:
        name = raw.decode() if isinstance(raw, bytes) else str(raw)
        if name.startswith("prim_"):
            idx[name] = prim_counter
            prim_counter += 1
            if prim_counter >= prim_size:
                break
    return idx


def read_phdf(filename, cgs=True):
    """Read a Parthenon .phdf snapshot into a dict of numpy arrays.

    Parameters
    ----------
    filename : str or Path
        Path to the .phdf file.
    cgs : bool
        True  → return all arrays in cgs (g, cm, s, erg, Gauss).
        False → return code-unit arrays unchanged.

    Returns
    -------
    dict with keys:

        t          float    simulation time (s if cgs, code-time if not)
        t_yr       float    simulation time in years
        nblocks    int      number of MeshBlocks
        levels     ndarray  (nblocks,) AMR level per block
        x1f, x2f, x3f  ndarray  face coordinates per block (cm or code)
        rho        ndarray  (nblocks, nz, ny, nx) density
        vx, vy, vz ndarray  velocity components
        |v|        ndarray  velocity magnitude
        P          ndarray  pressure
        Bx, By, Bz ndarray  B-field components  [MHD only]
        |B|        ndarray  |B| magnitude        [MHD only]
        B_microG   ndarray  |B| in microGauss    [MHD + cgs only]
        psi        ndarray  GLM scalar           [MHD only, may be None]
        phi_grav   ndarray  gravitational potential (erg/g or code)
        mhd        bool     True if B-field variables are present
        units      Units    the Units object used for conversion
    """
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py is required for read_phdf. Install with: pip install h5py")

    filename = Path(filename)
    run_dir  = filename.parent

    # Load units (required for cgs mode, optional otherwise)
    units = None
    if cgs:
        units = Units(run_dir)
    else:
        try:
            units = Units(run_dir)
        except Exception:
            pass

    out = {"units": units}

    with h5py.File(filename, "r") as f:
        info      = f["Info"].attrs
        prim_size = f["prim"].shape[1]
        comp_idx  = _resolve_components(info["ComponentNames"], prim_size)
        t_code    = float(info["Time"])

        # --- Time ---
        if cgs and units:
            out["t"]    = t_code * units.time
            out["t_yr"] = out["t"] / Units.YR
        else:
            out["t"]    = t_code
            out["t_yr"] = t_code / Units.YR if units else t_code

        # --- Block metadata ---
        out["nblocks"] = int(info["NumMeshBlocks"])
        out["levels"]  = np.asarray(f["Levels"][:])

        # --- Face coordinates (cell edges) ---
        scale_x = units.length if (cgs and units) else 1.0
        for ax in ("x1f", "x2f", "x3f"):
            if ax in f:
                out[ax] = np.asarray(f[ax][:]) * scale_x

        # --- Primitive variables ---
        # prim shape: (nblocks, ncomp, nz, ny, nx)
        prim = f["prim"][:]

        def get(name):
            return prim[:, comp_idx[name], ...]

        scale_v = units.velocity if (cgs and units) else 1.0
        scale_P = units.pressure if (cgs and units) else 1.0
        scale_r = units.density  if (cgs and units) else 1.0

        out["rho"] = get("prim_density")    * scale_r
        out["vx"]  = get("prim_velocity_1") * scale_v
        out["vy"]  = get("prim_velocity_2") * scale_v
        out["vz"]  = get("prim_velocity_3") * scale_v
        out["P"]   = get("prim_pressure")   * scale_P
        out["|v|"] = np.sqrt(out["vx"]**2 + out["vy"]**2 + out["vz"]**2)

        # --- MHD ---
        is_mhd = "prim_magnetic_field_1" in comp_idx
        out["mhd"] = is_mhd
        if is_mhd:
            scale_B = units.bfield if (cgs and units) else 1.0
            out["Bx"]  = get("prim_magnetic_field_1") * scale_B
            out["By"]  = get("prim_magnetic_field_2") * scale_B
            out["Bz"]  = get("prim_magnetic_field_3") * scale_B
            out["|B|"] = np.sqrt(out["Bx"]**2 + out["By"]**2 + out["Bz"]**2)
            if cgs:
                out["B_microG"] = out["|B|"] * 1e6
            # psi (GLM scalar) is sometimes a separate dataset, not in prim
            psi_in_prim = (
                "prim_magnetic_psi" in comp_idx
                and comp_idx["prim_magnetic_psi"] < prim.shape[1]
            )
            if psi_in_prim:
                out["psi"] = get("prim_magnetic_psi") * scale_B
            elif "prim_magnetic_psi" in f:
                out["psi"] = np.asarray(f["prim_magnetic_psi"][:]) * scale_B
            else:
                out["psi"] = None

        # --- Gravitational potential ---
        if "grav.phi" in f:
            phi = np.asarray(f["grav.phi"][:])
            if phi.ndim == 5 and phi.shape[1] == 1:
                phi = phi[:, 0, ...]          # drop singleton component axis
            scale_phi = (units.velocity ** 2) if (cgs and units) else 1.0
            out["phi_grav"] = phi * scale_phi

    return out


# ---------------------------------------------------------------------------
# CLI  (verify units  OR  demo snapshot read)
# ---------------------------------------------------------------------------

def _verify(run_dir):
    """Print Units and cross-check against the run log."""
    run_dir = Path(run_dir)
    u = Units(run_dir)
    print(u)
    print()

    log_files = (
        sorted(run_dir.glob("run.out"))
        or sorted(run_dir.glob("*.log"))
        or sorted(run_dir.glob("*.out"))
    )
    if not log_files:
        print("No log file found — skipping cross-check.")
        return 0

    log_path = log_files[0]
    text     = log_path.read_text(errors="ignore")

    def find(label, unit):
        pattern = (
            rf"(?m)^\s*{label}\s*:\s*([0-9.eE+\-]+)\s*"
            rf"\[\s*{re.escape(unit)}\s*\]"
        )
        m = re.search(pattern, text)
        return float(m.group(1)) if m else None

    log_vals  = {
        "mass [Msun]":   find("Mass",    "Msun"),
        "length [au]":   find("Length",  "au"),
        "time [yr]":     find("Time",    "yr"),
        "density [cgs]": find("Density", "g/cm^3"),   # plain string
    }
    json_vals = {
        "mass [Msun]":   u.mass    / u.MSUN,
        "length [au]":   u.length  / u.AU,
        "time [yr]":     u.time    / u.YR,
        "density [cgs]": u.density,
    }

    print(f"Cross-check vs {log_path.name}:")
    fail = False
    for k in log_vals:
        lv, jv = log_vals[k], json_vals[k]
        if lv is None:
            print(f"  {k:18s}  NOT FOUND IN LOG")
            continue
        rel    = abs(lv - jv) / max(abs(lv), 1e-30)
        status = "OK" if rel < 1e-3 else "MISMATCH"
        if rel >= 1e-3:
            fail = True
        print(f"  {k:18s}  log={lv:.6e}  json={jv:.6e}  rel_diff={rel:.2e}  [{status}]")

    print()
    print("RESULT: All checks passed." if not fail else "RESULT: MISMATCH — investigate above.")
    return 1 if fail else 0


def _demo_read(phdf_file):
    snap = read_phdf(phdf_file)
    u    = snap["units"]
    print(u)
    print()
    print(f"Snapshot:  {Path(phdf_file).name}")
    print(f"  Time:    {snap['t_yr']:.4f} yr  ({snap['t']:.4e} s)")
    print(f"  Blocks:  {snap['nblocks']}  levels {snap['levels'].min()}–{snap['levels'].max()}")
    print(f"  rho:     min={snap['rho'].min():.3e}  max={snap['rho'].max():.3e}  g/cm^3")
    print(f"  |v|:     min={snap['|v|'].min():.3e}  max={snap['|v|'].max():.3e}  cm/s")
    if snap["mhd"]:
        print(f"  |B|:     min={snap['B_microG'].min():.3f}  max={snap['B_microG'].max():.3f}  uG")
    if "phi_grav" in snap:
        print(f"  phi:     min={snap['phi_grav'].min():.3e}  max={snap['phi_grav'].max():.3e}  erg/g")
    if "x1f" in snap:
        L = snap["x1f"][0, -1] - snap["x1f"][0, 0]
        print(f"  box:     {L/u.AU:.1f} au per block (x1)")


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python3 athenapk_units.py <run_dir>         # verify units")
        print("  python3 athenapk_units.py read <phdf_file>  # snapshot read demo")
        sys.exit(2)
    if sys.argv[1] == "read":
        if len(sys.argv) < 3:
            print("Usage: python3 athenapk_units.py read <phdf_file>")
            sys.exit(2)
        _demo_read(sys.argv[2])
        sys.exit(0)
    sys.exit(_verify(sys.argv[1]))
