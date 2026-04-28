"""Verify the units patch works end-to-end.

Reads units.json from a run dir and cross-checks the values against the
"---  Normalization Units  ---" printout that collapse_be.cpp emits at
startup. Catches mismatches between code constants and what gets written.

Usage:
    python3 verify_units.py /path/to/run/dir
"""
import sys
import re
from pathlib import Path
from athenapk_units import Units


def parse_log(log_text):
    """Pull normalization values from the run log.

    Anchors to the '--- Normalization Units ---' block so it never
    accidentally picks up 'Total mass' or 'Central density' lines
    from the dimensional-parameters block above.
    """
    # Extract just the normalization block (between its header and the next ---)
    norm_match = re.search(
        r"---\s*Normalization Units\s*---(.+?)(?:---|$)",
        log_text,
        re.DOTALL | re.IGNORECASE,
    )
    if not norm_match:
        return {}

    block = norm_match.group(1)

    patterns = {
        "mass_msun":   r"Mass\s*:\s*([0-9.eE+\-]+)\s*\[Msun\]",
        "length_au":   r"Length\s*:\s*([0-9.eE+\-]+)\s*\[au\]",
        "time_yr":     r"Time\s*:\s*([0-9.eE+\-]+)\s*\[yr\]",
        "density_cgs": r"Density\s*:\s*([0-9.eE+\-]+)\s*\[g/cm\^3\]",
    }

    out = {}
    for key, pat in patterns.items():
        m = re.search(pat, block, re.IGNORECASE)
        if m:
            out[key] = float(m.group(1))

    return out


def main(run_dir):
    run_dir = Path(run_dir)

    # 1. Load and print units.json
    u = Units(run_dir)
    print("units.json contents:")
    print(u)
    print()

    # 2. Find log file
    log_files = sorted(
        list(run_dir.glob("*.out")) + list(run_dir.glob("*.log")),
        key=lambda p: p.stat().st_mtime,
        reverse=True,   # newest first
    )
    if not log_files:
        print("No .out / .log files found -- skipping log cross-check.")
        return 0

    log_file = log_files[0]
    log_text = log_file.read_text(errors="replace")  # handles any encoding

    # 3. Parse normalization block
    parsed = parse_log(log_text)
    if not parsed:
        print(f"WARNING: Could not find '--- Normalization Units ---' block in {log_file.name}.")
        print("  Cross-check skipped. Is this the right log file?")
        print(f"  Searched: {log_file}")
        return 0

    # 4. Cross-check
    print(f"Cross-check vs {log_file.name}:")
    checks = [
        ("mass [Msun]",   parsed.get("mass_msun"),   u.mass    / u.MSUN),
        ("length [au]",   parsed.get("length_au"),   u.length  / u.AU),
        ("time [yr]",     parsed.get("time_yr"),     u.time    / u.YR),
        ("density [cgs]", parsed.get("density_cgs"), u.density),
    ]
    fail = False
    for label, log_val, json_val in checks:
        if log_val is None:
            print(f"  {label:18s}  (not found in log)")
            continue
        rel = abs(log_val - json_val) / max(abs(log_val), 1e-30)
        status = "OK" if rel < 1e-4 else "MISMATCH"
        if rel >= 1e-4:
            fail = True
        print(
            f"  {label:18s}  log={log_val:.6e}  json={json_val:.6e}  "
            f"rel_diff={rel:.2e}  [{status}]"
        )

    print()
    if fail:
        print("RESULT: MISMATCH -- units.json disagrees with run log.")
    else:
        print("RESULT: All checks passed.")

    return 1 if fail else 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 verify_units.py <run_dir>")
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
