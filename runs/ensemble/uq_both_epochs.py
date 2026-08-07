#!/usr/bin/env python3
"""Run the ensemble UQ analysis at BOTH candidate measurement epochs off ONE IO sweep.

Why this exists: `uq.py` selects a member's matched snapshot by scanning every snapshot and
computing rho_max. Those snapshots are ~1.7 GB each and the density component is a strided read,
so a full pass over the 24-member campaign costs ~197 reads at ~17 s apiece (~55 min, IO-bound).
Running uq.py twice -- once per epoch -- pays that twice for *identical* rho_max values.

This driver computes rho_max once per snapshot, caches it to JSON, and then runs the ORIGINAL,
validated `uq.analyze` twice against the cache. Nothing about the selection logic, the tolerance
check, the degeneracy exclusion, or the statistics is reimplemented here -- `fr.rho_max_code` is
the only thing replaced, and only to serve a memoised value instead of re-reading the file.
`fr.measure_snapshot` is left completely untouched, so every reported number still comes from the
validated measurement path.

The cache is keyed by (path, size, mtime_ns), so a rewritten or extended snapshot invalidates its
own entry rather than silently serving a stale density.

Usage:  uq_both_epochs.py <design_dir> [--cache PATH] [--targets 1.829e5,1.829e6]
"""
import os, sys, json, glob, argparse
import numpy as np

sys.path.insert(0, "/beegfs/u/bbg6470/athenapk/runs")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flux_retention as fr
import uq

_real_rho_max_code = fr.rho_max_code


def build_cache(design_dir, cache_path):
    """One IO sweep: rho_max + time for every snapshot of every member."""
    cache = {}
    if os.path.exists(cache_path):
        try:
            cache = json.load(open(cache_path))
            print(f"[cache] loaded {len(cache)} entries from {cache_path}", flush=True)
        except Exception as e:
            print(f"[cache] ignoring unreadable cache ({e}); rebuilding", flush=True)
            cache = {}

    snaps = sorted(glob.glob(os.path.join(design_dir, "point*", "parthenon.out1.*.phdf")))
    print(f"[cache] {len(snaps)} snapshots to account for", flush=True)
    done = 0
    for s in snaps:
        st = os.stat(s)
        key = f"{s}|{st.st_size}|{st.st_mtime_ns}"
        if key in cache:
            continue
        try:
            rm, t = _real_rho_max_code(s)
            cache[key] = [rm, t]
        except Exception as e:
            # Record the failure so the sweep does not retry it every run, but store it as a
            # sentinel rather than a number -- a truncated file must never look like a density.
            cache[key] = None
            print(f"[cache] UNREADABLE {os.path.basename(s)}: {e}", flush=True)
        done += 1
        if done % 10 == 0:
            json.dump(cache, open(cache_path, "w"))
            print(f"[cache] {done} new / {len(snaps)} total ... {os.path.relpath(s, design_dir)}",
                  flush=True)
    json.dump(cache, open(cache_path, "w"))
    print(f"[cache] complete: {len(cache)} entries", flush=True)
    return cache


def install_cached_rho_max(cache):
    def cached(path):
        st = os.stat(path)
        key = f"{path}|{st.st_size}|{st.st_mtime_ns}"
        if key in cache:
            v = cache[key]
            if v is None:
                raise OSError(f"snapshot previously found unreadable: {path}")
            return float(v[0]), float(v[1])
        return _real_rho_max_code(path)      # cache miss -> fall back to the real read
    fr.rho_max_code = cached
    uq.fr.rho_max_code = cached


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("design_dir")
    ap.add_argument("--cache", default=None)
    ap.add_argument("--targets", default="1.829e5,1.829e6",
                    help="comma-separated target rho_max values [code]")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-log-dist", type=float, default=0.3)
    args = ap.parse_args()

    cache_path = args.cache or os.path.join(args.design_dir, "rho_max_cache.json")
    cache = build_cache(args.design_dir, cache_path)
    install_cached_rho_max(cache)

    for tgt in [float(x) for x in args.targets.split(",")]:
        cgs = tgt * 5.467e-19
        print("\n" + "=" * 78)
        print(f"EPOCH  target rho_max = {tgt:.4g} code = {cgs:.3e} g/cm^3 "
              f"({cgs / 1e-13:.1f}x rho_crit)")
        print("=" * 78, flush=True)
        uq.analyze(args.design_dir, tgt, args.stride, args.max_log_dist)
