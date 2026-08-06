#!/usr/bin/env python
"""Stop ensemble members that have reached the campaign's matched epoch but are still running.

WHY THIS EXISTS. `submit_point.sh` evaluates the rho >= STOP_CGS stop only BETWEEN chain slots,
at job start. A member that crosses the threshold mid-slot keeps integrating to the end of that
slot. That is the single most expensive compute in the run: past first-core formation dt collapses
and the cost per cycle climbs steeply, so the surplus is not a small tail. Measured on 2026-08-06:
with 4 h slots, point002 reached rho = 2.57e-11, i.e. **256x past** its own 1e-13 finish line, and
point000 20x. Shortening slots to 2 h bounds it (point003 overshot only 2.3x) but does not close it.

WHAT IT DOES. Every POLL_S seconds: read each member's newest READABLE snapshot; if rho_max is at
or past STOP_CGS and that member still has a running job, write STOP_CHAIN (so the queued successor
exits immediately at its own first check) and cancel the running slot. Data is never deleted --
snapshots and restarts stay on disk, and everything past the threshold is surplus to this campaign.

SAFETY, deliberately conservative -- this cancels real jobs:
  * Acts ONLY on a positive reading. An unreadable snapshot returns -1 and is ignored, never
    treated as "0, so keep going" (the bug that used to live in submit_point.sh) nor as a reason
    to cancel.
  * Writes STOP_CHAIN BEFORE cancelling, so a successor released by `afterany` sees the marker and
    exits cleanly instead of restarting the member.
  * MAX_ACTIONS_PER_PASS caps how much damage a mis-read could do in one sweep.
  * Idempotent: a member already carrying STOP_CHAIN is skipped.
  * Job -> member mapping is exact, not by name-guessing: original members are named
    `ens_point0NN`; chained successors are named plain `ens` and are identified by their WorkDir,
    which Parthenon's submit path sets to the member directory.

Prints one line per ACTION (or anomaly). Silence means everything is where it should be.
"""
import os
import subprocess
import sys
import time

import glob
import h5py
import numpy as np

DESIGN = "/beegfs/u/bbg6470/athenapk/runs/ensemble/design01"
RHO0 = 5.467e-19
# Overridable so the ACTION path can be exercised in DRYRUN against a live campaign by setting a
# threshold the running members have already passed. Production value is the campaign's 1e-13.
STOP_CGS = float(os.environ.get("STOP_CGS", "1.0e-13"))
POLL_S = 300
MAX_ACTIONS_PER_PASS = 6
USER = "bbg6470"
# DRYRUN=1 -> report the decisions it WOULD take and exit after one pass, writing nothing and
# cancelling nothing. Always run this before arming it against a live campaign.
DRYRUN = os.environ.get("DRYRUN", "") == "1"


def sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True,
                              timeout=120).stdout
    except Exception:
        return ""


def rho_max_cgs(d):
    """Newest READABLE snapshot. Returns (-1, 'none') if nothing can be read."""
    fs = sorted(glob.glob(os.path.join(d, "parthenon.out1.*.phdf")),
                key=os.path.getmtime, reverse=True)
    for f in fs:
        try:
            with h5py.File(f, "r") as h:
                return float(np.array(h["prim"][:, 0, ...]).max()) * RHO0, os.path.basename(f)
        except Exception:
            continue
    return -1.0, "none"


def running_jobs(states="RUNNING"):
    """member name -> [jobid], for the given SLURM states. Exact mapping; see docstring."""
    out = {}
    raw = sh(f"squeue -u {USER} -h -t {states} -o '%i|%j'")
    for line in raw.strip().splitlines():
        if "|" not in line:
            continue
        jid, name = line.split("|", 1)
        jid, name = jid.strip(), name.strip()
        member = None
        if name.startswith("ens_point"):
            member = name[len("ens_"):]
        elif name == "ens":
            wd = sh(f"scontrol show job {jid}")
            for tok in wd.split():
                if tok.startswith("WorkDir="):
                    member = os.path.basename(tok.split("=", 1)[1].rstrip("/"))
                    break
        if member and member.startswith("point"):
            out.setdefault(member, []).append(jid)
    return out


def main():
    # COST CONTROL. Reading rho_max means pulling the density field out of a ~1.8 GB snapshot, and
    # a first version did that for all 24 members every pass -- ~2 minutes of BeeGFS I/O per sweep,
    # for a campaign that runs for days. Two fixes: (a) only members with a RUNNING job are read,
    # since those are the only ones this can act on at all (members between slots are closed by the
    # successor's own check); (b) a (path, mtime) cache, so an unchanged snapshot is never re-read.
    # Per-member cache keyed on the newest snapshot's (path, mtime). An earlier version called
    # cache.clear() on every miss, so with 24 members only ONE entry ever survived and the cache
    # did nothing. Keyed per member, a member whose snapshot has not changed is never re-read.
    cache = {}

    def rho_cached(d, member):
        fs = sorted(glob.glob(os.path.join(d, "parthenon.out1.*.phdf")),
                    key=os.path.getmtime, reverse=True)
        if not fs:
            return -1.0, "none"
        key = (fs[0], os.path.getmtime(fs[0]))
        hit = cache.get(member)
        if hit and hit[0] == key:
            return hit[1]
        val = rho_max_cgs(d)
        cache[member] = (key, val)
        return val

    while True:
        try:
            jobs = running_jobs()
            actions = 0
            # Sweep EVERY member, not just running ones. A member that crossed the stop between
            # slots still needs STOP_CHAIN written: without it nothing marks it finished, its
            # queued successor sits in line only to start and immediately exit, and the campaign
            # tally reads 0-finished while members are in fact done. (Observed 2026-08-06: five
            # members past the stop, 0 of 24 carrying STOP_CHAIN.) The mtime cache keeps this
            # cheap -- an unchanged snapshot is never re-read.
            # A finished member can still have queued jobs: its own next slot and/or an `afterany`
            # successor submitted before it crossed the line. Each of those will start, allocate
            # 4 GPUs, read STOP_CHAIN, and exit -- a scheduling slot spent on nothing while real
            # members wait. Reap them. (Observed 2026-08-06: 7 such jobs against 7 finished
            # members.) Safe by construction: only members already carrying STOP_CHAIN are touched,
            # nothing is running, and no chain is broken because the member is done.
            queued = running_jobs(states="PENDING")

            for d in sorted(glob.glob(os.path.join(DESIGN, "point*/"))):
                member = os.path.basename(d.rstrip("/"))
                stop_file = os.path.join(d, "STOP_CHAIN")
                if os.path.exists(stop_file):
                    stale = queued.get(member, [])
                    if stale:
                        if DRYRUN:
                            print(f"WOULD REAP {member}: finished, but {len(stale)} job(s) still "
                                  f"queued: {','.join(stale)}", flush=True)
                        else:
                            sh(f"scancel {' '.join(stale)}")
                            print(f"WATCHDOG REAPED {member}: finished; cancelled "
                                  f"{len(stale)} queued job(s) {','.join(stale)}", flush=True)
                    continue                       # already closed
                rho, fname = rho_cached(d, member)
                if rho < 0:
                    continue                       # unreadable -> never act
                if rho < STOP_CGS:
                    continue                       # still climbing, legitimately
                over = rho / STOP_CGS
                jids = jobs.get(member, [])
                if not jids:
                    # Past the stop with nothing running: mark it done so the tally is honest and
                    # the queued successor exits at its very first check.
                    if DRYRUN:
                        print(f"WOULD MARK DONE {member}: rho={rho:.3e} cgs = {over:.1f}x past "
                              f"stop, from {fname} (no running job)", flush=True)
                    else:
                        with open(stop_file, "w") as fh:
                            fh.write(f"watchdog epoch stop rho={rho:.6e} cgs ({over:.1f}x past "
                                     f"{STOP_CGS:.0e}) from {fname} — no running job "
                                     f"{time.strftime('%F %T')}\n")
                        print(f"WATCHDOG DONE {member}: rho={rho:.3e} cgs = {over:.1f}x past "
                              f"stop; marked finished (nothing was running)", flush=True)
                    continue
                if actions >= MAX_ACTIONS_PER_PASS:
                    print(f"WATCHDOG cap reached ({MAX_ACTIONS_PER_PASS}); deferring {member}",
                          flush=True)
                    break
                if DRYRUN:
                    print(f"WOULD STOP {member}: rho={rho:.3e} cgs = {over:.1f}x past stop, "
                          f"from {fname}, would cancel {','.join(jids)}", flush=True)
                    actions += 1
                    continue
                with open(stop_file, "w") as fh:
                    fh.write(f"watchdog epoch stop rho={rho:.6e} cgs ({over:.1f}x past "
                             f"{STOP_CGS:.0e}) from {fname} {time.strftime('%F %T')}\n")
                sh(f"scancel {' '.join(jids)}")
                actions += 1
                print(f"WATCHDOG STOPPED {member}: rho={rho:.3e} cgs = {over:.1f}x past the "
                      f"1e-13 stop; cancelled job(s) {','.join(jids)} (data intact)", flush=True)
        except Exception as e:                     # never let a bad pass kill the watchdog
            print(f"WATCHDOG error (continuing): {type(e).__name__}: {e}", flush=True)
        if DRYRUN:
            print("DRYRUN complete — nothing written, nothing cancelled", flush=True)
            return 0
        time.sleep(POLL_S)


if __name__ == "__main__":
    sys.exit(main())
