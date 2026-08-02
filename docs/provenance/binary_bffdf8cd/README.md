# Provenance record — WP-20 candidate GPU binary `bffdf8cd`

**Created 2026-08-02.** Closes the action owed by Session A's `VALIDATION_PLAN.md` close-out
("create a provenance archive for `bffdf8cd` before it produces science").

> ## ⚠️ THIS RECORD IS INCOMPLETE — THE BINARY IS **NOT** BYTE-RECONSTRUCTABLE
>
> Unlike `binary_5ebddce0/`, this archive **cannot** be used to rebuild `bffdf8cd`. The source
> state that produced it was never committed and has since been overwritten in the working
> tree. What this record *does* give you is (a) an exact statement of what the binary contains,
> (b) an exact list of what is lost, and (c) the binary itself, preserved. Read
> "What is lost" before citing this binary for anything.

## The pair

| item | value |
|---|---|
| binary | `athenapk/build_gpu_wp20/bin/athenaPK` |
| **preserved copy** | `athenapk/build_gpu_wp20/bin/athenaPK_PRESERVED_bffdf8cd` |
| md5 | `bffdf8cd4b928578b47c1a39b497a51c` |
| size / mtime | 344,314,840 B, 2026-07-31 12:07:16 |
| AthenaPK branch | `flagship-phase2-ct` |
| AthenaPK HEAD | `29a7174` ("Phase 7: IC ensemble + UQ machinery") |
| Parthenon HEAD | `fe22627999c7f67a747b723c0db11fe6dfd138b7` |

**Note the build directory.** This binary lives in `build_gpu_wp20/`, **not** `build_gpu/`. It
never overwrote production: `build_gpu/bin/athenaPK` is still `5ebddce0`, mtime 2026-07-30
19:50, matching `docs/provenance/binary_5ebddce0/`. (An earlier session note claiming "the
`build_gpu` binary md5 changed" was wrong; corrected 2026-08-02.)

### The binary is preserved by hardlink

`athenaPK_PRESERVED_bffdf8cd` is a **hardlink** to the same inode (verified: identical inode,
`links=2`, md5 matches). A rebuild into `build_gpu_wp20` writes a *new* file and unlinks the
old name, so the preserved name keeps the inode alive. Cost is **zero bytes** — `du -ch` over
both names totals 329 M, not 658 M. **Do not delete it**; it is the only surviving artifact of
this binary's source state.

## What produced science on it

`grep -rl bffdf8cd runs/ --include=*.sh`:

- `runs/wp18_seed_ensemble/submit_seed.sh` — **WP-18's 12-seed IC ensemble**, the measurement
  that set the campaign's σ (~16 % `MEtor/MEpol` at matched state). This is the load-bearing one.
- `runs/wp0_verify/submit_verify.sh`
- `runs/root_ladder/submit_root.sh`

## What the binary contains — EXACT, from the object manifest

`object_manifest.txt` lists all **61** translation units linked into it, with compile
timestamps. This is direct evidence from the build tree, not inference. The consequential part:

| diagnostic | in `bffdf8cd`? | evidence |
|---|---|---|
| `grav_diag` | **yes** | object at 12:05:15 |
| `mag_diag` | **yes** | object at 12:05:14 |
| `cons_diag` | **NO** | no object; added to `src/CMakeLists.txt` at 14:50, after the 12:07 link |
| `angmom_diag` | **NO** | no object; same |

So `bffdf8cd` has **no angular-momentum and no conservation-budget history columns**. Any
`.hst` it produced cannot contain `am-*` or `cons-*` columns, and an analysis script that
expects them is reading the wrong run. Conversely, WP-4/WP-6 results can never have come from
this binary.

## What is lost, precisely

`DRIFT_after_build.txt` lists every source file modified after the 12:07 link. Seven files
drifted, but only these were **actually compiled into** the binary, so only these matter:

| file | modified | why it moved |
|---|---|---|
| `src/diagnostics/mag_diag.hpp` | 2026-07-31 13:41:07 | Session A added `mag_diag_rho_split` |
| `src/diagnostics/mag_diag.cpp` | 2026-07-31 13:41:18 | same |
| `src/hydro/hydro.cpp` | 2026-07-31 19:18:24 | WP-4/WP-6 diagnostic registration |
| `src/CMakeLists.txt` | 2026-07-31 14:50:50 | added `cons_diag` + `angmom_diag` |

The remaining three drifted files (`cons_diag.{cpp,hpp}`, `angmom_diag.{cpp,hpp}`) were **not**
in the binary at all, so their drift is irrelevant to reconstruction.

**Why the 12:07 content is unrecoverable:** those edits were uncommitted working-tree
modifications. `git log` shows no commit since before 2026-07-31, `git stash list` holds only an
unrelated old entry, and the reflog's newest entry is HEAD `29a7174` itself. Saving the tree now
would pin `bffdf8cd` to a source state that never produced it — which is exactly the failure this
record exists to prevent, so it is deliberately **not** done.

### What *is* exactly recoverable

- **The entire Parthenon half.** `parthenon_tracked_CURRENT.patch` is **byte-identical** to
  `docs/provenance/binary_5ebddce0/parthenon_tracked.patch` — the submodule has not been touched
  since. Applied on `fe22627`, this reproduces the Parthenon source exactly as built.
- All HEAD and submodule SHAs.
- The complete object manifest (what was linked, and when).

## Files

| file | what | trustworthy for `bffdf8cd`? |
|---|---|---|
| `binary.md5` | md5 + `ls -l` of the binary and its preserved hardlink | **yes** |
| `object_manifest.txt` | all 61 linked TUs with compile times | **yes — the key evidence** |
| `DRIFT_after_build.txt` | every source file touched after the link | **yes** |
| `git_head.txt`, `parthenon_head.txt`, `submodules.txt` | SHAs | **yes** |
| `parthenon_tracked_CURRENT.patch` | Parthenon `git diff HEAD` | **yes** — byte-identical to the 5ebddce0 archive |
| `tracked_CURRENT.patch` | AthenaPK `git diff HEAD`, **as of 2026-08-02** | **NO** — includes post-build edits; kept for delta analysis only |
| `untracked_src_diagnostics_CURRENT.tar.gz` | `src/diagnostics/`, **as of 2026-08-02** | **NO** — contains `cons_diag`/`angmom_diag`, which are *not* in this binary |

The two `_CURRENT` files are named that way on purpose: they describe **today's tree**, not the
binary. Do not apply them expecting `bffdf8cd`.

## How to cite this binary

For WP-18's σ result, cite the **binary** (md5 + preserved path), not a source state. The
ensemble is reproducible by re-running with the preserved binary; it is *not* reproducible by
rebuilding from source. If bit-level source provenance is required for publication, WP-18 must
be re-run on a properly frozen binary.

## The lesson, and the fix

The `5ebddce0` archive was made *because* a binary was about to produce science. `bffdf8cd`
produced science first and got archived second — by which time the tree had moved three times.

**Rule going forward: archive at build time, in the same script that builds.** A build that
cannot be archived should not be run. Add the freeze to `runs/submit_build_gpu.sh` so it is not
a step someone has to remember.

Heed the `5ebddce0` finding as well: `git status` in AthenaPK reports `M external/parthenon`, a
*pointer* change that does **not** carry the modified files inside the submodule. The submodule
diff must be archived separately or it is silently lost.
