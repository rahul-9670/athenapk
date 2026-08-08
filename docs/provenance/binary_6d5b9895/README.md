# Provenance record — binary `6d5b9895` (`build_gpu_v5/bin/athenaPK`)

## VERDICT: THIS BINARY CORRESPONDS TO NO COMMIT. DO NOT USE IT FOR PRODUCTION OR FOR ANY
## RESULT THAT HAS TO BE REPRODUCIBLE.

Reconstructed 2026-08-08 from build artefacts, because — unlike every other entry in this
directory — it was **never frozen by `freeze.sh` at build time**. There is therefore no
`tracked.patch`, no `git_head.txt`, and no recoverable source state. What follows is what the
artefacts on disk still prove, and nothing more.

| item | value |
|---|---|
| binary | `athenapk/build_gpu_v5/bin/athenaPK` |
| md5 | `6d5b9895c86f18eb9d8f8b48d10147bc` |
| size | 349984456 bytes |
| binary mtime | 2026-08-05 16:17:23 +0200 |
| compile window | 2026-08-05 16:12:10 → 16:17:12 (171 objects total, 63 under `src/`) |
| preserved copy | **NONE** — this is the scratch slot of `build_gpu_v5`, and a `make` there overwrites it |
| git HEAD at build time | not recorded; bracketed below |

## Why it matches no commit — the evidence

The build window 16:12:10–16:17:12 falls strictly between two commits:

```
9f406ce  2026-08-05 13:23:16  audit rounds 2 and 3: cover every module round 1 did not read (N1-N12)
                    <-- build_gpu_v5 built here, 16:12-16:17
4978a4e  2026-08-05 17:17:19  N13/N14: warn when an offline table's assumed units disagree with the run's
```

So the working tree at build time was `9f406ce` **plus uncommitted edits**. That alone would be
recoverable in principle. What makes it unrecoverable is *which* file was compiled last:

```
2026-08-05 16:16:44   pgen/turbulence.cpp.o
2026-08-05 16:17:12   hydro/hydro.cpp.o     <-- last object, 28 s after the previous one
```

`src/hydro/hydro.cpp` is **one of exactly two files that the next commit `4978a4e` changed**
(the other is `src/radiation/radiation.cpp`), and it was also changed by `4f9adff` later the same
evening. A file being recompiled last, alone, well after the rest of the tree finished, is the
signature of it being edited *while the build was running*. The `hydro.cpp` inside this binary is
therefore a mid-edit snapshot: it is not `9f406ce`'s version, and it is not `4978a4e`'s version,
and no commit in the history contains it.

This is the concrete instance of the hazard recorded in the memory note
`concurrent-session-edits-2026-08-05`: a second session was editing `src/` during this build.

## What this binary is NOT

* It is **not** a v4→v5 upgrade of the production lineage. The production chain is
  `f181c0a1` (v2) → `6b1fe753` (v3) → `869c1d34` (v4) → `84a6d248` → `967fced6`, all frozen.
* It is **not** the audit-batch binary. The audit fixes were committed at `4f9adff` and `0d3a559`,
  both *after* this build.
* Its `radiation.cpp` predates the N13/N14 unit-disagreement warning, and its `hydro.cpp` is in an
  indeterminate state with respect to it. Any statement about what warnings this binary emits is
  unverifiable.

## What still references it

Four scripts, all *test harnesses*, none production:

```
runs/deep_amr/submit_d1_leg.sh
runs/deep_amr/submit_d1_meminfo.sh
runs/audit_fix_regress/submit_gpu_gate.sh
runs/gate_v5/submit_gate_v5.sh
```

Any result those produced against this binary is a result against an unidentified source tree.
That does not automatically invalidate them — the D1 memory-scaling law, for instance, is a
hardware measurement that does not depend on which `hydro.cpp` variant ran — but it does mean
**no bit-level or physics claim may be anchored on this binary**.

## Disposition

Kept, labelled, quarantined. It is ~334 MiB with no hardlink protecting it, so the next
`make -C build_gpu_v5` destroys it silently; that is acceptable, because it has no reproducible
value to lose. If space is ever needed, this is the first binary to delete.

## The rule this incident produced

`docs/provenance/freeze.sh` must run **at build time**, from inside the build job, for every
binary that will be used for anything. A provenance record reconstructed afterwards — like this
one — can bracket a build but can never recover it.
