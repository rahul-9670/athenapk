# WP-0 acceptance — reproduction build: **CLOSED, PASS**

**2026-07-31, session A** (the archive itself was session B's; this closes its open item).

## Result

| | md5 | size (B) |
|---|---|---|
| original (archived) | `5ebddce0fcda54dfd133d934e4982468` | 341,615,640 |
| rebuild from this archive | `122e32ab1e671b17963d48bdd1f4b86a` | 341,619,704 |

**The md5 does NOT match. The behaviour does.**

Rebuilt at `/beegfs/u/bbg6470/wp0_repro/athenapk` (tree reset to `29a7174` + Parthenon
`fe22627`, both archived patches applied cleanly, `src/diagnostics/` restored from the tarball
— verified to contain exactly `mag_diag.{cpp,hpp}` and none of the later WP-5 work), built into
`build_gpu_repro` by `runs/submit_build_gpu_dir.sh` (job 2434896, exit 0).

Functional test: job **2434907** on g001, `runs/eos_smoke/fhc.in`, 12 cycles, 1 rank, 1 GPU,
identical deck and environment for both binaries.

> **`parthenon.out0.hst` BYTE-IDENTICAL.**

## Why md5 equality was never the right acceptance test

The binaries differ by **4,064 bytes (0.0012%)**, and the rebuild is the *larger* one. That is
consistent with longer embedded build-path strings — the reproduction tree lives at
`/beegfs/u/bbg6470/wp0_repro/athenapk/...` versus the original's
`/beegfs/u/bbg6470/athenapk/...` (+21 characters per embedded path, and `__FILE__` /
debug-path strings appear many times in a 341 MB binary). Compilation is also not
bit-reproducible in general here: no `-frandom-seed` pinning, no `-ffile-prefix-map`, and
link order can vary under `-j48`.

**Acceptance therefore = functional identity on a gate run, which PASSES.** The archive
reconstructs the production binary's behaviour exactly.

## What this licenses, and what it does not

- Licensed: every result produced by `5ebddce0` is reconstructable from this archive.
- NOT licensed: the *same md5*. Any future check must compare gate output, not checksums.
- If bit-reproducible builds are ever wanted, add `-ffile-prefix-map=$SRC=.` and a fixed
  `-frandom-seed`; not done, because functional identity is what the campaign needs.

## Same job also gated the successor binary

`build_gpu_wp20/bin/athenaPK` = `bffdf8cd4b928578b47c1a39b497a51c` (current tree: WP-20
`turb_ksample` + session B's WP-5 `grav_diag`, both default-OFF). Against the unmodified
eos_smoke deck it reproduced the original **byte-identical** on GPU — the OFF-state gate,
confirmed on real hardware rather than only on CPU — while `turb_ksample=k2` demonstrably
changes the result and prints
`k sampling : k2 (dN/dk ~ k^2 ...)  => E(k) ~ k^-1.667 [Kolmogorov = k^-1.667]`.

**A new provenance archive must be created for `bffdf8cd` before it produces science.**
