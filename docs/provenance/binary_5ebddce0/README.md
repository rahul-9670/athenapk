# Provenance record — production GPU binary `5ebddce0`

**Created 2026-07-31 (WP-0 of `VALIDATION_PLAN.md`).**

This directory freezes the exact source state that produced the GPU binary used by
the entire 2026-07-31 validation campaign. Every subsequent WP must cite this pair.

## The pair

| item | value |
|---|---|
| binary | `athenapk/build_gpu/bin/athenaPK` |
| md5 | `5ebddce0fcda54dfd133d934e4982468` |
| size / mtime | 341,615,640 B, 2026-07-30 19:50:52 |
| AthenaPK branch | `flagship-phase2-ct` |
| AthenaPK HEAD | `29a7174` ("Phase 7: IC ensemble + UQ machinery") |
| Parthenon HEAD | `fe22627999c7f67a747b723c0db11fe6dfd138b7` (v25.12-194) |
| Kokkos | `6739bc62` (4.7.02) — clean |

**The binary is NOT reconstructable from git alone.** It requires the two patches below
applied on top of those two SHAs, plus the untracked tarball.

## Files

| file | what |
|---|---|
| `git_head.txt` / `git_status.txt` | AthenaPK HEAD + 20 modified tracked files |
| `tracked.patch` | `git diff HEAD` in AthenaPK (247,798 B) |
| `untracked_src_diagnostics.tar.gz` | untracked `src/diagnostics/` (`cap_diag`, `mag_diag`) |
| `parthenon_head.txt` / `parthenon_status.txt` | Parthenon HEAD + 9 modified files |
| `parthenon_tracked.patch` | `git diff HEAD` inside the submodule (15,396 B) |
| `submodules.txt` | recursive submodule SHAs |
| `binary.md5` | md5 + `ls -l` of the binary |

### Finding: the submodule diff is invisible to the parent repo

`git status` in AthenaPK reports `M external/parthenon` — a *pointer* modification. That
line does **not** carry the 9 modified files inside Parthenon
(`bvals/comms/{boundary_communication,build_boundary_buffers,coalesced_buffers}.cpp`,
`coalesced_buffers.hpp`, `driver/driver.{cpp,hpp}`, `interface/update.cpp`, `mesh/mesh.cpp`,
`solvers/mg_solver.hpp` — the coalesced-comms CT fix and the MG solver work).
Archiving only the parent `git diff` would have silently lost them. Any future freeze of
this workspace must archive the submodule separately.

## Restoration recipe

```bash
source ~/athenapk_env.sh
cd /path/to/fresh/athenapk
git checkout 29a7174
git submodule update --init --recursive
git apply .../tracked.patch
tar xzf .../untracked_src_diagnostics.tar.gz
cd external/parthenon && git apply .../parthenon_tracked.patch
```

## Open acceptance item — rebuild reproduction

WP-0's acceptance also requires rebuilding from this state and confirming the md5
reproduces (or documenting why not — build nondeterminism is itself a finding).

**Deliberately NOT done on 2026-07-31.** The self-resuming submit scripts re-exec
`build_gpu/bin/athenaPK` at every slot boundary, so rebuilding into `build_gpu` while the
njeans ladder (jobs 2431681/82/83 + successors) and root ladder are live would swap the
binary mid-chain and destroy the campaign's internal consistency — exactly what
`VALIDATION_PLAN.md` §1 certifies is currently true.

The reproduction build must therefore go into a **separate** build directory
(`build_gpu_repro`), and only when GPU nodes are free. Tracked as the remaining half of WP-0.
