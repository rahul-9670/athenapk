# Provenance record — binary `f181c0a1`

**Frozen automatically at build time by `docs/provenance/freeze.sh` on 2026-08-03 09:31:57.**

| item | value |
|---|---|
| binary | `athenapk/build_gpu_v2/bin/athenaPK` |
| preserved copy | `athenapk/build_gpu_v2/bin/athenaPK_PRESERVED_f181c0a1` (hardlink, 0 bytes) |
| md5 | `f181c0a123c3093f0dd4e707faacbe00` |
| AthenaPK HEAD | `09f52bb` |
| Parthenon HEAD | `264a5e48d4930381ba05df2adc1c9b050f61c323` |

Because this was captured at build time, `tracked.patch` + `parthenon_tracked.patch` +
`untracked_src_diagnostics.tar.gz` applied on the two HEADs above **do** reproduce the source
state that produced this binary.

`object_manifest.txt` lists every translation unit linked in, with compile timestamps — use it
to settle what a given binary does and does not contain.

**Do not delete `athenaPK_PRESERVED_f181c0a1`.** It costs nothing (hardlink) and is the only
artifact that survives the next rebuild.
