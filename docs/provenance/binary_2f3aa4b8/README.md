# Provenance record — binary `2f3aa4b8`

**Frozen automatically at build time by `docs/provenance/freeze.sh` on 2026-08-06 17:15:56.**

| item | value |
|---|---|
| binary | `athenapk/build_gpu/bin/athenaPK` |
| preserved copy | `athenapk/build_gpu/bin/athenaPK_PRESERVED_2f3aa4b8` (hardlink, 0 bytes) |
| md5 | `2f3aa4b89ea7f53585e593518706759b` |
| AthenaPK HEAD | `a454fee` |
| Parthenon HEAD | `264a5e48d4930381ba05df2adc1c9b050f61c323` |

Because this was captured at build time, `tracked.patch` + `parthenon_tracked.patch` +
`untracked_src_diagnostics.tar.gz` applied on the two HEADs above **do** reproduce the source
state that produced this binary.

`object_manifest.txt` lists every translation unit linked in, with compile timestamps — use it
to settle what a given binary does and does not contain.

**Do not delete `athenaPK_PRESERVED_2f3aa4b8`.** It costs nothing (hardlink) and is the only
artifact that survives the next rebuild.
