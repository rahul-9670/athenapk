# Provenance record — binary `84a6d248`

**Frozen automatically at build time by `docs/provenance/freeze.sh` on 2026-08-06 09:43:23.**

| item | value |
|---|---|
| binary | `athenapk/build_gpu/bin/athenaPK` |
| preserved copy | `athenapk/build_gpu/bin/athenaPK_PRESERVED_84a6d248` (hardlink, 0 bytes) |
| md5 | `84a6d2486bdc0545efb76c9525374355` |
| AthenaPK HEAD | `4f9adff` |
| Parthenon HEAD | `264a5e48d4930381ba05df2adc1c9b050f61c323` |

Because this was captured at build time, `tracked.patch` + `parthenon_tracked.patch` +
`untracked_src_diagnostics.tar.gz` applied on the two HEADs above **do** reproduce the source
state that produced this binary.

`object_manifest.txt` lists every translation unit linked in, with compile timestamps — use it
to settle what a given binary does and does not contain.

**Do not delete `athenaPK_PRESERVED_84a6d248`.** It costs nothing (hardlink) and is the only
artifact that survives the next rebuild.
