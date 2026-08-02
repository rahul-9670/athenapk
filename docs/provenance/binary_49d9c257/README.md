# Provenance record — binary `49d9c257`

**Frozen automatically at build time by `docs/provenance/freeze.sh` on 2026-08-02 15:18:21.**

| item | value |
|---|---|
| binary | `athenapk/build_gpu/bin/athenaPK` |
| preserved copy | `athenapk/build_gpu/bin/athenaPK_PRESERVED_49d9c257` (hardlink, 0 bytes) |
| md5 | `49d9c25752843d25d2d9a9e44d110faa` |
| AthenaPK HEAD | `29a7174` |
| Parthenon HEAD | `fe22627999c7f67a747b723c0db11fe6dfd138b7` |

Because this was captured at build time, `tracked.patch` + `parthenon_tracked.patch` +
`untracked_src_diagnostics.tar.gz` applied on the two HEADs above **do** reproduce the source
state that produced this binary.

`object_manifest.txt` lists every translation unit linked in, with compile timestamps — use it
to settle what a given binary does and does not contain.

**Do not delete `athenaPK_PRESERVED_49d9c257`.** It costs nothing (hardlink) and is the only
artifact that survives the next rebuild.
