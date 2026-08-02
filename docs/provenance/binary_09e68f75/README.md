# Provenance record — binary `09e68f75`

**Frozen automatically at build time by `docs/provenance/freeze.sh` on 2026-08-02 09:01:43.**

| item | value |
|---|---|
| binary | `athenapk/build_cpu/bin/athenaPK` |
| preserved copy | `athenapk/build_cpu/bin/athenaPK_PRESERVED_09e68f75` (hardlink, 0 bytes) |
| md5 | `09e68f75f776d5c12dfb9374bb5a5059` |
| AthenaPK HEAD | `29a7174` |
| Parthenon HEAD | `fe22627999c7f67a747b723c0db11fe6dfd138b7` |

Because this was captured at build time, `tracked.patch` + `parthenon_tracked.patch` +
`untracked_src_diagnostics.tar.gz` applied on the two HEADs above **do** reproduce the source
state that produced this binary.

`object_manifest.txt` lists every translation unit linked in, with compile timestamps — use it
to settle what a given binary does and does not contain.

**Do not delete `athenaPK_PRESERVED_09e68f75`.** It costs nothing (hardlink) and is the only
artifact that survives the next rebuild.
