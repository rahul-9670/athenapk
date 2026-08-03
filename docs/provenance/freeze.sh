#!/bin/bash
# Freeze the source state that produced a binary, AT BUILD TIME.
#
# WHY THIS EXISTS. `binary_5ebddce0` was archived because a binary was about to produce
# science. `bffdf8cd` produced science FIRST and was archived second -- by which time the
# working tree had moved three times, and its source state is now permanently unrecoverable
# (the edits were never committed). See docs/provenance/binary_bffdf8cd/README.md.
#
# The fix is to make the freeze a step of the BUILD, not a thing someone remembers afterwards.
#
# Usage:  freeze.sh <build_dir>          e.g. freeze.sh build_gpu
# Writes: docs/provenance/binary_<md5-8>/
#
# Safe to re-run: if the archive for that md5 already exists it is left alone and the script
# exits 0 (same binary => same source state => nothing to add). Never overwrites an archive.
# Always exits 0 -- a provenance failure must not fail the build that produced the binary.

set -uo pipefail
REPO=/beegfs/u/bbg6470/athenapk
BUILD=${1:-build_gpu}
BIN=$REPO/$BUILD/bin/athenaPK

if [ ! -f "$BIN" ]; then
  echo "freeze.sh: no binary at $BIN -- nothing to freeze (build probably failed)"; exit 0
fi

# GIT MUST BE ON PATH, and it is NOT there by default on this cluster -- `git` lives in the
# pkgsrc bin that ~/athenapk_env.sh adds (or `module load git/2.48.1`). Every git call below
# redirects stderr to /dev/null, so without this check a git-less shell writes SIX EMPTY FILES
# and the archive looks complete while carrying no source state at all. That is exactly the
# failure the whole freeze mechanism exists to prevent, and it happened to binary f181c0a1
# (2026-08-03) when freeze.sh was run by hand from a shell that had not sourced the env.
if ! command -v git >/dev/null 2>&1; then
  echo "freeze.sh: **GIT NOT ON PATH** -- refusing to write an archive with empty git state."
  echo "freeze.sh: run 'source ~/athenapk_env.sh' (or 'module load git/2.48.1') and retry:"
  echo "freeze.sh:   bash $REPO/docs/provenance/freeze.sh $BUILD"
  exit 0
fi

MD5=$(md5sum "$BIN" | awk '{print $1}')
SHORT=${MD5:0:8}
A=$REPO/docs/provenance/binary_$SHORT

# 1. Preserve the binary itself by HARDLINK -- zero bytes, and it survives the next rebuild
#    (a rebuild writes a new file and unlinks the old name; the inode lives on under this one).
#    NOTE: this MUST happen before the "archive already exists" early-return below. An earlier
#    version returned first, so re-running against an already-archived binary left it
#    unpreserved -- and the next rebuild would have destroyed it. That nearly cost the
#    archived production binary 5ebddce0 (2026-08-02).
PRES="$REPO/$BUILD/bin/athenaPK_PRESERVED_$SHORT"
if [ ! -e "$PRES" ]; then
  ln "$BIN" "$PRES" 2>/dev/null && echo "freeze.sh: preserved $PRES"
fi

if [ -d "$A" ]; then
  echo "freeze.sh: $A already exists for md5 $MD5 -- leaving the archive untouched"; exit 0
fi
mkdir -p "$A" || { echo "freeze.sh: cannot create $A"; exit 0; }

cd "$REPO" || exit 0

# 2. Identity of the binary.
{ md5sum "$BIN"; [ -e "$PRES" ] && md5sum "$PRES"; ls -l --time-style=+%F_%T "$BIN"; } > "$A/binary.md5"

# 3. Git state -- parent repo AND the submodule. The submodule diff is invisible to the
#    parent: `git status` shows only `M external/parthenon`, a POINTER change, which does not
#    carry the modified files inside Parthenon. Archive it separately or it is silently lost.
git rev-parse HEAD                        > "$A/git_head.txt"        2>/dev/null
git status --porcelain                    > "$A/git_status.txt"      2>/dev/null
git diff HEAD                             > "$A/tracked.patch"       2>/dev/null
git -C external/parthenon rev-parse HEAD  > "$A/parthenon_head.txt"  2>/dev/null
git -C external/parthenon status --porcelain > "$A/parthenon_status.txt" 2>/dev/null
git -C external/parthenon diff HEAD       > "$A/parthenon_tracked.patch"  2>/dev/null
git submodule status --recursive          > "$A/submodules.txt"      2>/dev/null

# 4. Untracked sources that are compiled in but invisible to git diff.
tar czf "$A/untracked_src_diagnostics.tar.gz" src/diagnostics/ 2>/dev/null

# 5. What was ACTUALLY linked, with compile times. This is direct evidence from the build
#    tree rather than inference, and it is what made the bffdf8cd post-mortem possible at all
#    (it proved cons_diag/angmom_diag were absent from that binary).
find "$BUILD/src/CMakeFiles/athenaPK.dir" -name "*.o" -printf "%T+  %p\n" 2>/dev/null \
  | sed "s|$BUILD/src/CMakeFiles/athenaPK.dir/||" | sort > "$A/object_manifest.txt"

cat > "$A/README.md" <<EOF
# Provenance record — binary \`$SHORT\`

**Frozen automatically at build time by \`docs/provenance/freeze.sh\` on $(date +%F\ %T).**

| item | value |
|---|---|
| binary | \`athenapk/$BUILD/bin/athenaPK\` |
| preserved copy | \`athenapk/$BUILD/bin/athenaPK_PRESERVED_$SHORT\` (hardlink, 0 bytes) |
| md5 | \`$MD5\` |
| AthenaPK HEAD | \`$(git rev-parse --short HEAD 2>/dev/null)\` |
| Parthenon HEAD | \`$(git -C external/parthenon rev-parse HEAD 2>/dev/null)\` |

Because this was captured at build time, \`tracked.patch\` + \`parthenon_tracked.patch\` +
\`untracked_src_diagnostics.tar.gz\` applied on the two HEADs above **do** reproduce the source
state that produced this binary.

\`object_manifest.txt\` lists every translation unit linked in, with compile timestamps — use it
to settle what a given binary does and does not contain.

**Do not delete \`athenaPK_PRESERVED_$SHORT\`.** It costs nothing (hardlink) and is the only
artifact that survives the next rebuild.
EOF

echo "freeze.sh: wrote provenance archive $A"
exit 0
