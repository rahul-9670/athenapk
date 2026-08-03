#!/bin/bash
#SBATCH --job-name=build_gpu
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=01:30:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/build_gpu_compile_%j.out
set -o pipefail
source ~/athenapk_env.sh
module load cuda/12.5.1
echo "nvcc: $(which nvcc)"
cd /beegfs/u/bbg6470/athenapk
# NOTE(2026-07-11): no "make clean" here — clean deletes generated/provenance.cpp (a
# BYPRODUCT of ParthenonAlwaysCheckGit) and the parthenon target compiles it without a
# dependency on that custom target, so -j48 races and fails (job 2321949). Incremental
# builds are correct for source edits; after CMake changes, re-run cmake (which
# regenerates provenance.cpp) instead of cleaning.
make -C build_gpu athenaPK -j48
echo "BUILD_DONE exit=$?"
ls -la /beegfs/u/bbg6470/athenapk/build_gpu/bin/athenaPK
# Freeze source->binary provenance IMMEDIATELY, while the tree still matches the binary.
# NOTE(2026-08-02): this is here because `bffdf8cd` produced science (WP-18's 12-seed
# ensemble) and was archived only afterwards -- by which time the working tree had moved
# three times and its source state was permanently unrecoverable, since the edits were never
# committed. See docs/provenance/binary_bffdf8cd/README.md. freeze.sh is idempotent (an
# existing archive for the same md5 is left alone) and always exits 0, so it can never fail
# the build that produced the binary.
/beegfs/u/bbg6470/athenapk/docs/provenance/freeze.sh build_gpu
