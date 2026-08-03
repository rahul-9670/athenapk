#!/bin/bash
#SBATCH --job-name=build_gpu_dir
#SBATCH --account=banerjee_std
#SBATCH --partition=std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=03:00:00
#SBATCH --output=/beegfs/u/bbg6470/athenapk/runs/build_%x_%j.out
set -o pipefail
#
# Parameterized GPU build into an ARBITRARY source tree + build dir. Exists because
# runs/submit_build_gpu.sh hardcodes /beegfs/u/bbg6470/athenapk/build_gpu, and we must NOT
# rebuild in place: the self-resuming submit scripts re-exec build_gpu/bin/athenaPK at every
# slot boundary, and the linker cannot even write an executable that a running job has mapped
# (ETXTBSY). Env:
#   SRC   = source tree (default the main workspace)
#   BDIR  = build directory to configure+build (created if absent)
#
# CMake args replicate build_gpu/CMakeCache.txt exactly (checked 2026-07-31):
#   Release, nvcc_wrapper, CUDA on, HOPPER90, OpenMP OFF, sparse/examples/HDF5-compression off.
source ~/athenapk_env.sh
module load cuda/12.5.1
SRC="${SRC:-/beegfs/u/bbg6470/athenapk}"
BDIR="${BDIR:?set BDIR}"
echo "SRC=$SRC"; echo "BDIR=$BDIR"; echo "nvcc: $(which nvcc)"; date

if [ ! -f "$BDIR/CMakeCache.txt" ]; then
  echo "=== configuring fresh ==="
  cmake -S "$SRC" -B "$BDIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER="$SRC/external/parthenon/external/Kokkos/bin/nvcc_wrapper" \
    -DHDF5_ROOT=/beegfs/u/bbg6470/hdf5_parallel_install \
    -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON \
    -DKokkos_ENABLE_CUDA_LAMBDA=ON -DKokkos_ENABLE_OPENMP=OFF -DKokkos_ENABLE_SERIAL=ON \
    -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=ON -DKokkos_ENABLE_COMPLEX_ALIGN=ON \
    -DPARTHENON_DISABLE_OPENMP=ON -DPARTHENON_DISABLE_EXAMPLES=ON \
    -DPARTHENON_DISABLE_SPARSE=ON -DPARTHENON_DISABLE_HDF5_COMPRESSION=ON \
    -DPARTHENON_COPYRIGHT_CHECK_DEFAULT=OFF \
    -DPARTHENON_ENABLE_PYTHON_MODULE_CHECK=OFF \
    || { echo "CMAKE_FAILED"; exit 1; }
    # PARTHENON_ENABLE_PYTHON_MODULE_CHECK=OFF is REQUIRED here: Parthenon's configure
    # checks for numpy/unyt/matplotlib/h5py/scipy in the *system* python3.12 (which has
    # none of them -- this workspace keeps them in venvs/analysis_env) and hard-fails at
    # cmake/PythonModuleCheck.cmake:44. The existing build_gpu cache carries this flag OFF;
    # omitting it fails configure in 17 s (job 2434894). It only disables regression tests.
fi

# NOTE: no "make clean" -- clean deletes generated/provenance.cpp (a byproduct of
# ParthenonAlwaysCheckGit) and the parthenon target compiles it without depending on that
# custom target, so -j races and fails (job 2321949, 2026-07-11).
make -C "$BDIR" athenaPK -j48
RC=$?
echo "BUILD_DONE exit=$RC"; date
ls -la "$BDIR/bin/athenaPK" 2>/dev/null && md5sum "$BDIR/bin/athenaPK"
