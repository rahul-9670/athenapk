# Minimal reproducer: multi-rank GPU divergence of GMG-preconditioned BiCGSTAB Poisson solve

Relates to parthenon-hpc-lab/parthenon#1400 and #1364.

## Symptom
On >=2 GPU ranks, the GMG-preconditioned BiCGSTAB Poisson solve diverges (residual
stalls high, grav.phi goes NaN). On 1 GPU rank, or on CPU, or with
CUDA_LAUNCH_BLOCKING=1, it converges. Same binary, same input, same rank count;
only the execution mode differs. This is the async-execution race described in #1364,
reached via the multi-stage MG/BiCGSTAB communication path.

## Config
- AthenaPK self-gravity port, branch rahul/self-gravity-port
- Parthenon 287313a7 (v25.12), Kokkos 4.7.02
- 64^3 root grid, 16^3 meshblocks (64 base blocks), AMR Jeans criterion numlevel=2
- Dirichlet ("zero") gravity BCs, BiCGSTAB + Multigrid preconditioner
- No MHD turbulence (turb_mach=0), minimal deterministic setup

## Build (H100 / Hopper, CUDA 12.5.1)
Configure with the Parthenon-nested Kokkos nvcc_wrapper (not a top-level one):

    cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON \
          -DCMAKE_CXX_COMPILER=<repo>/external/parthenon/external/Kokkos/bin/nvcc_wrapper \
          -DHDF5_ROOT=<parallel-hdf5-install> \
          -DPARTHENON_DISABLE_EXAMPLES=ON -DPARTHENON_ENABLE_PYTHON_MODULE_CHECK=OFF ..
    make -j athenaPK

## Reproduce
DIVERGES (2 GPU ranks, no blocking):

    mpirun -n 2 ./athenaPK -i repro.in

CONVERGES (2 GPU ranks, with blocking):

    CUDA_LAUNCH_BLOCKING=1 mpirun -n 2 ./athenaPK -i repro.in

(print_per_step = true prints the per-iteration residual.)

## Observed (2x H100)
- no blocking:  3.06 -> 7.6 -> 7.4 -> ... (stalls, diverges)
- blocking=1:   3.06 -> 1.28 -> 0.099 -> 5.9e-3 -> ... -> 1.3e-11 (converges)

## Notes
- parthenon/mesh/do_coalesced_comms=true (the #1364 mitigation) cannot be applied here:
  it aborts with "coalesced comms and multiple communication stages can't be used
  concurrently" (mesh.hpp), because the BiCGSTAB+MG solver uses multiple comm stages.
- Operator-level and Jacobi-smoother Kokkos::fence() were each insufficient; only
  global serialization (CUDA_LAUNCH_BLOCKING=1) works.
- Unmodified Artemis self-gravity runs clean on the same stack, so this is specific
  to the AthenaPK comm path, not the environment.
