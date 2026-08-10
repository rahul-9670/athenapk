# AthenaPK - a performance portable block structured AMR MHD code
# Copyright (c) 2020-2023, Athena Parthenon Collaboration. All rights reserved.
# Licensed under the 3-Clause License (the "LICENSE");

cmake_minimum_required(VERSION 3.13)

# Imports machine-specific configuration
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/parthenon/CMakeLists.txt)
  include(external/parthenon/cmake/MachineCfg.cmake)
else()
  message(WARNING
    "Could not find MachineCfg.cmake to process machine config file."
    "If you're using Parthenon as submodule and see this message, please report this bug on GitLab."
    "Using a prebuilt Parthenon library with AthenaPK is currently not tested."
    "For further directions open an issue on GitLab.")
endif()

project(athenaPK VERSION 25.05 LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

# make paths in target_source() absolute
# used for sources in src/pgen folder
# new behavior introduces in CMake 3.13 that's why this is the min. version above
cmake_policy(SET CMP0076 NEW)

# Don't allow in-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR
      "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
      "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# If the user doesn't specify a build type, prefer Release
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

include(cmake/Format.cmake)

# Kokkos itself will actually be discovered by the Parthenon package.
# Just need to set a hint to use the AthenaPK submodule and not the Parthenon one.
set(Kokkos_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/Kokkos)

# The MPI package found (and it's variables) by the Parthenon package are not available
# in the parent scope here so we need to look for MPI here as it is used in the CMake
# logic for setting up tests.
set(NUM_MPI_PROC_TESTING "4" CACHE STRING "Number of mpi processors to use when running tests with MPI")
if (NOT PARTHENON_DISABLE_MPI)
  find_package(MPI COMPONENTS CXX)
  if (NOT MPI_FOUND)
    message(FATAL_ERROR "MPI is required but couldn't be found. "
    "If you want to build Parthenon without MPI, please rerun CMake with -DPARTHENON_DISABLE_MPI=ON")
  endif()
endif()

option(AthenaPK_ENABLE_TESTING "Enable AthenaPK test" ON)
set(PARTHENON_ENABLE_PYTHON_MODULE_CHECK ${AthenaPK_ENABLE_TESTING} CACHE BOOL "Check if local python version contains all modules required for running tests.")

set(PARTHENON_ENABLE_TESTING OFF CACHE BOOL "Disable Parthenon testing.")
set(PARTHENON_DISABLE_OPENMP ON CACHE BOOL "Disable OpenMP")
set(PARTHENON_DISABLE_EXAMPLES ON CACHE BOOL "Don't build Parthenon examples.")
set(PARTHENON_DISABLE_SPARSE ON CACHE BOOL "Disable sparse (not used in AthenaPK yet)")

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/parthenon/CMakeLists.txt)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/parthenon parthenon)
else()
  find_package(parthenon REQUIRED)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_subdirectory(src)

# The upstream regression suite moved to /beegfs/u/bbg6470/validation/tst on 2026-08-10,
# together with the problem generators and input decks it drives. It is not part of the
# flagship build (the flagship compiles only the collapse_be pgen, so most of the suite
# has nothing to run against). Point AthenaPK_REGRESSION_DIR at a restored tst/regression
# to re-enable it. See /beegfs/u/bbg6470/validation/README.md.
set(AthenaPK_REGRESSION_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tst/regression"
    CACHE PATH "Path to the regression test suite (moved out of tree 2026-08-10)")
if (AthenaPK_ENABLE_TESTING)
  if (EXISTS ${AthenaPK_REGRESSION_DIR}/CMakeLists.txt)
    include(CTest)
    add_subdirectory(${AthenaPK_REGRESSION_DIR} tst_regression)
  else()
    message(STATUS
      "AthenaPK_ENABLE_TESTING is ON but no regression suite at "
      "${AthenaPK_REGRESSION_DIR} -- skipping. The suite lives in "
      "/beegfs/u/bbg6470/validation/tst.")
  endif()
endif()

