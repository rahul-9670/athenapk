#ifndef PGEN_PGEN_HPP_
#define PGEN_PGEN_HPP_
//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2020, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
//
// Flagship build: the only problem generator compiled in is `collapse_be`, the
// magnetized Bonnor-Ebert sphere / first-hydrostatic-core collapse used by every
// production deck (runs/root_ladder, runs/ensemble).
//
// The 23 upstream and validation problem generators that used to be declared here
// were moved to /beegfs/u/bbg6470/validation/code/pgen/ on 2026-08-10, together
// with their input decks and the regression suite.  Nothing was deleted: git tag
// `validation-complete-2026-08-10` holds the full tree.  Re-enabling one means
// restoring its .cpp, re-declaring its namespace here, adding it to
// pgen/CMakeLists.txt and adding its branch to main.cpp -- see
// /beegfs/u/bbg6470/validation/README.md.

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

namespace collapse_be {
void ProblemGenerator(parthenon::MeshBlock *pmb, parthenon::ParameterInput *pin);
// WP-13: registers the collapse_be_* Params. MUST be wired in main.cpp -- it runs on both a
// fresh start and a restart, whereas ProblemGenerator does NOT run on restart, which silently
// dropped the ApplyBarotropicCooling task (and its outside-sphere momentum BC) on every resume.
void ProblemInitPackageData(parthenon::ParameterInput *pin,
                            parthenon::StateDescriptor *hydro_pkg);
parthenon::TaskStatus ApplyBarotropicCooling(parthenon::MeshData<parthenon::Real> *md,
                                              const parthenon::SimTime &tm,
                                              const parthenon::Real dt);
} // namespace collapse_be

#endif // PGEN_PGEN_HPP_
