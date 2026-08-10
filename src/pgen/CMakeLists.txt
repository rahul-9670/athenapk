# Copyright (c) 2020-2021, Athena Parthenon Collaboration. All rights reserved.
# Licensed under the 3-Clause License (the "LICENSE");

# Flagship problem generator only.
#
# The 23 upstream/validation problem generators (advection, blast, cloud, cluster,
# cpaw, cshock, diffusion, field_loop, jeans, kh, linear_wave, linear_wave_mhd,
# lw_implode, orszag_tang, poisson_test, polytrope, rad_pulse, rad_shadow,
# rad_shock, rand_blast, shattering, sod, turbulence) were moved out of the build
# on 2026-08-10 -- see /beegfs/u/bbg6470/validation/code/pgen/ for the sources and
# /beegfs/u/bbg6470/validation/README.md for how to put one back.  They are also
# preserved in git at tag `validation-complete-2026-08-10`.

target_sources(athenaPK PRIVATE
    collapse_be.cpp
    )
