# Copyright (c) 2020-2021, AthenaPK Collaboration. All rights reserved.
# Licensed under the 3-Clause License (the "LICENSE");

add_executable(
    athenaPK
        main.cpp
        eos/adiabatic_glmmhd.cpp
        units.hpp
        units/physical_units.hpp
        units/ionization_environment.hpp
        bvals/boundary_conditions_apk.hpp
        hydro/diffusion/ambipolar.cpp
        hydro/diffusion/conduction.cpp
        hydro/diffusion/diffusion.cpp
        hydro/diffusion/diffusion.hpp
        hydro/diffusion/hall.cpp
        hydro/diffusion/resistivity.cpp
        hydro/diffusion/viscosity.cpp
        diagnostics/angmom_diag.cpp
        diagnostics/angmom_diag.hpp
        diagnostics/cons_diag.cpp
        diagnostics/cons_diag.hpp
        diagnostics/grav_diag.cpp
        diagnostics/grav_diag.hpp
        diagnostics/mag_diag.cpp
        diagnostics/mag_diag.hpp
        hydro/hydro_driver.cpp
        hydro/hydro.cpp
        hydro/ct/ct.cpp
        hydro/ct/ct.hpp
        hydro/glmmhd/dedner_source.cpp
        hydro/prolongation/custom_ops.hpp
        refinement/gradient.cpp
        refinement/other.cpp
        refinement/jeans.cpp
        refinement/nonideal.cpp
        tracers/tracers.cpp
        tracers/tracers.hpp
        sinks/sinks.cpp
        sinks/sinks.hpp
        utils/few_modes_ft.cpp
        utils/few_modes_ft.hpp
	self_gravity/self_gravity.cpp
	self_gravity/self_gravity_driver.cpp
	radiation/radiation.cpp
	radiation/radiation_moments.cpp
	chemistry/chemistry.cpp
	chemistry/chemistry.hpp
	dust/dust.cpp
	dust/dust.hpp
	dust/dust_pkg.hpp
	chemistry/network_h2.hpp
	chemistry/network_gow17_reduced.hpp
)

add_subdirectory(pgen)

target_link_libraries(athenaPK PRIVATE parthenon)
