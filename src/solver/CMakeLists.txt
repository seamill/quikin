
add_library(qk_solver SHARED
    swap.cpp
    print.cpp
    ssprk2_advection.cpp
    gaussian.cpp
    bc_periodic.cpp
    bc_copy_out.cpp
    bc_conducting_wall.cpp
    bc_no_slip.cpp
    bc_free_slip.cpp
    ssprk3.cpp
    rk1.cpp
    euler.cpp
    maxwell.cpp
    lorentz_force.cpp
    current_source.cpp
    fill.cpp
    shock_tube.cpp
    radial_shock_tube.cpp
    magnetic_pulse.cpp
    brio_wu_shock_tube.cpp
    gem_challenge.cpp
    electromagnetic_wave_1D.cpp
)
    
target_link_libraries(qk_solver qk)
target_link_libraries(qk_solver qk_grid)
target_link_libraries(qk_solver qk_data)
target_link_libraries(qk_solver qk_variable)
target_link_libraries(qk_solver qk_basis)