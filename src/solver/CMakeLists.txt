
add_library(qk_solver SHARED
    swap.cpp
    ssprk2_advection.cpp
    gaussian.cpp
    bc_periodic.cpp
)
    
target_link_libraries(qk_solver qk)
target_link_libraries(qk_solver qk_grid)
target_link_libraries(qk_solver qk_data)
target_link_libraries(qk_solver qk_variable)
target_link_libraries(qk_solver qk_basis)