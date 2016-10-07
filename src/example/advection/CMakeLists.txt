
add_library(qk_example_advection SHARED
    advect_gaussian.cpp
)
    
target_link_libraries(qk_example_advection qk)
target_link_libraries(qk_example_advection qk_grid)
target_link_libraries(qk_example_advection qk_variable)
target_link_libraries(qk_example_advection qk_solver)
