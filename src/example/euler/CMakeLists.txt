
add_library(qk_example_euler SHARED
    sod_shock_tube.cpp
)
    
target_link_libraries(qk_example_euler qk)
target_link_libraries(qk_example_euler qk_grid)
target_link_libraries(qk_example_euler qk_variable)
target_link_libraries(qk_example_euler qk_solver)
