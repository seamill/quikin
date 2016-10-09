
add_library(qk_example_maxwell SHARED
    gaussian_pulse.cpp
    wave_1D.cpp
)
    
target_link_libraries(qk_example_maxwell qk)
target_link_libraries(qk_example_maxwell qk_grid)
target_link_libraries(qk_example_maxwell qk_variable)
target_link_libraries(qk_example_maxwell qk_solver)
