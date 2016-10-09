
add_library(qk_example_plasma SHARED
    two_fluid_brio_wu.cpp
    harris_current_sheet.cpp
    gem_challenge.cpp
)
    
target_link_libraries(qk_example_plasma qk)
target_link_libraries(qk_example_plasma qk_grid)
target_link_libraries(qk_example_plasma qk_variable)
target_link_libraries(qk_example_plasma qk_solver)
