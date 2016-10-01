
add_library(qk_solvers SHARED
    solver.cpp
    advection.cpp
)
    
target_link_libraries(qk_solvers qk)
target_link_libraries(qk_solvers qk_grid)
target_link_libraries(qk_solvers qk_data)