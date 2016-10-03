
add_library(qk_spatial_solver SHARED
    advection.cpp
)
    
target_link_libraries(qk_spatial_solver qk)
target_link_libraries(qk_spatial_solver qk_grid)
target_link_libraries(qk_spatial_solver qk_data)
target_link_libraries(qk_spatial_solver qk_variable)
target_link_libraries(qk_spatial_solver qk_basis)