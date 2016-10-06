
add_subdirectory(lib)
add_subdirectory(solver)
add_subdirectory(spatial_solver)
add_subdirectory(temporal_solver)
add_subdirectory(data)
add_subdirectory(variable)
add_subdirectory(grid)
add_subdirectory(basis)

add_executable(quikin main.cpp)

link_directories(lib solvers spatial_solver temporal_solver data variable grid basis)

target_link_libraries(quikin qk)
target_link_libraries(quikin qk_solver)
target_link_libraries(quikin qk_spatial_solver)
target_link_libraries(quikin qk_temporal_solver)
target_link_libraries(quikin qk_data)
target_link_libraries(quikin qk_variable)
target_link_libraries(quikin qk_grid)
target_link_libraries(quikin qk_basis)
