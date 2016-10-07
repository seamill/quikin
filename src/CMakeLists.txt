
add_subdirectory(lib)
add_subdirectory(solver)
add_subdirectory(spatial_solver)
add_subdirectory(temporal_solver)
add_subdirectory(data)
add_subdirectory(variable)
add_subdirectory(grid)
add_subdirectory(basis)
add_subdirectory(example)

add_executable(quikin main.cpp)

#link_directories(lib solvers spatial_solver temporal_solver data variable grid basis)

target_link_libraries(quikin qk_example_advection)
target_link_libraries(quikin qk_example_euler)
target_link_libraries(quikin qk_example_maxwell)
target_link_libraries(quikin qk_example_plasma)
