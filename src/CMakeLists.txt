
add_subdirectory(lib)
add_subdirectory(solvers)
add_subdirectory(data)
add_subdirectory(grid)

add_executable(quikin main.cpp)

link_directories(lib solvers data grid)

target_link_libraries(quikin qk)
target_link_libraries(quikin qk_grid)
target_link_libraries(quikin qk_data)
target_link_libraries(quikin qk_solvers)