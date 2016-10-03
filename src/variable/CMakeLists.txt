
add_library(qk_variable SHARED
    variable.cpp
    variable_id.cpp
    variable_manager.cpp
)
    
target_link_libraries(qk_variable qk)
target_link_libraries(qk_variable qk_grid)
target_link_libraries(qk_variable qk_basis)
target_link_libraries(qk_variable qk_data)