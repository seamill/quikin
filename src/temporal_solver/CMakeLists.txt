
add_library(qk_temporal_solver SHARED
    ssprk2.cpp
)
    
target_link_libraries(qk_temporal_solver qk)
target_link_libraries(qk_temporal_solver qk_data)