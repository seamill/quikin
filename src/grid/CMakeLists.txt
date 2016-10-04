
add_library(qk_grid SHARED
    grid.cpp
    rectilinear.cpp
    structured.cpp
)

target_link_libraries(qk_grid qk)
target_link_libraries(qk_grid qk_basis)