
add_library(qk_solvers SHARED
    qksolver_advection.cpp
#    qksolver_boltzmannbgk_1d2v_shocktube_var.cpp
#    qksolver_boltzmannbgk_1d2v_shocktube.cpp
    qksolver.cpp
)
    
target_link_libraries(qk_solvers qk)
target_link_libraries(qk_solvers qk_grid)
target_link_libraries(qk_solvers qk_data)