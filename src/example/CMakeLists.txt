
add_subdirectory(advection)
add_subdirectory(euler)
add_subdirectory(maxwell)
add_subdirectory(plasma)

link_directories(advection euler maxwell plasma)