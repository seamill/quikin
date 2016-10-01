
add_library(qk_data SHARED
    dataset.cpp
    datachunk.cpp
    extended_datachunk.cpp
    shared_datachunk.cpp
)

target_link_libraries(qk_data qk)