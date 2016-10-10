cmake_minimum_required(VERSION 3.6)
project(quikin)

set(quikin_VERSION_MAJOR 0)
set(quikin_VERSION_MINOR 0)

set(CMAKE_MACOSX_RPATH 1)

option(USE_OPENMP "Use OpenMP" OFF)

#if(OPENMP_FOUND)
#    message("OPENMP FOUND")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#else(not OPENMP_FOUND)
#    message("OPENMP NOT FOUND")
#endif(OPENMP_FOUND)

if(USE_OPENMP)
    message("-- Using OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif(USE_OPENMP)


if(CMAKE_BUILD_TYPE MATCHES DEBUG)
    add_definitions(-D_QK_RANGE_CHECK_=1)
elseif(CMAKE_BUILD_TYPE MATCHES RELEASE)
    add_compile_flag("-O3")
endif(CMAKE_BUILD_TYPE MATCHES DEBUG)

include_directories(src)
add_subdirectory(src)

enable_testing()

include_directories(test)
add_subdirectory(test)


