cmake_minimum_required(VERSION 3.6)
project(quikin)

set(quikin_VERSION_MAJOR 0)
set(quikin_VERSION_MINOR 0)

set(CMAKE_MACOSX_RPATH 1)

add_definitions(-D_QK_RANGE_CHECK_=1)

include_directories(src)
add_subdirectory(src)

enable_testing()

include_directories(test)
add_subdirectory(test)


