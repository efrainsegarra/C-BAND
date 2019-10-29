# Setup file:
cmake_minimum_required (VERSION 2.6)
find_package(ROOT REQUIRED)


project (C-BAND)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_subdirectory( Hipo4 )
