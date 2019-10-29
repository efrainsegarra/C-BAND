# CMakeLists for Calibration BAND Project
project(C-BAND)

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)
find_package(LZ4 REQUIRED)
#set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")

# set globals for use of all packages
file( GLOB HEADERS ${CMAKE_SOURCE_DIR}/*.h )

# Add directories
add_subdirectory(Hipo4)
add_subdirectory(Banks4)
add_subdirectory(Scripts_Spring19_Hipo4)


