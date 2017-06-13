cmake_minimum_required(VERSION 3.7)

set(CMAKE_Fortran_COMPILER "/opt/local/bin/gfortran")

project(PhD CXX Fortran)

set(CMAKE_CXX_STANDARD 14)

include_directories(program program_test frameworks/boost_1_64_0)

set(SOURCE_FILES main.cpp)
add_executable(program_run ${SOURCE_FILES})

add_subdirectory(program)
add_subdirectory(program_test)

target_link_libraries(program_run program)