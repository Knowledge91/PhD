cmake_minimum_required(VERSION 3.7)

set(CMAKE_Fortran_COMPILER "/opt/local/bin/gfortran")

project(PhD CXX Fortran)
set(CMAKE_CXX_STANDARD 14)


# SET UP MINUIT2 http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/
file(GLOB ROOT_O
        "/Users/Knowledge/Developer/PhD/TauFitter/frameworks/Minuit2-5.34.14/src/*.o")
include_directories(program program_test frameworks/boost_1_64_0
        /Users/Knowledge/Developer/PhD/TauFitter/frameworks/Minuit2-5.34.14/inc)
# Minuit needs OpenMP, otherways missing symbols
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCE_FILES main.cpp)
add_executable(program_run ${SOURCE_FILES})

add_subdirectory(program)
add_subdirectory(program_test)

target_link_libraries(program_run program
        ${ROOT_O})
