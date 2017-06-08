project(program_tests)

include_directories(basic_tests)

add_subdirectory(lib/gtest-1.7.0)
add_subdirectory(basic_tests)