include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(SOURCE_FILES ConstantsFixture.cpp runAlpha_check.cpp numericalMethods_check.cpp AdlerFunction_check.cpp.cpp)

add_executable(runBasicTests ${SOURCE_FILES})

target_link_libraries(runBasicTests gtest gtest_main)
target_link_libraries(runBasicTests program)