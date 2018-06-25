
if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(FATAL_ERROR "Coverage only available when compiling with GCC.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

add_custom_target(coverage
                  COMMAND ctest
                  COMMAND lcov --capture --directory ${CMAKE_BINARY_DIR} --output-file coverage.info
                  # Remove the external libraries to get coverage for neon source only
                  COMMAND lcov --remove coverage.info '/usr/*'
                                                      '${CMAKE_BINARY_DIR}/catch*'
                                                      '${CMAKE_BINARY_DIR}/eigen3*'
                                                      '${CMAKE_BINARY_DIR}/json*'
                                                      '${CMAKE_BINARY_DIR}/range-v3*'
                                                      '${CMAKE_BINARY_DIR}/termcolor*'
                                                      -o coverage.info)
                  # COMMAND genhtml coverage.info --output-directory coverage_output
                  # COMMAND firefox coverage_output/index.html)
