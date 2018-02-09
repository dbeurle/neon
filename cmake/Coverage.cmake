
if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(FATAL_ERROR "Coverage only available when compiling with GCC.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

add_custom_target(coverage
                  COMMAND ${CMAKE_MAKE_PROGRAM} -j4
                  COMMAND ctest
                  COMMAND lcov --capture --directory ${CMAKE_BINARY_DIR} --output-file coverage.info
                  # Remove the external libraries to get coverage for source only
                  COMMAND lcov --remove coverage.info '/usr/*' -o coverage.info
                  COMMAND lcov --remove coverage.info '/build/blaze/*' -o coverage.info
                  COMMAND lcov --remove coverage.info '/build/catch/*' -o coverage.info
                  COMMAND lcov --remove coverage.info '/build/eigen3/*' -o coverage.info
                  COMMAND lcov --remove coverage.info '/build/json/*' -o coverage.info
                  COMMAND lcov --remove coverage.info '/build/range-v3/*' -o coverage.info
                  COMMAND lcov --remove coverage.info '/build/termcolor/*' -o coverage.info
                  # Generate
                  COMMAND genhtml coverage.info --output-directory coverage_output
                  COMMAND xdg-open coverage_output/index.html)
