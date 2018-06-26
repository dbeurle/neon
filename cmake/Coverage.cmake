
add_custom_target(coverage
                  COMMAND lcov --capture --directory ${CMAKE_BINARY_DIR} --output-file coverage.info
                  # Remove the external libraries to get coverage for neon source only
                  COMMAND lcov --remove coverage.info '/usr/*'
                                                      '${CMAKE_BINARY_DIR}/catch*'
                                                      '${CMAKE_BINARY_DIR}/eigen3*'
                                                      '${CMAKE_BINARY_DIR}/json*'
                                                      '${CMAKE_BINARY_DIR}/range-v3*'
                                                      '${CMAKE_BINARY_DIR}/termcolor*'
                                                      '${CMAKE_SOURCE_DIR}/tests*' -o coverage.info)
                  # COMMAND genhtml coverage.info --output-directory coverage_output
                  # COMMAND firefox coverage_output/index.html)
