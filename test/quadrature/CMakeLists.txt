
set(test_names minimum_degree line triangle quadrilateral hexahedron prism pyramid sphere tetrahedron)

foreach(test_name IN LISTS test_names)

    add_executable(quadrature_${test_name}_test ${test_name}.cpp)
    add_dependencies(quadrature_${test_name}_test neon catch_main cube_fixture)
    target_link_libraries(quadrature_${test_name}_test PRIVATE catch_main
                                                    cube_fixture
                                                    neon
                                                    Catch2::Catch2
                                                    OpenMP::OpenMP_CXX)

    target_include_directories(quadrature_${test_name}_test PUBLIC ${CMAKE_SOURCE_DIR}/src
                                                        ${EIGEN_INCLUDE_DIR})

    set_target_properties(quadrature_${test_name}_test PROPERTIES CXX_STANDARD 17
                                                       CXX_STANDARD_REQUIRED YES
                                                       CXX_EXTENSIONS NO
                                            COMPILE_FLAGS "-Wall -Wextra")

    if(ENABLE_COVERAGE)
        set_target_properties(quadrature_${test_name}_test PROPERTIES COMPILE_FLAGS "--coverage")
        target_link_libraries(quadrature_${test_name}_test PRIVATE --coverage)
    endif()

    add_test(quadrature_${test_name}_test ${CMAKE_CURRENT_BINARY_DIR}/quadrature_${test_name}_test)

endforeach()
