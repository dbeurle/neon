
add_library(catch_main STATIC test_main.cpp)
target_link_libraries(catch_main PUBLIC Catch2::Catch)

add_library(cube_fixture STATIC fixtures/cube_mesh.cpp)

set(test_names boundary
               constitutive
               compute_device
               dof_numbering
               gauss_quadrature
               geometry_profile
               linear_beam_theory
               linear_solvers
               eigenvalue_solvers
               svd_solvers
               jacobian_determinant
               material
               material_coordinates
               matrix_system
               mesh
               nodal_coordinates
               nodal_variables
               sequence
               tensor
               time_stepping
               time_integrators
               trapezoidal)

foreach(test_name IN LISTS test_names)

    add_executable(${test_name}_test ${test_name}.cpp)
    add_dependencies(${test_name}_test neon catch_main cube_fixture)
    target_link_libraries(${test_name}_test PRIVATE catch_main
                                               cube_fixture
                                               neon
                                               Catch2::Catch
                                               OpenMP::OpenMP_CXX)

    target_include_directories(${test_name}_test PUBLIC ${CMAKE_SOURCE_DIR}/src
                                                   fixtures
                                                   ${EIGEN_INCLUDE_DIR}
                                                   ${VTK_INCLUDE_DIRS}
                                                   ${RV3_INCLUDE_DIR})

    set_target_properties(${test_name}_test PROPERTIES CXX_STANDARD 17
                                                  CXX_STANDARD_REQUIRED YES
                                                  CXX_EXTENSIONS NO
                                       COMPILE_FLAGS "-Wall")

    if(ENABLE_COVERAGE)
        set_target_properties(${test_name}_test PROPERTIES COMPILE_FLAGS "--coverage")
        target_link_libraries(${test_name}_test PRIVATE --coverage)
    endif()

    add_test(${test_name}_test ${CMAKE_CURRENT_BINARY_DIR}/${test_name}_test)

endforeach()
