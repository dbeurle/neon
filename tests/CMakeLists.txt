
foreach(test
             nodal_coordinates
             material_coordinates
             gauss_quadrature
             AdaptiveIncrementTest
             BoundaryTest
             ConstitutiveTest
             DofAllocation
             GeometryTest
             LinearSolverTest
             MaterialTest
             MatrixSystemTest
             ModuleTest
             TensorTest
             )
    add_executable(${test} ${test}.cpp)
    add_dependencies(${test} catch)
    target_link_libraries(${test} LINK_PUBLIC neon)
    add_test(${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
endforeach()
