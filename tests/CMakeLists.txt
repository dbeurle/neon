
foreach(test GeometryTest
             MaterialTest
             ConstitutiveTest
             QuadratureTest
             MatrixSystemTest
             AdaptiveIncrementTest
             BoundaryTest
             TensorTest
             LinearSolverTest)
    add_executable(${test} ${test}.cpp)
    add_dependencies(${test} catch)
    target_link_libraries(${test} LINK_PUBLIC neon)
    add_test(${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
endforeach()
