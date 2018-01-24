
foreach(test
             adaptive_increment
             dof_numbering
             gauss_quadrature
             material_coordinates
             nodal_coordinates
             tensor
             linear_solvers
             material
             boundary
             constitutive
             mesh
             # MatrixSystemTest
             )
    add_executable(${test} ${test}.cpp)
    add_dependencies(${test} catch)
    target_link_libraries(${test} LINK_PUBLIC neon)
    add_test(${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
endforeach()
