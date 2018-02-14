
foreach(test
             boundary
             constitutive
             dof_numbering
             gauss_quadrature
             linear_solvers
             material
             material_coordinates
             matrix_system
             mesh
             nodal_coordinates
             tensor
             time_stepping
             )
    add_executable(${test} ${test}.cpp)
    add_dependencies(${test} catch)
    target_link_libraries(${test} LINK_PUBLIC neon)
    add_test(${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
endforeach()
