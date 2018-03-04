
foreach(test
             boundary
             constitutive
             dof_numbering
             gauss_quadrature
             linear_solvers
             eigen_solvers
             material
             material_coordinates
             matrix_system
             mesh
             nodal_coordinates
             tensor
             time_stepping)

    add_executable(${test} ${test}.cpp)
    add_dependencies(${test} catch neon)

    target_link_libraries(${test} LINK_PUBLIC neon)
    target_include_directories(${test} PUBLIC ${CMAKE_SOURCE_DIR}/src
                                              OpenMP::OpenMP_CXX
                                              ${CATCH_INCLUDE_DIR}
                                              ${EIGEN_INCLUDE_DIR}
                                              ${CATCH_INCLUDE_DIR}
                                              ${VTK_INCLUDE_DIRS}
                                              ${RV3_INCLUDE_DIR})

    set_target_properties(${test} PROPERTIES CXX_STANDARD 17
                                             CXX_STANDARD_REQUIRED YES
                                             CXX_EXTENSIONS NO
                                  COMPILE_FLAGS "-Wall")

    add_test(${test} ${CMAKE_CURRENT_BUILD_DIR}/${test})
endforeach()
