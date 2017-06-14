
add_executable(GeometryUnitTest GeometryTest.cpp)
target_link_libraries(GeometryUnitTest LINK_PUBLIC neon)
add_test(GeometryUnitTest ${CMAKE_CURRENT_BUILD_DIR}/GeometryUnitTest)

add_executable(MaterialTest MaterialTest.cpp)
target_link_libraries(MaterialTest LINK_PUBLIC neon)
add_test(MaterialTest ${CMAKE_CURRENT_BUILD_DIR}/MaterialTest)

add_executable(ConstitutiveTest ConstitutiveTest.cpp)
target_link_libraries(ConstitutiveTest LINK_PUBLIC neon)
add_test(ConstitutiveTest ${CMAKE_CURRENT_BUILD_DIR}/ConstitutiveTest)

# add_executable(NumericalQuadrature QuadratureUnitTest.cpp)
# target_link_libraries(NumericalQuadrature LINK_PUBLIC neon)
# add_test(NumericalQuadrature ${CMAKE_CURRENT_BUILD_DIR}/NumericalQuadrature)
#
# add_executable(TensorTest TensorUnitTest.cpp)
# target_link_libraries(TensorTest LINK_PUBLIC neon)
# add_test(TensorTest ${CMAKE_CURRENT_BUILD_DIR}/TensorTest)
#
# add_executable(BasicMeshTest BasicMeshTest.cpp)
# target_link_libraries(BasicMeshTest LINK_PUBLIC neon)
# add_test(BasicMeshTest ${CMAKE_CURRENT_BUILD_DIR}/BasicMeshTest)
#
# add_executable(femSolidTest femSolidTest.cpp)
# target_link_libraries(femSolidTest LINK_PUBLIC neon)
# add_test(femSolidTest ${CMAKE_CURRENT_BUILD_DIR}/femSolidTest)
