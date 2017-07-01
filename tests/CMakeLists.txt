
add_executable(GeometryUnitTest GeometryTest.cpp)
target_link_libraries(GeometryUnitTest LINK_PUBLIC neon)
add_test(GeometryUnitTest ${CMAKE_CURRENT_BUILD_DIR}/GeometryUnitTest)

add_executable(MaterialTest MaterialTest.cpp)
target_link_libraries(MaterialTest LINK_PUBLIC neon)
add_test(MaterialTest ${CMAKE_CURRENT_BUILD_DIR}/MaterialTest)

add_executable(ConstitutiveTest ConstitutiveTest.cpp)
target_link_libraries(ConstitutiveTest LINK_PUBLIC neon)
add_test(ConstitutiveTest ${CMAKE_CURRENT_BUILD_DIR}/ConstitutiveTest)

add_executable(QuadratureTest QuadratureTest.cpp)
target_link_libraries(QuadratureTest LINK_PUBLIC neon)
add_test(QuadratureTest ${CMAKE_CURRENT_BUILD_DIR}/QuadratureTest)

add_executable(SystemTest MatrixSystemTest.cpp)
target_link_libraries(SystemTest LINK_PUBLIC neon)
add_test(SystemTest ${CMAKE_CURRENT_BUILD_DIR}/SystemTest)

# add_executable(BasicMeshTest BasicMeshTest.cpp)
# target_link_libraries(BasicMeshTest LINK_PUBLIC neon)
# add_test(BasicMeshTest ${CMAKE_CURRENT_BUILD_DIR}/BasicMeshTest)
#
# add_executable(femSolidTest femSolidTest.cpp)
# target_link_libraries(femSolidTest LINK_PUBLIC neon)
# add_test(femSolidTest ${CMAKE_CURRENT_BUILD_DIR}/femSolidTest)
