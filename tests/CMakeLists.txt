
add_executable(GeometryUnitTest GeometryTest.cpp)
add_dependencies(GeometryUnitTest catch)
target_link_libraries(GeometryUnitTest LINK_PUBLIC neon)
add_test(GeometryUnitTest ${CMAKE_CURRENT_BUILD_DIR}/GeometryUnitTest)

add_executable(MaterialTest MaterialTest.cpp)
add_dependencies(MaterialTest catch)
target_link_libraries(MaterialTest LINK_PUBLIC neon)
add_test(MaterialTest ${CMAKE_CURRENT_BUILD_DIR}/MaterialTest)

add_executable(ConstitutiveTest ConstitutiveTest.cpp)
add_dependencies(ConstitutiveTest catch)
target_link_libraries(ConstitutiveTest LINK_PUBLIC neon)
add_test(ConstitutiveTest ${CMAKE_CURRENT_BUILD_DIR}/ConstitutiveTest)

add_executable(QuadratureTest QuadratureTest.cpp)
add_dependencies(QuadratureTest catch)
target_link_libraries(QuadratureTest LINK_PUBLIC neon)
add_test(QuadratureTest ${CMAKE_CURRENT_BUILD_DIR}/QuadratureTest)

add_executable(SystemTest MatrixSystemTest.cpp)
add_dependencies(SystemTest catch)
target_link_libraries(SystemTest LINK_PUBLIC neon)
add_test(SystemTest ${CMAKE_CURRENT_BUILD_DIR}/SystemTest)

add_executable(AdaptiveIncrement AdaptiveIncrementTest.cpp)
add_dependencies(AdaptiveIncrement catch)
target_link_libraries(AdaptiveIncrement LINK_PUBLIC neon)
add_test(AdaptiveIncrement ${CMAKE_CURRENT_BUILD_DIR}/AdaptiveIncrement)

add_executable(BoundaryTest BoundaryTest.cpp)
add_dependencies(BoundaryTest catch)
target_link_libraries(BoundaryTest LINK_PUBLIC neon)
add_test(BoundaryTest ${CMAKE_CURRENT_BUILD_DIR}/BoundaryTest)
