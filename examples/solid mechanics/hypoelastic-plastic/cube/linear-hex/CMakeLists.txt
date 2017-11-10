
# Copy across the required input and mesh files
configure_file("cube.json"
               "${CMAKE_CURRENT_BINARY_DIR}/cube.json"
               COPYONLY)

configure_file("cube.mesh"
               "${CMAKE_CURRENT_BINARY_DIR}/cube.mesh"
               COPYONLY)

add_test(J2PlasticityHex neonfe cube.json)
