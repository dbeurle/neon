
# Copy across the mesh asset
configure_file("cube.mesh"
               "${CMAKE_CURRENT_BINARY_DIR}/cube.mesh"
               COPYONLY)

# Run the simulation
add_test(NAME neohooke_cube
         COMMAND ${CMAKE_BINARY_DIR}/neonfe ${CMAKE_CURRENT_SOURCE_DIR}/cube.json)

# Verify using known solution
add_test(NAME neohooke_cube_verify
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/verify.py visualisation/cube_1.vtu displacement)

set_tests_properties(neohooke_cube_verify PROPERTIES DEPENDS neohooke_cube)
