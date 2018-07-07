
# Copy across the required input and mesh files
configure_file("cube.json"
               "${CMAKE_CURRENT_BINARY_DIR}/cube.json"
               COPYONLY)

configure_file("cube.mesh"
               "${CMAKE_CURRENT_BINARY_DIR}/cube.mesh"
               COPYONLY)

configure_file("verify.py"
               "${CMAKE_CURRENT_BINARY_DIR}/verify.py"
               COPYONLY)

add_test(neohooke_cube
         ${CMAKE_BINARY_DIR}/neonfe cube.json)

# add_test(neohooke_cube_verify
#          ${PYTHON_EXECUTABLE} check.py visualisation/cube_1.vtu displacement)
