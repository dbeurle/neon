
# Copy across the required input and mesh files
configure_file("specimen.json"
               "${CMAKE_CURRENT_BINARY_DIR}/specimen.json"
               COPYONLY)

configure_file("specimen.mesh"
              "${CMAKE_CURRENT_BINARY_DIR}/specimen.mesh"
              COPYONLY)

add_test(gaussian_affine_microsphere_ageing_specimen neonfe specimen.json)
