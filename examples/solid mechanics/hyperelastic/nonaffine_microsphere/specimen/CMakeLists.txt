
# Copy across the required input and mesh files
configure_file("ElastomerTestSpecimen.json"
               "${CMAKE_CURRENT_BINARY_DIR}/ElastomerTestSpecimen.json"
               COPYONLY)

configure_file("ElastomerTestSpecimen.mesh"
              "${CMAKE_CURRENT_BINARY_DIR}/ElastomerTestSpecimen.mesh"
              COPYONLY)

add_test(nonaffine_microsphere_specimen neonfe ElastomerTestSpecimen.json)
