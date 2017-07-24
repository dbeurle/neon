
file(GLOB_RECURSE NEON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}" "*.cpp")

# Remove the main file to make a library
list(REMOVE_ITEM NEON_SOURCES "neon.cpp")

add_library(neon ${NEON_SOURCES})

target_link_libraries(neon
                      ${Boost_LIBRARIES}
                      ${CMAKE_THREAD_LIBS_INIT}
                      ${JSONCPP_LIBRARIES}
                      ${MUMPS_LIBRARY}
                      ${BLAS_LIBRARIES}
                      ${VTK_LIBRARIES}
                      pastix
                      metis
                      scotch
                      scotcherr
                      hwloc)
