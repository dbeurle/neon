
file(GLOB_RECURSE NEON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}" "*.cpp")

# Remove the main file to make a library
list(REMOVE_ITEM NEON_SOURCES "neon.cpp")

add_library(neon ${NEON_SOURCES})

target_compile_options(neon PUBLIC $<$<CONFIG:RELEASE>:-O3 -DNDEBUG>)

add_dependencies(neon eigen3)
add_dependencies(neon range-v3)
add_dependencies(neon termcolor)

target_link_libraries(neon
                      ${Boost_LIBRARIES}
                      ${JSONCPP_LIBRARIES}
                      ${MUMPS_LIBRARY}
                      ${BLAS_LIBRARIES}
                      ${VTK_LIBRARIES}
                      ${Hwloc_LIBRARIES}
                      ${CUDA_LIBRARIES}
                      ${CUDA_CUBLAS_LIBRARIES}
                      ${CUDA_cusparse_LIBRARY}
                      pastix
                      scotch
                      scotcherr)
