
file(GLOB_RECURSE NEON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}" "*.cpp")

# Remove the main file to make a library
list(REMOVE_ITEM NEON_SOURCES "neon.cpp")

cuda_add_library(dmatvec solver/linear/dmatrix_vector_product.cu STATIC)

target_link_libraries(dmatvec ${CUDA_LIBRARIES})

add_library(neon ${NEON_SOURCES})

target_compile_options(neon PUBLIC $<$<CONFIG:RELEASE>:-O3 -DNDEBUG>)

add_dependencies(neon eigen3 range-v3 termcolor blaze dmatvec)

target_link_libraries(neon
                      ${Boost_LIBRARIES}
                      ${JSONCPP_LIBRARY}
                      ${MUMPS_LIBRARY}
                      ${BLAS_LIBRARIES}
                      ${VTK_LIBRARIES}
                      ${Hwloc_LIBRARIES}
                      ${CUDA_LIBRARIES}
                      ${CUDA_CUBLAS_LIBRARIES}
                      ${CUDA_cusparse_LIBRARY}
                      dmatvec
                      pastix
                      scotch
                      scotcherr)
