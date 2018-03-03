
file(GLOB_RECURSE NEON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}" "*.cpp")

# Remove the main file to make a library
list(REMOVE_ITEM NEON_SOURCES "neon.cpp")

add_library(neon ${NEON_SOURCES})

target_compile_options(neon PUBLIC $<$<CONFIG:RELEASE>:-O3 -DNDEBUG>)

set_target_properties(neon PROPERTIES CXX_STANDARD 17
                                      CXX_STANDARD_REQUIRED YES
                                      CXX_EXTENSIONS NO
                           COMPILE_FLAGS "-Wall -Wno-unused-variable ${OpenMP_CXX_FLAGS}")

include(CheckIPOSupported)
check_ipo_supported(RESULT is_ipo_supported)
if(is_ipo_supported)
   set_target_properties(neon PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

if (ENABLE_NATIVE)
   set_target_properties(neon COMPILE_FLAGS "-march=native -mtune=native")
endif()

target_include_directories(neon PUBLIC ${CMAKE_SOURCE_DIR}/src
                                       ${Boost_INCLUDE_DIR}
                                       ${JSON_INCLUDE_DIR}
                                       ${MUMPS_INCLUDE_DIR}
                                       ${EIGEN_INCLUDE_DIR}
                                       ${BLAZE_INCLUDE_DIR}
                                       ${VTK_INCLUDE_DIRS}
                                       ${RV3_INCLUDE_DIR}
                                       ${TERMCOLOR_INCLUDE_DIR})

add_dependencies(neon eigen3 range-v3 termcolor json)

if (ENABLE_CUDA)
    add_library(dmatvec STATIC ${CMAKE_SOURCE_DIR}/src/solver/linear/dmatrix_vector_product.cu)
    target_link_libraries(neon dmatvec)
    add_dependencies(neon dmatvec)
endif()

target_link_libraries(neon
                      OpenMP::OpenMP_CXX
                      ${Boost_LIBRARIES}
                      ${MUMPS_LIBRARY}
                      ${BLAS_LIBRARIES}
                      ${VTK_LIBRARIES}
                      ${Hwloc_LIBRARIES}
                      ${TBB_LIBRARIES}
                      pastix
                      jsoncpp
                      scotch
                      scotcherr
                      ${ARPACK_LIBRARY})
