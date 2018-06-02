
include(ExternalProject)

ExternalProject_Add(
    range-v3
    URL https://github.com/ericniebler/range-v3/archive/0.3.6.tar.gz
    PREFIX ${CMAKE_BINARY_DIR}/range-v3
    # Disable the install and update steps
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(range-v3 source_dir)

set(RV3_INCLUDE_DIR "${source_dir}/include" CACHE INTERNAL "range-v3 include directory")
