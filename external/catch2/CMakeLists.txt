
include(ExternalProject)

ExternalProject_Add(
    catch2-external
    PREFIX ${CMAKE_BINARY_DIR}/external/catch2
    URL https://github.com/catchorg/Catch2/archive/v2.4.2.tar.gz
    # Disable the install and update steps
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(catch2-external source_dir)

add_library(Catch2 INTERFACE)

add_dependencies(Catch2 catch2-external)

target_include_directories(Catch2 INTERFACE ${source_dir}/single_include)

add_library(Catch2::Catch2 ALIAS Catch2)
