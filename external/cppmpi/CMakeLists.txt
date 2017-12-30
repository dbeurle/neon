cmake_minimum_required(VERSION 3.0)
project(cppmpi_builder CXX)
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    cppmpi
    PREFIX ${CMAKE_BINARY_DIR}/cppmpi
    GIT_REPOSITORY https://github.com/dbeurle/cppmpi.git
    TIMEOUT 10
    # UPDATE_COMMAND ${HG_EXECUTABLE} pull
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
   )

ExternalProject_Get_Property(cppmpi source_dir)
set(CPPMPI_INCLUDE_DIR ${source_dir}/include CACHE INTERNAL "Path to include folder for cppmpi")
