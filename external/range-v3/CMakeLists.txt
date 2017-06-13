cmake_minimum_required(VERSION 3.0)
project(range-v3_builder CXX)
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    range-v3
    PREFIX ${CMAKE_BINARY_DIR}/range-v3
    GIT_REPOSITORY https://github.com/ericniebler/range-v3.git
    GIT_SHALLOW 1
    GIT_PROGRESS 1
    TIMEOUT 10
    # UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
   )

# Expose required variable (RANGE-V3_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(range-v3 source_dir)
set(RV3_INCLUDE_DIR ${source_dir}/include CACHE INTERNAL "Path to include folder for Range-v3")
