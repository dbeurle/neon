cmake_minimum_required(VERSION 3.0)
project(blaze_builder CXX)
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    blaze
    PREFIX ${CMAKE_BINARY_DIR}/blaze
    GIT_REPOSITORY https://bitbucket.org/blaze-lib/blaze.git
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

# Expose required variable (CATCH_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(blaze source_dir)
set(BLAZE_INCLUDE_DIR ${source_dir}/blaze CACHE INTERNAL "Path to include folder for Blaze")
