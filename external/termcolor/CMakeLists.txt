
include(ExternalProject)

find_package(Git REQUIRED)

ExternalProject_Add(
    termcolor
    PREFIX ${CMAKE_BINARY_DIR}/termcolor
    GIT_REPOSITORY https://github.com/ikalnytskyi/termcolor.git
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

# Expose required variable (TERMCOLOR_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(termcolor source_dir)
set(TERMCOLOR_INCLUDE_DIR ${source_dir}/include CACHE INTERNAL "Path to include folder for termcolor")
