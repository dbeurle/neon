
include(ExternalProject)

ExternalProject_Add(json

    URL https://github.com/nlohmann/json/releases/download/v3.1.2/include.zip
    URL_HASH SHA256=495362ee1b9d03d9526ba9ccf1b4a9c37691abe3a642ddbced13e5778c16660c
    PREFIX ${CMAKE_BINARY_DIR}/json
    # Disable the install and update steps
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
INSTALL_COMMAND ""
)

ExternalProject_Get_Property(json source_dir)

set(JSON_INCLUDE_DIR "${source_dir}" CACHE INTERNAL "json include directory")
