
include(ExternalProject)

ExternalProject_Add(json
		GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_SHALLOW 1
        PREFIX ${CMAKE_BINARY_DIR}/json
		# Disable the install and update steps
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
		INSTALL_COMMAND ""
)

ExternalProject_Get_Property(json source_dir)

set(JSON_INCLUDE_DIR "${source_dir}" CACHE INTERNAL "json include directory")
