
# cmake_minimum_required(VERSION 3.0)
# project(jsoncpp_builder CXX)
include(ExternalProject)

ExternalProject_Add(json-cpp
		GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
        GIT_SHALLOW 1
		CMAKE_ARGS
			-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
			-DJSONCPP_WITH_TESTS:BOOL=OFF
            -DBUILD_SHARED_LIBS:BOOL=ON
            -DBUILD_STATIC_LIBS:BOOL=OFF
		# Disable the install and update steps
        UPDATE_COMMAND ""
		INSTALL_COMMAND ""
	)

ExternalProject_Get_Property(json-cpp source_dir)
ExternalProject_Get_Property(json-cpp binary_dir)

set(JSONCPP_INCLUDE_DIR "${source_dir}/include" CACHE INTERNAL "JSONCPP include")
set(JSONCPP_LIBRARY_DIR "${binary_dir}/src/lib_json")
set(JSONCPP_LIBRARY jsoncpp CACHE INTERNAL "JSONCPP library")

add_library(jsoncpp SHARED IMPORTED)

set_property(TARGET jsoncpp PROPERTY IMPORTED_LOCATION ${JSONCPP_LIBRARY_DIR}/jsoncpp.so)

add_dependencies(jsoncpp json-cpp)

link_directories(${JSONCPP_LIBRARY_DIR})
