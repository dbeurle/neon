cmake_minimum_required(VERSION 3.0)

project(neonfe CXX)

set(EXT_PROJECTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")
add_subdirectory(${EXT_PROJECTS_DIR}/catch)
add_subdirectory(${EXT_PROJECTS_DIR}/eigen3)
add_subdirectory(${EXT_PROJECTS_DIR}/range-v3)
include_directories(${CATCH_INCLUDE_DIR} ${RV3_INCLUDE_DIR}
                    ${EIGEN_INCLUDE_DIR} ${COMMON_INCLUDES})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 17)
add_definitions("-Wall")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Boost 1.40 COMPONENTS system filesystem log thread timer REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(BLAS REQUIRED)

find_package(MUMPS REQUIRED)
include_directories(${MUMPS_INCLUDE_DIR})

find_package(Jsoncpp REQUIRED)

include_directories("${CMAKE_SOURCE_DIR}/src")

include(CTest)
add_subdirectory(src)
add_subdirectory(tests)

add_executable(neonfe src/neon.cpp)

add_dependencies(neon catch)
add_dependencies(neon eigen3)
add_dependencies(neon range-v3)

# Linker statements
target_link_libraries(neonfe neon)
