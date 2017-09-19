cmake_minimum_required(VERSION 3.8)

project(neonfe C CXX)

set(CMAKE_CXX_STANDARD 17)
add_definitions("-Wall -Wno-sign-compare -Wno-unused-variable -Wno-uninitialized -Wno-inconsistent-missing-override -Wno-unknown-pragmas")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(ENABLE_CUDA "Enable CUDA linear solver" OFF)
option(ENABLE_NATIVE "Optimise for current architecture" OFF)
option(ENABLE_COVERAGE "Set compiler flag for coverage analysis" OFF)

if (ENABLE_CUDA)
    find_package(CUDA REQUIRED)
    add_definitions(-DENABLE_CUDA)
endif()

if (ENABLE_NATIVE)
    add_definitions("-march=native -mtune=native")
endif()

if (ENABLE_COVERAGE)
    include(Coverage)
endif()

set(EXTERNAL_PROJECTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")

add_subdirectory(${EXTERNAL_PROJECTS_DIR}/catch)
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/eigen3)
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/range-v3)
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/termcolor)

include_directories(${CATCH_INCLUDE_DIR}
                    ${RV3_INCLUDE_DIR}
                    ${EIGEN_INCLUDE_DIR}
                    ${TERMCOLOR_INCLUDE_DIR}
                    ${COMMON_INCLUDES})

find_package(VTK COMPONENTS vtkIOXML REQUIRED)
include(${VTK_USE_FILE})

find_package(Boost 1.40 COMPONENTS filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(MUMPS REQUIRED)
include_directories(${MUMPS_INCLUDE_DIR})

find_package(Jsoncpp REQUIRED)
include_directories(${JSONCPP_INCLUDE_DIR})

find_package(HwLoc REQUIRED)

include_directories("${CMAKE_SOURCE_DIR}/src")

include(CTest)
add_subdirectory(src)
add_subdirectory(tests)

add_executable(neonfe src/neon.cpp)

target_link_libraries(neonfe neon)

install(TARGETS neonfe RUNTIME DESTINATION bin)
