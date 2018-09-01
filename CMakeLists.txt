
cmake_minimum_required(VERSION 3.9)

project(neonfe VERSION 0.1 LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(ENABLE_LTO "Enable link time optimisations" OFF)
option(ENABLE_CUDA "Enable CUDA linear solver" OFF)
option(ENABLE_OCL "Enable OpenCL accelerated linear solver" OFF)
option(ENABLE_TESTS "Enable the test suite of the unit tests examples" ON)
option(ENABLE_NATIVE "Optimise for current architecture" OFF)
option(ENABLE_PROFILE "Set compiler flag for profiling" OFF)
option(ENABLE_COVERAGE "Set compiler flag for coverage analysis" OFF)

set(EXTERNAL_PROJECTS_DIR "${PROJECT_SOURCE_DIR}/external")
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/json)
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/eigen3)
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/range-v3)
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/termcolor)

find_package(ARPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(Boost 1.40 COMPONENTS filesystem REQUIRED)
find_package(Catch2 REQUIRED)
find_package(HwLoc REQUIRED)
find_package(MUMPS REQUIRED)
find_package(OpenMP REQUIRED)
find_package(TBB REQUIRED)
find_package(VTK COMPONENTS vtkIOXML REQUIRED)

find_package(OpenCL REQUIRED)
find_package(ViennaCL CONFIG)
add_definitions(-DENABLE_OCL)

if (ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDA REQUIRED)
    add_definitions(-DENABLE_CUDA)
endif()

if (ENABLE_COVERAGE)
    include(Coverage)
endif()

if (ENABLE_TESTS)
    # Required for integration tests
    include(FindPythonInterp)

    include(CTest)
    add_subdirectory(test)
    add_subdirectory(examples)
endif()

add_subdirectory(src)

add_executable(neonfe src/neon.cpp)

target_link_libraries(neonfe neon)

install(TARGETS neonfe RUNTIME DESTINATION bin)
