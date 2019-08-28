cmake_minimum_required(VERSION 3.9)

project(neonfe VERSION 0.1 LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}"
      CACHE STRING "Choose the type of build."
      FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE
               PROPERTY STRINGS
                        "Debug"
                        "Release"
                        "MinSizeRel"
                        "RelWithDebInfo")
endif()

# Performance options
option(ENABLE_LTO "Enable link time optimisations" OFF)
option(ENABLE_NATIVE "Optimise for current architecture" OFF)
option(ENABLE_FAST_MATH "Compiler optimisations for floating point math" OFF)
# Metrics
option(ENABLE_TESTS "Enable the test suite of the unit tests examples" ON)
option(ENABLE_PROFILE "Set compiler flag for profiling" OFF)
option(ENABLE_COVERAGE "Set compiler flag for coverage analysis" OFF)
# Specialised accelerators
option(ENABLE_CUDA "Enable acceleratedCUDA solvers" OFF)
option(ENABLE_OPENCL "Enable OpenCL accelerated solvers" OFF)

set(EXTERNAL_PROJECTS_DIR "${PROJECT_SOURCE_DIR}/external")
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/json)
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/eigen3)
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/range-v3)
add_subdirectory(${EXTERNAL_PROJECTS_DIR}/termcolor)

find_package(ARPACK REQUIRED)
find_package(Pastix REQUIRED)
find_package(BLAS REQUIRED)
find_package(Boost 1.40 COMPONENTS filesystem REQUIRED)
find_package(MUMPS REQUIRED)
find_package(OpenMP REQUIRED)
find_package(TBB REQUIRED)
find_package(VTK COMPONENTS vtkIOXML REQUIRED)

if(ENABLE_OPENCL)
  find_package(OpenCL REQUIRED)
  find_package(ViennaCL CONFIG)
  add_definitions(-DENABLE_OPENCL)
endif()

if(ENABLE_CUDA)
  enable_language(CUDA)
  find_package(CUDA REQUIRED)
  add_definitions(-DENABLE_CUDA)
endif()

if(ENABLE_COVERAGE)
  include(Coverage)
endif()

if(ENABLE_TESTS)
  find_package(Catch2 2.3)

  if(NOT ${Catch2_FOUND})
    message(STATUS "System Catch2 not found.  Building...")
    add_subdirectory(${EXTERNAL_PROJECTS_DIR}/catch2)
  endif()

  # Required for integration tests
  include(FindPythonInterp)

  enable_testing()
  add_subdirectory(test)
  add_subdirectory(examples)
endif()

add_subdirectory(src)

add_executable(neonfe src/neon.cpp)

target_link_libraries(neonfe neon)

install(TARGETS neonfe RUNTIME DESTINATION bin)
