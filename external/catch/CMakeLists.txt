
include(ExternalProject)

ExternalProject_Add(
    catch
    URL https://github.com/catchorg/Catch2/releases/download/v2.2.1/catch.hpp
    DOWNLOAD_NO_EXTRACT 1
    PREFIX ${CMAKE_BINARY_DIR}/catch
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

# Expose required variable (CATCH_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(catch DOWNLOAD_DIR)
set(CATCH_INCLUDE_DIR ${DOWNLOAD_DIR} CACHE INTERNAL "Path to include folder for Catch")
