
# Find PASTIX
#
# Find the PASTIX includes and library
#
# if you need to add a custom library search path, do it via via CMAKE_PREFIX_PATH
#
# This module defines
#  PASTIX_INCLUDE_DIRS, where to find header, etc.
#  PASTIX_LIBRARIES, the libraries needed to use PASTIX.
#  PASTIX_FOUND, If false, do not try to use PASTIX.
#  PASTIX_INCLUDE_PREFIX, include prefix for PASTIX

# only look in default directories
find_path(PASTIX_INCLUDE_DIR NAMES pastix/pastix.h pastix.h DOC "PaStiX include dir")

find_library(PASTIX_LIBRARY NAMES pastix DOC "PaStiX library")

set(PASTIX_INCLUDE_DIRS ${PASTIX_INCLUDE_DIR})
set(PASTIX_LIBRARIES ${PASTIX_LIBRARY})

# find PASTIX_INCLUDE_PREFIX
find_path(PASTIX_INCLUDE_PREFIX NAMES pastix.h PATH_SUFFIXES pastix)
if (${PASTIX_INCLUDE_PREFIX} MATCHES "pastix")
	set(PASTIX_INCLUDE_DIR "${PASTIX_INCLUDE_DIR}/pastix")
endif()

# handle the QUIETLY and REQUIRED arguments and set PASTIX_FOUND to TRUE
# if all listed variables are TRUE, hide their existence from configuration view
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(PASTIX DEFAULT_MSG PASTIX_INCLUDE_DIR PASTIX_LIBRARY)

mark_as_advanced(PASTIX_INCLUDE_DIR PASTIX_LIBRARY)
