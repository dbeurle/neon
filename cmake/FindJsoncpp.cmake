
#Copyright (C) 2011-2016 Peter Spiess-Knafl

#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to
#use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
#the Software, and to permit persons to whom the Software is furnished to do so,
#subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
#PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
#TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
#OR OTHER DEALINGS IN THE SOFTWARE.

# Taken from Cinemast, from user debris.

# Find jsoncpp
#
# Find the jsoncpp includes and library
#
# if you nee to add a custom library search path, do it via via CMAKE_PREFIX_PATH
#
# This module defines
#  JSONCPP_INCLUDE_DIRS, where to find header, etc.
#  JSONCPP_LIBRARIES, the libraries needed to use jsoncpp.
#  JSONCPP_FOUND, If false, do not try to use jsoncpp.
#  JSONCPP_INCLUDE_PREFIX, include prefix for jsoncpp

# Only look in default directories
find_path(JSONCPP_INCLUDE_DIRS
	      json/json.h
          PATH_SUFFIXES jsoncpp)

find_library(JSONCPP_LIBRARY
	         NAMES jsoncpp
	         DOC "jsoncpp library")

set(JSONCPP_LIBRARIES ${JSONCPP_LIBRARY})
set(JSONCPP_INCLUDE_DIRS ${JSONCPP_INCLUDE_DIRS})

# find JSONCPP_INCLUDE_PREFIX
find_path(JSONCPP_INCLUDE_PREFIX
	      json/json.h
          PATH_SUFFIXES jsoncpp)

if (${JSONCPP_INCLUDE_PREFIX} MATCHES "jsoncpp")
	set(JSONCPP_INCLUDE_PREFIX "jsoncpp/")
else()
	set(JSONCPP_INCLUDE_PREFIX "json/")
endif()

# handle the QUIETLY and REQUIRED arguments and set JSONCPP_FOUND to TRUE
# if all listed variables are TRUE, hide their existence from configuration view
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(jsoncpp DEFAULT_MSG
JSONCPP_INCLUDE_DIRS JSONCPP_LIBRARY)
mark_as_advanced (JSONCPP_INCLUDE_DIRS JSONCPP_LIBRARY)
