# - Find a tvmet installation
# The following variables are set if tvmet is found.  If tvmet is not
# found, tvmet_FOUND is set to false.
#  tvmet_FOUND         - Set to true when tvmet is found.
#  tvmet_INCLUDE_DIRS  - Include directories for tvmet
#
# The following cache entries must be set by the user to locate tvmet:
#  tvmet_DIR  - This is the directory that contains sub-folder include 
#				which should contain another subfolder called tvmet which
#               contains the headers
#=============================================================================
# Copyright 2001-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

# Assume not found.
SET(tvmet_FOUND 0)

# Construct consitent error messages for use below.
SET(tvmet_DIR_DESCRIPTION "directory that contains sub-folder include which should contain another subfolder called tvmet which contains the headers")

SET(tvmet_DIR_MESSAGE "tvmet not found.  Set the tvmet_DIR cmake cache entry to the ${tvmet_DIR_DESCRIPTION}")

# Look for Usetvmet.cmake in build trees or under <prefix>/include/tvmet.
FIND_PATH(tvmet_DIR
   NAMES tvmet.h
   PATH_SUFFIXES "include/tvmet/"
   HINTS $ENV{LIB}
   DOC "The ${tvmet_DIR_DESCRIPTION}"
   )

   IF(tvmet_DIR)
	MESSAGE(STATUS "tvmet headers located under ${tvmet_DIR}")
	set(tvmet_FOUND 1)
	set(tvmet_INCLUDE_DIRS "${tvmet_DIR}/include")
   ENDIF(tvmet_DIR)
   
#-----------------------------------------------------------------------------
IF(tvmet_FOUND)
ELSE(tvmet_FOUND)
  # tvmet not found, explain to the user how to specify its location.
  IF(tvmet_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${tvmet_DIR_MESSAGE})
  ELSE(tvmet_FIND_REQUIRED)
    IF(NOT tvmet_FIND_QUIETLY)
      MESSAGE(STATUS ${tvmet_DIR_MESSAGE})
    ENDIF(NOT tvmet_FIND_QUIETLY)
  ENDIF(tvmet_FIND_REQUIRED)
ENDIF(tvmet_FOUND)