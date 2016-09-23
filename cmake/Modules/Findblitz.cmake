# - Find a tvmet installation
# The following variables are set if blitz is found.  If blitz is not
# found, blitz_FOUND is set to false.
#  blitz_FOUND         - Set to true when blitz is found.
#  blitz_INCLUDE_DIRS  - Include directories for blitz
#  blitz_LIBRARY_DIRS  - Link directories for blitz libraries
#
# The following cache entries must be set by the user to locate blitz:
#  blitz_DIR  - This is the directory that contains sub-folder include 
#				which should contain another subfolder called blitz which
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
SET(blitz_FOUND 0)

# Construct consitent error messages for use below.
SET(blitz_DIR_DESCRIPTION "directory that contains sub-folder include which should contain another subfolder called blitz which contains the headers")

SET(blitz_DIR_MESSAGE "blitz not found.  Set the blitz_DIR cmake cache entry to the ${blitz_DIR_DESCRIPTION}")

# Look for Useblitz.cmake in build trees or under <prefix>/include/blitz.
FIND_PATH(blitz_DIR
   NAMES blitz.h
   PATH_SUFFIXES "include/blitz/"
   HINTS $ENV{LIB}
   DOC "The ${blitz_DIR_DESCRIPTION}"
   )

   IF(blitz_DIR)
   
	MESSAGE(STATUS "blitz headers located under ${blitz_DIR}")
	set(blitz_FOUND 1)
	set(blitz_INCLUDE_DIRS "${blitz_DIR}/include")
	set(blitz_LIBRARY_DIRS "${blitz_DIR}/lib")
	
	FIND_LIBRARY(blitz_LIBRARIES 
	NAMES blitz
	PATHS ${blitz_DIR}/lib)
	
	IF(blitz_LIBRARIES)
		MESSAGE(STATUS "blitz libraries found ${blitz_LIBRARIES}")
	ELSEIF(blitz_LIBRARIES)
		MESSAGE(FATAL_ERROR "blitz libraries could not found.")
	ENDIF(blitz_LIBRARIES)
	
   ENDIF(blitz_DIR)
   
   
#-----------------------------------------------------------------------------
IF(blitz_FOUND)
ELSE(blitz_FOUND)
  # blitz not found, explain to the user how to specify its location.
  IF(blitz_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${blitz_DIR_MESSAGE})
  ELSE(blitz_FIND_REQUIRED)
    IF(NOT blitz_FIND_QUIETLY)
      MESSAGE(STATUS ${blitz_DIR_MESSAGE})
    ENDIF(NOT blitz_FIND_QUIETLY)
  ENDIF(blitz_FIND_REQUIRED)
ENDIF(blitz_FOUND)