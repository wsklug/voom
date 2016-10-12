# - Find a tvmet installation
# The following variables are set if blitz is found.  If blitz is not
# found, blitz_FOUND is set to false.
#  blitz_FOUND         - Set to true when blitz is found.
#  blitz_INCLUDE_DIRS  - Include directories for blitz
#  blitz_LIBRARY_DIRS  - Link directories for blitz libraries
#
# The following cache entries must be set by the user to locate blitz:
#  blitz_ROOT  - The expected directory structure is
#				 ${blitz_ROOT}/include/blitz/blitz.h
#				 ${blitz_ROOT}/lib/blitz.lib
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

SET(blitz_DIR_MESSAGE "blitz not found.  Set the blitz_ROOT cmake cache entry to the ${blitz_DIR_DESCRIPTION}")

FIND_PATH(blitz_INCLUDE_DIRS
   NAMES blitz.h
   PATHS ${blitz_ROOT}
   PATH_SUFFIXES "include/blitz/"
   HINTS $ENV{LD_LIBRARY_PATH}
   DOC "The ${blitz_DIR_DESCRIPTION}"
   )
STRING( REGEX REPLACE "/blitz$" "" blitz_INCLUDE_DIRS ${blitz_INCLUDE_DIRS})
STRING( REGEX REPLACE "/include$" "" blitz_DIR ${blitz_INCLUDE_DIRS})
MESSAGE(STATUS "blitz_DIR found : ${blitz_DIR}")
MESSAGE(STATUS "blitz_INCLUDE_DIRS found : ${blitz_INCLUDE_DIRS}")

   IF(blitz_DIR)
	set(blitz_FOUND 1)
	set(blitz_INCLUDE_DIRS "${blitz_DIR}/include")
	
	FIND_LIBRARY(blitz_LIBRARIES 
	NAMES blitz
	PATHS ${blitz_DIR}/lib)
	
	IF(blitz_LIBRARIES)
		MESSAGE(STATUS "blitz libraries found ${blitz_LIBRARIES}")
		set(blitz_LIBRARY_DIRS "${blitz_DIR}/lib")
	ELSEIF(NOT blitz_LIBRARIES)
		MESSAGE(FATAL_ERROR "blitz libraries not found.")
	ENDIF()
	
   ENDIF(blitz_DIR)
   
mark_as_advanced (blitz_DIR)
   
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
