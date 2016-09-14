project(src)

set(SUBD Body Elements Geometry Materials 
	Mesh Model Quadrature Shape Solvers 
	VoomMath Applications) 

IF(tvmet_FOUND)
	set(HAVE_TVMET 1)
ENDIF(tvmet_FOUND)

IF(blitz_FOUND)
	set(HAVE_BLITZ 1)
ENDIF(blitz_FOUND)

MESSAGE(STATUS "Checking if uint needs to be typedef")
include(${CMAKE_SOURCE_DIR}/cmake/Modules/CheckTypedefUint.cmake)
configure_file(${CMAKE_SOURCE_DIR}/cmake/Modules/CheckUint.cc ${CMAKE_BINARY_DIR}/CMakeFiles)
check_typedef_uint()

IF(UINT_TYPEDEF_NEEDED)
	set(USE_UINT_TYPEDEF 1)
	MESSAGE(STATUS "Checking if uint needs to be typedef - YES")
ELSE(UINT_TYPEDEF_NEEDED)
	MESSAGE(STATUS "Checking if uint needs to be typedef - NO")
ENDIF(UINT_TYPEDEF_NEEDED)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

foreach(CURRD IN LISTS SUBD)
	add_subdirectory(${CURRD})
endforeach(CURRD)
