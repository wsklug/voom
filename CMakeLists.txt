cmake_minimum_required(VERSION 3.0.2)

project(voom C CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#To help CMake find the below packages, user should set the following variables
#when configuring with CMake
#tvmet_DIR e.g. tvmet_DIR = "/usr/local"
#blitz_DIR e.g. blitz_DIR = "/usr/local"
#VTK_DIR e.g. VTK_DIR = "/usr/local/VTK-7.0"
#GSL_ROOT_DIR e.g. GSL_ROOT_DIR =  "/usr/local"

# set(VTK_DIR "/usr/lib/cmake/vtk-7.0" CACHE PATH "Directory containing VTKConfig.cmake" FORCE)
# set(tvmet_DIR "/usr/local" CACHE PATH "Directory containing include/tvmet/tvmet.h" FORCE)
# set(blitz_DIR "/usr/local" CACHE PATH "Directory containing include/blitz/blitz.h" FORCE)
# set(GSL_ROOT_DIR "/usr" CACHE PATH "Directory containing lib/libgsl.a" FORCE)

#enable_language(Fortran)

find_package(tvmet REQUIRED)
find_package(blitz REQUIRED)
find_package(VTK REQUIRED)
find_package(GSL REQUIRED)

set( VOOM_USE_OPENMP FALSE CACHE BOOL "Flag to turn on/off use of OpenMP" )
mark_as_advanced (VOOM_USE_OPENMP)

#If you want to use OpenMP then add "-D VOOM_USE_OPENMP=ON" to the cmake command argument list
if (VOOM_USE_OPENMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
endif()

include(${VTK_USE_FILE})

add_subdirectory(src)

