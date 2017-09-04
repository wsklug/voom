cmake_minimum_required(VERSION 3.0.2)

project(voom)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#To help CMake find the below packages, user should set the following variables
#when configuring with CMake
#tvmet_DIR e.g. tvmet_DIR = "/usr/local"
#blitz_DIR e.g. blitz_DIR = "/usr/local"
#VTK_DIR e.g. VTK_DIR = "/usr/local/VTK-7.0"
#GSL_ROOT_DIR e.g. GSL_ROOT_DIR =  "/usr/local"

#set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose build type, options are: Debug|Release|MinSizeRel|RelWithDebInfo" FORCE)

# set(VTK_DIR "/u/local/apps/vtk/5.8.0/lib/vtk-5.8" CACHE PATH "Directory containing VTKConfig.cmake" FORCE)
# set(tvmet_DIR "/u/local/apps/tvmet/current" CACHE PATH "Directory containing include/tvmet/tvmet.h" FORCE)
# set(blitz_DIR "/u/local/apps/blitz/current" CACHE PATH "Directory containing include/blitz/blitz.h" FORCE)
# set(GSL_ROOT_DIR "/usr" CACHE PATH "Directory containing lib/libgsl.a" FORCE)

cmake_policy(SET CMP0053 NEW)
cmake_policy(SET CMP0054 NEW)

enable_language(Fortran)

find_package(tvmet REQUIRED)
find_package(blitz REQUIRED)
find_package(VTK REQUIRED)
find_package(GSL REQUIRED)
find_package(Eigen REQUIRED)
IF(USE_MULTIPRECISION)
	find_package(MPFR REQUIRED)
	find_package(GMP REQUIRED)
ENDIF(USE_MULTIPRECISION)

include(${VTK_USE_FILE})

if(WIN32)
	include_directories(${CMAKE_SOURCE_DIR}/cmake/Windows)
	add_definitions(-DNOMINMAX)
endif(WIN32)

add_subdirectory(src)

