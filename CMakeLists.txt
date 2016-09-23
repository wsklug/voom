cmake_minimum_required(VERSION 3.0.2)

project(voom)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#To help CMake find the below packages, user should set the following variables
#when configuring with CMake
#tvmet_DIR e.g. tvmet_DIR = "/usr/local"
#blitz_DIR e.g. blitz_DIR = "/usr/local"
#VTK_DIR e.g. VTK_DIR = "/usr/local/VTK-7.0"
#GSL_ROOT_DIR e.g. GSL_ROOT_DIR =  "/usr/local"

set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose build type, options are: Debug|Release|MinSizeRel|RelWithDebInfo")

# set(VTK_DIR "C:/Software/gcc/VTK-7.0/lib/cmake/vtk-7.0" CACHE PATH "Directory containing VTKConfig.cmake" FORCE)
# set(tvmet_DIR "C:/Software/tvmet" CACHE PATH "Directory containing include/tvmet/tvmet.h" FORCE)
# set(blitz_DIR "C:/Software/gcc/blitz" CACHE PATH "Directory containing include/blitz/blitz.h" FORCE)
# set(GSL_ROOT_DIR "C:/Software/mingw-w64/mingw64/x86_64-w64-mingw32" CACHE PATH "Directory containing lib/libgsl.a" FORCE)

set(VTK_DIR "/usr/lib/cmake/vtk-7.0" CACHE PATH "Directory containing VTKConfig.cmake" FORCE)
set(tvmet_DIR "/usr/local" CACHE PATH "Directory containing include/tvmet/tvmet.h" FORCE)
set(blitz_DIR "/usr/local" CACHE PATH "Directory containing include/blitz/blitz.h" FORCE)
set(GSL_ROOT_DIR "/usr" CACHE PATH "Directory containing lib/libgsl.a" FORCE)

enable_language(Fortran)

find_package(tvmet REQUIRED)
find_package(blitz REQUIRED)
find_package(VTK REQUIRED)
find_package(GSL REQUIRED)

include(${VTK_USE_FILE})

add_subdirectory(src)

