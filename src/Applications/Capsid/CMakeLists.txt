project(Capsid)

if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
	configure_file(${CMAKE_SOURCE_DIR}/cmake/Windows/getopt.c ${PROJECT_BINARY_DIR} COPYONLY)
	configure_file(${CMAKE_SOURCE_DIR}/cmake/Windows/getopt.h ${PROJECT_BINARY_DIR} COPYONLY)
	configure_file(${CMAKE_SOURCE_DIR}/cmake/Windows/unistd.h ${PROJECT_BINARY_DIR} COPYONLY)
endif()

add_library(HelperFunctions OBJECT HelperFunctions.cc)
target_include_directories(HelperFunctions PUBLIC ${VTK_INCLUDE_DIRS}
	${tvmet_INCLUDE_DIRS} ${blitz_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/src	${CMAKE_SOURCE_DIR}/src/Node ${CMAKE_BINARY_DIR}/src
	${PROJECT_BINARY_DIR})
	
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	add_library (getopt OBJECT getopt.c)
endif()

#set(ENV{LIB} $ENV{LIB} "C:/Software/OpenBLAS/lib" PATH)
find_package(BLAS REQUIRED)

#OpenBLAS library provides some lapack functions as well. Otherwise we need to link
# lapack separately
set(LINK_LAPACK TRUE)
if(BLAS_FOUND)
	if(BLAS_openblas_LIBRARY)
		set(LINK_LAPACK FALSE)
	endif(BLAS_openblas_LIBRARY)
endif(BLAS_FOUND)

if(LINK_LAPACK)
	find_package(LAPACK)
endif(LINK_LAPACK)

find_package(Threads)

set(voom_LIBRARIES Body Solvers Materials Elements Quadrature Shape VoomMath Geometry Mesh Model)

set(capsid_targets asphNoRescale annealing asphBrownian asphHarmonic dual relaxMorse 
    asphmodVsFVKcont indentWparticles edgeStrainAndValence harmonic2D Morse2D Morse2DPeriodic)

foreach(curr_target ${capsid_targets})

	if(CMAKE_SYSTEM_NAME STREQUAL "Windows" AND curr_target STREQUAL dual)
		add_executable(${curr_target} $<TARGET_OBJECTS:getopt> ${curr_target}.cc)
	else()
		add_executable(${curr_target} $<TARGET_OBJECTS:HelperFunctions> ${curr_target}.cc)
	endif()
	target_link_libraries(${curr_target} ${blitz_LIBRARIES})
	
	foreach(CURRL ${voom_LIBRARIES})
		target_include_directories(${curr_target} PRIVATE $<TARGET_PROPERTY:${CURRL},INTERFACE_INCLUDE_DIRECTORIES>)
		if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
			target_include_directories(${curr_target} PRIVATE ${PROJECT_BINARY_DIR})
		endif()
		target_link_libraries(${curr_target} ${CURRL})
	endforeach(CURRL)
	
	target_link_libraries(${curr_target} ${VTK_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
	
endforeach(curr_target)

set (helper_targets densityPlot rotatePoints createTfromOutline rmsd cellDataToPointData sphHarmCoeff)

foreach(curr_target ${helper_targets})

	add_executable(${curr_target} ${curr_target}.cc)
	target_include_directories(${curr_target} PRIVATE ${PROJECT_BINARY_DIR} ${tvmet_INCLUDE_DIRS})
	target_link_libraries(${curr_target} ${VTK_LIBRARIES})
	
endforeach(curr_target)


#Using BOOST header files
#set(BOOST_ROOT "C:/Software/intel/boost_1_62_0")
find_package(Boost REQUIRED)
target_include_directories(sphHarmCoeff PRIVATE ${Boost_INCLUDE_DIRS})