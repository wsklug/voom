project(Capsid)

set(voom_LIBRARIES Body Solvers Materials Elements Quadrature Shape VoomMath Geometry Mesh Model)

add_library(HelperFunctions OBJECT HelperFunctions.cc)
target_include_directories(HelperFunctions PUBLIC ${VTK_INCLUDE_DIRS}
	${tvmet_INCLUDE_DIRS} ${blitz_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/src	${CMAKE_SOURCE_DIR}/src/Node ${CMAKE_BINARY_DIR}/src
	${PROJECT_BINARY_DIR})

foreach(CURRL ${voom_LIBRARIES})
	target_include_directories(HelperFunctions PRIVATE $<TARGET_PROPERTY:${CURRL},INTERFACE_INCLUDE_DIRECTORIES>)
		if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
			target_include_directories(HelperFunctions PRIVATE ${PROJECT_BINARY_DIR})
                endif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
endforeach(CURRL)

if(NOT BLAS_FOUND)
	find_package(BLAS REQUIRED)
endif(NOT BLAS_FOUND)
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

#find_package(Threads)

set(capsid_targets asphericityMorse annealing brownian asphericityHarm dual asphericityCont 
    indentWparticles harmonic2D Morse2D Morse2DPeriodic reOrient twoMeshAnneal indentContinuum
    twoMeshIndent OPSAsphericity OPSBrownianAnneal dual2 testCode)

IF(WIN32 AND BLA_VENDOR MATCHES "Intel[0-9]*_64lp_seq")
	link_directories("${MKL_ROOT}/lib/intel64")
ENDIF()
	
foreach(curr_target ${capsid_targets})
	add_executable(${curr_target} $<TARGET_OBJECTS:HelperFunctions> ${curr_target}.cc)
	target_link_libraries(${curr_target} PRIVATE ${blitz_LIBRARIES})
	
	foreach(CURRL ${voom_LIBRARIES})
		target_include_directories(${curr_target} PRIVATE $<TARGET_PROPERTY:${CURRL},INTERFACE_INCLUDE_DIRECTORIES>)
		if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
			target_include_directories(${curr_target} PRIVATE ${PROJECT_BINARY_DIR})
		endif()
		target_link_libraries(${curr_target} PRIVATE ${CURRL})
	endforeach(CURRL)
	IF(WIN32 AND BLA_VENDOR MATCHES "Intel[0-9]*_64lp_seq")
		target_link_libraries(${curr_target} PRIVATE mkl_intel_ilp64.lib mkl_sequential.lib mkl_core.lib)		
	ELSE()
		target_link_libraries(${curr_target} PUBLIC ${LAPACK_LIBRARIES} PUBLIC ${BLAS_LIBRARIES})
	ENDIF()
	target_link_libraries(${curr_target} PRIVATE ${VTK_LIBRARIES})
	
endforeach(curr_target)

if(NOT OPENMP_FOUND)
    find_package(OpenMP)
endif(NOT OPENMP_FOUND)
if(OPENMP_FOUND)
    target_compile_options(dual2 PRIVATE ${OPENMP_CXX_FLAGS})
    target_link_libraries(dual2 PRIVATE "gomp")
endif(OPENMP_FOUND)

set (helper_targets densityPlot rotatePoints createTfromOutline rmsd sphHarmCoeff)

foreach(curr_target ${helper_targets})

	add_executable(${curr_target} ${curr_target}.cc)
	target_include_directories(${curr_target} PRIVATE ${PROJECT_BINARY_DIR} ${tvmet_INCLUDE_DIRS})
	target_link_libraries(${curr_target} ${VTK_LIBRARIES})
	
endforeach(curr_target)

#Using BOOST header files
find_package(Boost REQUIRED)
target_include_directories(sphHarmCoeff PRIVATE ${Boost_INCLUDE_DIRS})
