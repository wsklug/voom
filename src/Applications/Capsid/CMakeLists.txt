project(Capsid)

add_library(HelperFunctions OBJECT HelperFunctions.cc)
target_include_directories(HelperFunctions PUBLIC ${VTK_INCLUDE_DIRS}
	${tvmet_INCLUDE_DIRS} ${blitz_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/src	${CMAKE_SOURCE_DIR}/src/Node ${CMAKE_BINARY_DIR}/src)

set(ENV{LIB} $ENV{LIB} "C:/Software/OpenBLAS/lib" PATH)
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

set(capsid_targets asphNoRescale annealing asphBrownian asphHarmonic dual relaxMorse asphmodVsFVKcont indentWparticles)

foreach(curr_target ${capsid_targets})

	add_executable(${curr_target} $<TARGET_OBJECTS:HelperFunctions> ${curr_target}.cc)
	target_link_libraries(${curr_target} ${blitz_LIBRARIES})
	
	foreach(CURRL ${voom_LIBRARIES})
		target_include_directories(${curr_target} PRIVATE $<TARGET_PROPERTY:${CURRL},INTERFACE_INCLUDE_DIRECTORIES>)
		target_link_libraries(${curr_target} ${CURRL})
	endforeach(CURRL)
	
	# if(VTK_BUILD_SHARED_LIBS STREQUAL ON)
		# if(MINGW)
			
			# foreach(vtk_lib ${VTK_LIBRARIES})
				# target_link_libraries(${curr_target} "${vtk_lib}")
			# endforeach(vtk_lib)
			
		# endif(MINGW)
	# else(VTK_BUILD_SHARED_LIBS STREQUAL ON)
		# target_link_libraries(${curr_target} ${VTK_LIBRARIES})
	# endif(VTK_BUILD_SHARED_LIBS STREQUAL ON)
	
	target_link_libraries(${curr_target} ${VTK_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
	
endforeach(curr_target)
