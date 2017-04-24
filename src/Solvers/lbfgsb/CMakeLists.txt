project(lbfgsb)
enable_language(Fortran)

IF(WIN32 AND MSVC)
	IF(NOT ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel"))
		MESSAGE("Please set MINGW_GFORTRAN to mingw-gfortran.exe location.")
		set (MINGW_GFORTRAN "C:/Software/mingw-w64/mingw64/bin/gfortran.exe" CACHE FILEPATH "Please provide full path to mingw-gfortran.exe")
	ENDIF()
	#IF(BLA_VENDOR MATCHES "Intel[0-9]*_64lp_seq")
	#	link_directories("${MKL_ROOT}/lib/intel64")
	#ENDIF()
ENDIF()
IF(VOOM_USE_OPENMP)	
	add_library(lbfgsb STATIC lbfgsb-routines-parallel.f)
ELSE()
	add_library(lbfgsb STATIC lbfgsb-routines.f)
ENDIF()

# #set(ENV{LIB} $ENV{LIB} "C:/Software/OpenBLAS/lib" PATH)
# find_package(BLAS REQUIRED)

# #OpenBLAS library provides some lapack functions as well. Otherwise we need to link
# # lapack separately
# set(LINK_LAPACK TRUE)
# if(BLAS_FOUND)
	# if(BLAS_openblas_LIBRARY)
		# set(LINK_LAPACK FALSE)
	# endif(BLAS_openblas_LIBRARY)
# endif(BLAS_FOUND)

# if(LINK_LAPACK)
	# find_package(LAPACK)
# endif(LINK_LAPACK)

# find_package(Threads)

# IF(WIN32)
	# IF(BLA_VENDOR MATCHES "Intel[0-9]*_64lp_seq")
		# MESSAGE("Please ensure that CMake Cache variable MKL_ROOT is set correctly.")
		# MESSAGE("Using MKL_ROOT = ${MKL_ROOT}")		
		# target_include_directories(lbfgsb PRIVATE "${MKL_ROOT}/include")
		# target_link_libraries(lbfgsb PUBLIC mkl_blas95_ilp64.lib mkl_lapack95_ilp64.lib mkl_intel_ilp64.lib mkl_sequential.lib mkl_core.lib)
	# ENDIF()
# ENDIF()