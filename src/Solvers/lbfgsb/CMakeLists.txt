project(lbfgsb)

if(VOOM_USE_OPENMP)
	add_library(lbfgsb lbfgsb-routines-parallel.f)
else()
	add_library(lbfgsb lbfgsb-routines.f)
endif()