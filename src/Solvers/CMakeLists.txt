project(Solvers)

enable_language(Fortran)

set(SRCS DirectLinearSolver.cc 	
	ConjugateGradientWSK.cc 
	ViscousRelaxation.cc 	
	Lbfgsb.cc		
	lbfgsb-routines.f	
	CGfast.cc 		
	SimulatedAnnealing.cc	
	MontecarloTwoStages.cc	
	CgDescent.cc		
	cg_descent.c		
	BrownianDynamics.cc	
	BrownianDynamics3D.cc	
	MontecarloProtein.cc    
        KMCprotein.cc)
set(INCDIRS Geometry Node Materials Model Elements Shape Quadrature VoomMath Body Elements)

add_library(Solvers ${SRCS})

target_include_directories(Solvers PUBLIC ${PROJECT_SOURCE_DIR} )
target_include_directories(Solvers PUBLIC ${PROJECT_SOURCE_DIR}/.. )
target_include_directories(Solvers PUBLIC ${PROJECT_BINARY_DIR}/..)


foreach(CURRD IN LISTS INCDIRS)
	target_include_directories(Solvers PUBLIC ${PROJECT_SOURCE_DIR}/../${CURRD} )
endforeach(CURRD)

target_include_directories(Solvers PUBLIC ${tvmet_INCLUDE_DIRS} )
target_include_directories(Solvers PUBLIC ${blitz_INCLUDE_DIRS} )