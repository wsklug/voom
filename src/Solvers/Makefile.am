## Makefile.am -- Process this file with automake to produce Makefile.in
AM_CPPFLAGS=-I$(srcdir)			\
	-I$(srcdir)/..			\
	-I$(srcdir)/../Body/            \
	-I$(srcdir)/../Model/           \
        -I$(srcdir)/../VoomMath/	\
        -I$(srcdir)/../Node/            \
        -I$(srcdir)/../Materials/       \
        -I$(srcdir)/../Geometry/        \
        -I$(srcdir)/../Quadrature/      \
	-I$(srcdir)/../Elements/        \
	-I$(srcdir)/../Body/            \
        -I$(srcdir)/../Shape/		\
	-I$(blitz_includes)		\
	-I$(tvmet_includes)
lib_LIBRARIES=libSolvers.a
libSolvers_a_SOURCES=		\
	DirectLinearSolver.cc 	\
	ConjugateGradientWSK.cc \
	ViscousRelaxation.cc 	\
	Lbfgsb.cc		\
	lbfgsb-routines.f	\
	CGfast.cc 		\
	SimulatedAnnealing.cc	\
	MontecarloTwoStages.cc	\
	CgDescent.cc		\
	cg_descent.c		\
	BrownianDynamics.cc	\
	BrownianDynamics3D.cc	\
	MontecarloProtein.cc
