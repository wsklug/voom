bin_PROGRAMS    = FiniteContractionApp
INCLUDES        =-I$(srcdir)/../../       		\
		-I$(srcdir)/../../Body/		\
		-I$(srcdir)/../../Model/		\
		-I$(srcdir)/../../Geometry/       	\
		-I$(srcdir)/../../Node/           	\
		-I$(srcdir)/../../Shape/          	\
		-I$(srcdir)/../../Materials/       	\
		-I$(srcdir)/../../Elements/       	\
		-I$(srcdir)/../../VoomMath/      	\
		-I$(srcdir)/../../Quadrature/		\
		-I$(srcdir)/../../Solvers/		\
		-I$(blitz_includes) -I$(tvmet_includes)
FiniteContractionApp_SOURCES    =  FiniteContractionApp.cc 
AM_LDFLAGS    = -L$(blitz_libraries) \
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Elements/		\
	-L../../VoomMath/			\
	-L/usr/lib/gcc/i486-linux-gnu/4.2/		
LDADD      = -lblitz              \
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lElements			\
	-lVoomMath				\
	-llapack -lblas			\
	-lgfortran
