bin_PROGRAMS    = LMEhomog
INCLUDES        =-I$(srcdir)/../../       		\
		-I$(srcdir)/../../Body/			\
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
LMEhomog_SOURCES    =  LMEhomog.cc 
AM_LDFLAGS    = -L$(blitz_libraries) 	\
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Elements/		\
	-L../../VoomMath/			\
	-L/usr/lib/gcc/i486-linux-gnu/4.2/	\
	-L../../Body/	
LDADD      = -lblitz              	\
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lElements			\
	-lBody				\
	-lVoomMath			\
	-llapack -lblas			\
	-lgfortran
