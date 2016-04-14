bin_PROGRAMS    = archaea archaeaS archaeaT4 archaeaThermal archaeaPr minArchaea
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
		-I$(srcdir)/../../Mesh/		\
		-I$(blitz_includes) -I$(tvmet_includes)
archaea_SOURCES = archaea.cc
archaeaS_SOURCES = archaeaS.cc
archaeaT4_SOURCES = archaeaT4.cc
archaeaThermal_SOURCES = archaeaThermal.cc
archaeaPr_SOURCES = archaeaPr.cc
minArchaea_SOURCES = minArchaea.cc
AM_LDFLAGS    = -L$(blitz_libraries) \
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Body/                   \
        -L../../Elements/		\
	-L../../VoomMath/		\
        -L../../Mesh/                   \
	-L/usr/lib/gcc/i486-linux-gnu/4.2/	\
	-L../../Body/	
LDADD      = -lblitz              \
        -lBody				\
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lElements			\
	-lVoomMath			\
        -lMesh                          \
        -lpthread                       \
        -ldl                            \
	-llapack -lblas			\
	-lgfortran
