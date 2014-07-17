bin_PROGRAMS    = axi domain shock indent
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
axi_SOURCES    = axi.cc
domain_SOURCES    = domain.cc
shock_SOURCES	= shock.cc
indent_SOURCES	= indent.cc
LDFLAGS    = -L$(blitz_libraries) \
        -L../                          \
        -L../../VoomMath/                  \
        -L../../Materials/            \
        -L../../Geometry/             \
	-L../../Model/		\
	-L../../Body/			\
	-L../../Elements/		\
        -L../../Quadrature/           \
        -L../../Shape/			\
        -L../../Solvers/			\
	-L../../Mesh/			\
	-L/sw/lib/gcc4.2/lib/
LDADD      = -lblitz              	\
        -lVoomMath  	\
	-lModel		\
	-lBody		\
        -lMaterials 	\
        -lGeometry  	\
        -lQuadrature	\
	-lElements  	\
        -lShape		\
        -lSolvers	\
	-llapack	\
	-lblas		\
	-lgfortran
