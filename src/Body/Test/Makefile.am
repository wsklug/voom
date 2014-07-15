bin_PROGRAMS    = test axi
INCLUDES        =-I$(srcdir)/..				\
		-I$(srcdir)/../../       		\
		-I$(srcdir)/../../Geometry/       	\
		-I$(srcdir)/../../Node/           	\
		-I$(srcdir)/../../Shape/          	\
		-I$(srcdir)/../../Materials/       	\
		-I$(srcdir)/../../Elements/       	\
		-I$(srcdir)/../../VoomMath/      	\
		-I$(srcdir)/../../Quadrature/		\
		-I$(blitz_includes) -I$(tvmet_includes)
test_SOURCES    = testCapsid.cpp
axi_SOURCES    = axi.cc
LDFLAGS    = -L$(blitz_libraries) \
        -L../                          \
        -L../../VoomMath/                  \
        -L../../Node/                  \
        -L../../Materials/            \
        -L../../Geometry/             \
	-L../../Elements/		\
        -L../../Quadrature/           \
        -L../../Shape/			\
	-L../../Mesh/
LDADD      = -lblitz              	\
        -lVoomMath                         	\
        -lMaterials                    	\
        -lGeometry                     	\
        -lQuadrature                   	\
	-lElements			\
        -lShape				\
	-lMesh
