bin_PROGRAMS    = test3dEl testMem testMemMod testPotential
AM_CPPFLAGS     =                 	\
	-I$(blitz_includes)            	\
	-I$(tvmet_includes)            	\
	-I$(srcdir)/../          	\
	-I$(srcdir)/../../          	\
	-I$(srcdir)/../../VoomMath/	\
	-I$(srcdir)/../../Node/         \
	-I$(srcdir)/../../Materials/    \
	-I$(srcdir)/../../Geometry/     \
	-I$(srcdir)/../../Quadrature/   \
	-I$(srcdir)/../../Shape/    
test3dEl_SOURCES = Element3D_TEST.cc
testMem_SOURCES = MembraneTest.cc
testMemMod_SOURCES = MembraneModTest.cc
testPotential_SOURCES = TestPotential.cc
LDFLAGS    = -L$(blitz_libraries) \
	-L../                     \
	-L../../VoomMath/         \
	-L../../Materials/        \
	-L../../Geometry/         \
	-L../../Quadrature/       \
	-L../../Shape/                 
LDADD      = 	        \
	-lQuadrature	\
        -lElements	\
	-lMaterials	\
	-lGeometry	\
	-lShape		\
	-lVoomMath	\
	-lblitz		\
	-lm
