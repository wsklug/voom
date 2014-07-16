bin_PROGRAMS    = nonlinElasQuadTet axi contract angle entropicSpring
AM_CPPFLAGS        =                 	\
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
nonlinElasQuadTet_SOURCES    = nonlinElasQuadTet.cc
contract_SOURCES = contract.cc
axi_SOURCES    = axi.cc
angle_SOURCES    = angle.cc
entropicSpring_SOURCES    = entropicSpring.cc
LDFLAGS    = -L$(blitz_libraries) \
	-L../                          \
	-L../../VoomMath/                  \
	-L../../Materials/            \
	-L../../Geometry/             \
	-L../../Quadrature/           \
	-L../../Shape/                 
LDADD      = 	\
	-lVoomMath	\
	-lMaterials	\
	-lGeometry	\
	-lQuadrature	\
        -lElements	\
	-lShape		\
	-lblitz

