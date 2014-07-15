bin_PROGRAMS    = nonlinElasQuadTet axi contract angle entropicSpring c0membranegl c0membrane c0_membrane_configurational springshell
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
c0membranegl_SOURCES    = c0membranegl.cc
c0membrane_SOURCES    = c0membrane.cc
c0_membrane_configurational_SOURCES = c0_membrane_configurational.cc
springshell_SOURCES = springshell.cc
LDFLAGS    = -L$(blitz_libraries) \
	-L../                          \
	-L../../VoomMath/                  \
	-L../../Materials/            \
	-L../../Geometry/             \
	-L../../Quadrature/           \
	-L../../Shape/                 
LDADD      = 	\
	-lQuadrature	\
        -lElements	\
	-lMaterials	\
	-lGeometry	\
	-lShape		\
	-lVoomMath	\
	-lblitz		\
	-lm

