## Makefile.am -- Process this file with automake to produce Makefile.in
AM_CPPFLAGS=	-I$(srcdir)/..			\
		-I$(srcdir)/../Geometry/       	\
		-I$(srcdir)/../Node/           	\
		-I$(srcdir)/../Shape/          	\
		-I$(srcdir)/../Materials/       \
		-I$(srcdir)/../VoomMath/      	\
		-I$(srcdir)/../Quadrature/	\
		-I$(blitz_includes) -I$(tvmet_includes)
lib_LIBRARIES=libElements.a
libElements_a_SOURCES= \
	Element.cc NobleElem.cc LoopGhostBC.cc Dirichlet.cc \
	C0MembraneGL.cc S2.cc MembraneGLImplicitMass.cc AffinityElement.cc \
	SpringShell.cc RigidHemisphereAL.cc RigidPlateAL.cc BrickElement.cc \
	Element3D.cc C0MembraneShear.cc C0MembraneStretch.cc PotentialElement.cc \
	BrownianKick.cc RadialSpring.cc
