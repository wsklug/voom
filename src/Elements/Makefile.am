## Makefile.am -- Process this file with automake to produce Makefile.in
AM_CPPFLAGS=	-I$(srcdir)/..			\
		-I$(srcdir)/../Geometry/       	\
		-I$(srcdir)/../Node/           	\
		-I$(srcdir)/../Shape/          	\
		-I$(srcdir)/../Material/       	\
		-I$(srcdir)/../VoomMath/      	\
		-I$(srcdir)/../Quadrature/	\
		-I$(blitz_includes) -I$(tvmet_includes)
lib_LIBRARIES=libElements.a
libElements_a_SOURCES=Element.cc NobleElem.cc LoopGhostBC.cc AffinityElement.cc
