## Makefile.am -- Process this file with automake to produce Makefile.in
AM_CPPFLAGS =	-I$(srcdir)/..			\
		-I$(srcdir)/../Geometry/       	\
		-I$(srcdir)/../Node/           	\
		-I$(srcdir)/../Shape/          	\
		-I$(srcdir)/../Materials/     	\
		-I$(srcdir)/../Elements/       	\
		-I$(srcdir)/../VoomMath/      	\
		-I$(srcdir)/../Quadrature/	\
		-I$(blitz_includes) 		\
		-I$(tvmet_includes)		\
		-I$(vtk_includes)		

lib_LIBRARIES = libBody.a
libBody_a_SOURCES = Body.cc GenericBody.cc PotentialBody.cc ViscosityBody.cc ProteinBody.cc
