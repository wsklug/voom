bin_PROGRAMS    = coneproblem
AM_CPPFLAGS =-I$(srcdir)/../../       		\
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
#
AM_LDFLAGS = -L$(blitz_libraries) \
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Body/		\
        -L../../Elements/		\
	-L../../VoomMath/		\
	-L../../Mesh/
#
LDADD	= -lblitz              		\
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
        -lShape				\
	-lBody			\
	-lElements			\
        -lGeometry                      \
	-lVoomMath				\
	-lMesh				\
	-lvtkverdict			\
	-lvtktiff -lvtkpng -lvtkjpeg	\
	-lvtkNetCDF -lvtkNetCDF_cxx	\
	-lvtkexpat			\
	-lvtksqlite			\
	-lvtkzlib			\
	-lvtkmetaio			\
	-lvtkDICOMParser		\
	-lvtkIO 			\
	-lvtkGraphics			\
	-lvtkGenericFiltering 		\
	-lvtkFiltering 			\
	-lvtkCommon 			\
	-lvtksys			\
	-ldl				\
	-llapack			\
	-lblas				\
	-lpthread			\
	-lgfortran
#
coneproblem_SOURCES = coneproblem.cc
