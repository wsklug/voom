bin_PROGRAMS    = relax countValence refine disassemble grow dual indentWparticles indent asphVsFVK asphVsFVKcont asphmodVsFVK asphmodVsFVKcont asphVsFVKMorse
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
# capsid_lbfgsb_SOURCES    =  capsid-lbfgsb.cc 
countValence_SOURCES    =  countValence.cc 
#
relax_SOURCES    =  relax.cc 
#
refine_SOURCES    =  refine.cc 
#
dual_SOURCES    =  dual.cc 
#
disassemble_SOURCES    =  disassemble.cc 
#
grow_SOURCES    =  grow.cc C0MembraneGLoutput.cc
#
indentWparticles_SOURCES    =  indentWparticles.cc 
#
indent_SOURCES = indent.cc
#
asphVsFVK_SOURCES = asphVsFVK.cc
#
asphVsFVKcont_SOURCES = asphVsFVKcont.cc
#
asphmodVsFVK_SOURCES = asphmodVsFVK.cc
#
asphmodVsFVKcont_SOURCES = asphmodVsFVKcont.cc
#
asphVsFVKMorse_SOURCES = asphVsFVKMorse.cc