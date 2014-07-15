bin_PROGRAMS    = brownian box boxnl relax periodicbox translate shear triangulate trianglestrain
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
        -L../../Body/             	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Elements/		\
	-L../../VoomMath/		\
	-L../../Mesh/			
LDADD=	-lblitz              		\
        -lBody                    	\
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lElements			\
	-lVoomMath				\
	-lMesh				\
	-lvtkIO 			\
	-lvtkGraphics			\
	-lvtkGenericFiltering 		\
	-lvtkFiltering 			\
	-lvtkRendering			\
	-lvtkCommon 			\
	-lvtksys			\
        -lvtkIO                         \
        -lvtkGraphics                   \
        -lvtkGenericFiltering           \
        -lvtkFiltering                  \
        -lvtkRendering                  \
        -lvtkCommon                     \
        -lvtksys                        \
	-llapack			\
	-lblas				\
	-lifport -lifcore
#
brownian_SOURCES    =  brownian.cc 
#
box_SOURCES    =  box.cc
#
periodicbox_SOURCES    =  periodicbox.cc GelOutput.cc
#
boxnl_SOURCES    =  boxnl.cc 
#
relax_SOURCES    =  relax.cc GelOutput.cc
#
translate_SOURCES    =  translate.cc GelOutput.cc
#
shear_SOURCES    =  shear.cc GelOutput.cc
#
triangulate_SOURCES = triangulate.cc

