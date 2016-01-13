bin_PROGRAMS    = brownian box boxnl relax periodicbox translate shear affinetest nematictest bareBoneGel gelpropertytester rodliquid inputdriver
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
	-L../../Mesh/			\
	-L../../Elements/
#
LDADD	= -lblitz              		\
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lElements			\
	-lVoomMath				\
	-lMesh				\
	-llapack			\
	-lblas				\
	-lgfortran			\
	-lpthread			\
	-ldl
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
affinetest_SOURCES    =  affinetest.cc GelOutput.cc
#
nematictest_SOURCES    =  nematictest.cc SemiflexibleInput.cc GelOutput.cc
#
bareBoneGel_SOURCES    =  bareBoneGel.cc SemiflexibleInput.cc  GelOutput_2.cc
#
gelpropertytester_SOURCES    =  gelpropertytester.cc GelOutput.cc
#
rodliquid_SOURCES    =  rodliquid.cc GelOutput.cc
#
inputdriver_SOURCES    =  inputdriver.cc SemiflexibleInput.cc GelOutput.cc
