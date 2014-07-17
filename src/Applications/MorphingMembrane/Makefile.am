bin_PROGRAMS    = membrane meshMaker
#CXXFLAGS= -g -ggdb -Wall -w
INCLUDES        =-I ./                 \
        -I$(blitz_includes)            \
        -I$(tvmet_includes)            \
        -I../                          \
        -I../../                       \
        -I../../../                    \
        -I../../../../                 \
        -I../../Math/                  \
        -I../../Node/                  \
        -I ../../Materials/            \
        -I ../../Geometry/             \
	-I/home/software/              \
	-I../../AdsExt                 \
        -I ../../Quadrature/           \
	-I ../../Elements/             \
	-I ../../Potentials/           \
	-I ../../Constraints/          \
	-I ../../Body/                 \
	-I ../../Model/			\
	-I ../../Solvers/			\
        -I ../../Shape/
membrane_SOURCES    = membrane.cc
membrane_LDFLAGS    = -L$(blitz_libraries) \
        -L../                          \
        -L../../Math/                  \
        -L../../Node/                  \
        -L../../Materials/            \
        -L../../Geometry/             \
	-L../../Elements/		\
        -L../../Quadrature/           \
	-L../../Model/			\
	-L../../Solvers/		\
        -L../../Shape/			\
	-L/home/software/ifc-9.0/lib/
membrane_LDADD      = -lblitz              	\
        -lMath                         	\
        -lMaterials                    	\
        -lGeometry                     	\
        -lQuadrature                   	\
	-lElements			\
	-lModel				\
	-lSolvers			\
        -lShape				\
	-llapack -lblas			\
	-lifport -lifcore
meshMaker_SOURCES    =  mesh.c brep.c refine.c delaunay.f meshMaker.c
meshMaker_LDFLAGS    = 	-L/home/software/ifc-9.0/lib/  -nofor_main
meshMaker_LDADD      = 	-lifport -lifcore
