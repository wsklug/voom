bin_PROGRAMS    = tether
#CXXFLAGS= -g -ggdb -W -Wall
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
	-I ../../Constraints/          \
	-I ../../Potentials/           \
	-I ../../Body/                 \
	-I ../../Model/                \
	-I ../../Solvers/              \
        -I ../../Shape/
#
tether_SOURCES    =  tether.cc 
tether_LDFLAGS    = -L$(blitz_libraries) \
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Elements/		\
	-L../../Math/
tether_LDADD      = -lblitz              \
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lElements			\
	-lMath				\
	-llapack			\
	-lblas				\
	-lmpe

