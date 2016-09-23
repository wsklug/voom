bin_PROGRAMS    = tetCapsid_lbfgsb_quad
#CXXFLAGS= -g -ggdb -W -Wall
INCLUDES        =-I ./                 \
        -I$(blitz_includes)            \
        -I$(tvmet_includes)            \
        -I../                          \
        -I../../                       \
        -I../../../                    \
        -I../../../../                 \
        -I../../VoomMath/                  \
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
tetCapsid_lbfgsb_quad_SOURCES    =  tetCapsid-lbfgsb-quad.cc 
tetCapsid_lbfgsb_quad_LDFLAGS    = -L$(blitz_libraries) \
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Elements/		\
	-L../../VoomMath/			
tetCapsid_lbfgsb_quad_LDADD      = -lblitz              \
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lElements			\
	-lVoomMath				\
	-llapack			\
	-lblas	                        \
	-lgfortran			
