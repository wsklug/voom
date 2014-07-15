bin_PROGRAMS    = test
CXXFLAGS= -g -ggdb -W -Wall
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
        -I ../../Quadrature/           \
	-I ../../Elements/             \
	-I ../../Body/                 \
	-I ../../Model/                \
	-I ../../Solvers/              \
        -I ../../Shape/
test_SOURCES    = cgDescent.cc
test_LDFLAGS    = -L$(blitz_libraries) \
        -L../                          \
        -L../../Math/                  \
        -L../../Node/                  \
        -L../../Materials/             \
        -L../../Geometry/              \
	-L../../Solvers/               \
	-L../../Model/                 \
        -L../../Shape/
test_LDADD      = -lblitz              \
        -lMath                         \
        -lMaterials                    \
        -lGeometry                     \
	-lModel                        \
	-lSolvers                      \
        -lShape



