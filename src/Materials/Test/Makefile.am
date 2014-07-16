bin_PROGRAMS 	= test sceTest
INCLUDES	= -I ./                             \
                  -I ../                            \
                  -I ../../                         \
                  -I ../../Geometry/                \
                  -I$(blitz_includes)               \
                  -I$(tvmet_includes)               \
                  -I ../../Math/
test_SOURCES 	= test.cpp 
test_LDFLAGS 	= -L$(blitz_libraries)              \
                  -L../                             \
                  -L../../Geometry/                 \
                  -L../../Math/
test_LDADD	= -lblitz                           \
                  -lMaterials                       \
                  -lGeometry                        \
                  -lMath
sceTest_SOURCES = sce-test.cc
sceTest_LDFLAGS = -L$(blitz_libraries)              \
                  -L../                             \
                  -L../../Geometry/                 \
                  -L../../Math/
sceTest_LDADD	= -lblitz                           \
                  -lMaterials                       \
                  -lGeometry                        \
                  -lMath

