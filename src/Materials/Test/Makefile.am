bin_PROGRAMS 	= test sceTest gl 
AM_CPPFLAGS	= -I ./                             \
                  -I ../                            \
                  -I ../../                         \
                  -I ../../Geometry/                \
                  -I$(blitz_includes)               \
                  -I$(tvmet_includes)               \
                  -I ../../VoomMath/
AM_LDFLAGS 	= -L$(blitz_libraries)              \
                  -L../                             \
                  -L../../Geometry/                 \
                  -L../../VoomMath/
LDADD		= -lblitz                           \
                  -lMaterials                       \
                  -lGeometry                        \
                  -lVoomMath
test_SOURCES 	= test.cpp 
sceTest_SOURCES = sce-test.cc
gl_SOURCES	= gl.cc

