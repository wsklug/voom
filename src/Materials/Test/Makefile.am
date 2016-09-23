bin_PROGRAMS 	= PotentialTest 
AM_CPPFLAGS	= -I ./                             \
                  -I ../                            \
                  -I ../../                         \
                  -I ../../Geometry/                \
                  -I$(blitz_includes)               \
                  -I$(tvmet_includes)               \
                  -I ../../VoomMath/		    \
		  -I ../../Node/
AM_LDFLAGS 	= -L$(blitz_libraries)              \
                  -L../                             \
                  -L../../Geometry/                 \
                  -L../../VoomMath/ 		    \
		  -L../../Node/
LDADD		= -lblitz                           \
                  -lMaterials                       \
                  -lGeometry                        \
                  -lVoomMath
PotentialTest_SOURCES = PotentialTest.cc 

