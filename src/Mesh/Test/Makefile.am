bin_PROGRAMS    = test
INCLUDES        =-I ./                 \
	-I $(blitz_includes)            \
	-I $(tvmet_includes)            \
	-I ../                          \
	-I ../../                       
test_SOURCES    = test.cc
test_LDFLAGS    = -L$(blitz_libraries) \
	-L../                          
test_LDADD      = -lblitz -lMesh
