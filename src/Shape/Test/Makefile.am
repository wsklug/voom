bin_PROGRAMS    = test testHermite
INCLUDES        = -I$(blitz_includes) -I $(tvmet_includes) \
	-I$(srcdir)/../        	\
	-I$(srcdir)/../../      \
	-I$(srcdir)/../../VoomMath/ \
	-I$(srcdir)/../../Node/ \
	-I$(srcdir)/../../Shape/    
test_SOURCES    = test.cpp
testHermite_SOURCES    = testHermite.cc
LDFLAGS    = -L$(blitz_libraries) \
	-L../                          \
	-L../../VoomMath/                  \
	-L../../Node/                  \
	-L../../Shape/                 
LDADD      = -lShape -lblitz
