bin_PROGRAMS    = testLME testLMEtet
INCLUDES        = -I$(blitz_includes) -I $(tvmet_includes) \
	-I$(srcdir)/../     		\
	-I$(srcdir)/../../      	\
	-I$(srcdir)/../../VoomMath/ 	\
	-I$(srcdir)/../../Node/ 	\
	-I$(srcdir)/../../Shape/    
testLME_SOURCES    = testLME.cpp
testLMEtet_SOURCES = testLMEtet.cpp
LDFLAGS    = -L$(blitz_libraries) 	\
	-L../                          	\
	-L../../VoomMath/               \
	-L../../Node/                  	\
	-L../../Shape/                 
LDADD      = -lShape		\
	     -lblitz 		\
	     -ldl               \
	     -llapack -lblas	\
	     -lgfortran