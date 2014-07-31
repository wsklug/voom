bin_PROGRAMS 	= test
AM_CPPFLAGS=-I$(srcdir)/.. -I$(srcdir)/../.. -I$(blitz_includes) -I$(tvmet_includes) 
test_SOURCES 	= test.cpp
test_LDFLAGS 	= -L$(blitz_libraries) -L../ 
test_LDADD	=  -lQuadrature -lblitz

