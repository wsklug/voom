bin_PROGRAMS 	= test
INCLUDES	=-I ./ -I ../ -I ../../Math/   -I$(blitz_includes) -I$(tvmet_includes) 
test_SOURCES 	= testlib.cpp
test_LDFLAGS 	= -L$(blitz_libraries) -L../ -L../../Math/
test_LDADD	= -lblitz -lFEMMath

