bin_PROGRAMS    = app
AM_CPPFLAGS =-I$(srcdir)/../../       		\
		-I$(srcdir)/../../Node/           	\
		-I$(srcdir)/../../Shape/          	\
		-I$(srcdir)/../../Quadrature/          	\
		-I$(srcdir)/../../Materials/       	\
		-I$(srcdir)/../../Elements/       	\
		-I$(srcdir)/../../VoomMath/      	\
		-I$(blitz_includes) -I$(tvmet_includes)
#
app_SOURCES    =  app.cc 
app_LDFLAGS    = -L$(blitz_libraries) \
        -L../                          	\
	-L../../Quadrature/		\
        -L../../Shape/			\
        -L../../Elements/		
app_LDADD      = -lblitz              \
	-lQuadrature                   	\
	-lElements			\
        -lShape				