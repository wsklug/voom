bin_PROGRAMS    = minimize minimize_stretch minstretch_series minStretchSeriesSymm minStretchSeriesNONSymm MontecarloHexagonal MontecarloTvirus min_str_LargeT Frustration_LargeT
INCLUDES        =-I$(srcdir)/../../       		\
		-I$(srcdir)/../../Body/		\
		-I$(srcdir)/../../Model/		\
		-I$(srcdir)/../../Geometry/       	\
		-I$(srcdir)/../../Node/           	\
		-I$(srcdir)/../../Shape/          	\
		-I$(srcdir)/../../Materials/       	\
		-I$(srcdir)/../../Elements/       	\
		-I$(srcdir)/../../VoomMath/      	\
		-I$(srcdir)/../../Quadrature/		\
		-I$(srcdir)/../../Solvers/		\
		-I$(srcdir)/../../Mesh/		\
		-I$(blitz_includes) -I$(tvmet_includes)
minimize_SOURCES = minimize.cc
minimize_stretch_SOURCES = minimize_stretch.cc
minstretch_series_SOURCES = minstretch_series.cc
minStretchSeriesSymm_SOURCES = minStretchSeriesSymm.cc
minStretchSeriesNONSymm_SOURCES = minstretch_NonSymmSeries.cc
MontecarloHexagonal_SOURCES = MontecarloHexagonal.cc
MontecarloTvirus_SOURCES = MontecarloTvirus.cc
min_str_LargeT_SOURCES = min_str_LargeT.cc
Frustration_LargeT_SOURCES = Frustration_LargeT.cc
AM_LDFLAGS    = -L$(blitz_libraries) \
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Body/                   \
        -L../../Elements/		\
	-L../../VoomMath/		\
        -L../../Mesh/                   \
	-L/usr/lib/gcc/i486-linux-gnu/4.2/	\
	-L../../Body/	
LDADD      = -lblitz              \
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lElements			\
	-lBody				\
	-lVoomMath			\
        -lMesh                          \
        -lpthread                       \
        -ldl                            \
	-llapack -lblas			\
	-lgfortran
