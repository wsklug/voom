## Makefile.am -- Process this file with automake to produce Makefile.in
AM_CPPFLAGS=-I$(srcdir)/.. -I$(srcdir)/../Node -I$(srcdir)/../Geometry -I$(srcdir)/../VoomMath -I$(blitz_includes) -I$(tvmet_includes) 
lib_LIBRARIES=libMaterials.a
libMaterials_a_SOURCES= 	\
	SCElastic.cc 		\
	EvansElastic.cc		\
	TwoPhaseElastic.cc	\
	GLElastic.cc            \
	FVK.cc			\
	StVenant.cc 		\
	LinearizedMaterial.cc 	\
	Hookean2D.cc 		\
	Hookean.cc              \
	CompNeoHookean.cc       \
        ModEvansElastic.cc      \
	EvansElastic_Skewed.cc  \
	EvansElastic_SkewedMin.cc	\
	EvansElastic_Stretch.cc \
        HomogMP.cc		\
	Potential.cc 		\
	LennardJones.cc		\
	LennardJonesFT.cc       \
	SpringPotential.cc	\
	SpringPotentialSQ.cc    \
        Morse.cc		\
	ProteinMorse.cc		\
        ProteinLennardJones.cc
