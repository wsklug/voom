project(Elements)

set(SRCS Element.cc NobleElem.cc LoopGhostBC.cc Dirichlet.cc
	C0MembraneGL.cc S2.cc MembraneGLImplicitMass.cc AffinityElement.cc
	SpringShell.cc RigidHemisphereAL.cc RigidPlateAL.cc BrickElement.cc
	Element3D.cc C0MembraneShear.cc C0MembraneStretch.cc PotentialElement.cc
	BrownianKick.cc BrownianKickOPS.cc RadialSpring.cc)
	
set(INCDIRS Geometry Node Materials Shape Quadrature VoomMath)

add_library(Elements ${SRCS})

target_include_directories(Elements PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(Elements PUBLIC ${PROJECT_SOURCE_DIR}/..)
target_include_directories(Elements PUBLIC ${PROJECT_BINARY_DIR}/..)


foreach(CURRD ${INCDIRS})
	target_include_directories(Elements PUBLIC ${PROJECT_SOURCE_DIR}/../${CURRD} )
endforeach(CURRD)

target_include_directories(Elements PUBLIC ${tvmet_INCLUDE_DIRS} )
target_include_directories(Elements PUBLIC ${blitz_INCLUDE_DIRS} )
target_include_directories(Elements PUBLIC ${EIGEN_INCLUDE_DIRS})
#target_include_directories(Elements PUBLIC ${VTK_INCLUDE_DIRS} )