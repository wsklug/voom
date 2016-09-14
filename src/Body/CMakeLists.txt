project(Body)

set(SRCS Body.cc GenericBody.cc PotentialBody.cc ViscosityBody.cc ProteinBody.cc)
set(INCDIRS Geometry Node Materials Elements Shape Quadrature VoomMath)

add_library(Body ${SRCS})

#MESSAGE(STATUS ${PROJECT_SOURCE_DIR})

target_include_directories(Body PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(Body PUBLIC ${PROJECT_SOURCE_DIR}/..)
target_include_directories(Body PUBLIC ${PROJECT_BINARY_DIR}/..)



foreach(CURRD ${INCDIRS})
	target_include_directories(Body PUBLIC ${PROJECT_SOURCE_DIR}/../${CURRD} )
endforeach(CURRD)

target_include_directories(Body PUBLIC ${tvmet_INCLUDE_DIRS} )
target_include_directories(Body PUBLIC ${blitz_INCLUDE_DIRS} )
#target_include_directories(Body PUBLIC ${VTK_INCLUDE_DIRS} )