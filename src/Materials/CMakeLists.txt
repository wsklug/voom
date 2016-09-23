project(Materials)

set(SRCS SCElastic.cc EvansElastic.cc TwoPhaseElastic.cc	
	GLElastic.cc FVK.cc	StVenant.cc	LinearizedMaterial.cc 	
	Hookean2D.cc	Hookean.cc 	CompNeoHookean.cc       
    ModEvansElastic.cc EvansElastic_Skewed.cc  
	EvansElastic_SkewedMin.cc	EvansElastic_Stretch.cc 
    HomogMP.cc	Potential.cc 	LennardJones.cc		
	LennardJonesFT.cc  SpringPotential.cc	SpringPotentialSQ.cc    
    Morse.cc	ProteinMorse.cc ProteinLennardJones.cc)
	
set(INCDIRS Geometry Node VoomMath)

add_library(Materials ${SRCS})

target_include_directories(Materials PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(Materials PUBLIC ${PROJECT_SOURCE_DIR}/.. )
target_include_directories(Materials PUBLIC ${PROJECT_BINARY_DIR}/..)


foreach(CURRD IN LISTS INCDIRS)
	target_include_directories(Materials PUBLIC ${PROJECT_SOURCE_DIR}/../${CURRD} )
endforeach(CURRD)

target_include_directories(Materials PUBLIC ${tvmet_INCLUDE_DIRS} )
target_include_directories(Materials PUBLIC ${blitz_INCLUDE_DIRS} )