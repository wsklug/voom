project(Shape)

set(SRCS LoopShellShape.cc ShapeTri3.cc ShapeTri6.cc
	ShapeQ4.cc ShapeTet4CP.cc ShapeTet10.cc ShapeHex8.cc
	Hermite.cc MFShape.cc ShapeBrick6.cc ShapeBrick9.cc
	LMEshape.cc LMEshape2.cc LMEtet.cc)
	
set(INCDIRS Node VoomMath)

add_library(Shape ${SRCS})

target_include_directories(Shape PUBLIC ${PROJECT_SOURCE_DIR} )
target_include_directories(Shape PUBLIC ${PROJECT_SOURCE_DIR}/.. )
target_include_directories(Shape PUBLIC ${PROJECT_BINARY_DIR}/..)


foreach(CURRD IN LISTS INCDIRS)
	target_include_directories(Shape PUBLIC ${PROJECT_SOURCE_DIR}/../${CURRD} )
endforeach(CURRD)

target_include_directories(Shape PUBLIC ${tvmet_INCLUDE_DIRS} )
target_include_directories(Shape PUBLIC ${blitz_INCLUDE_DIRS} )