project(Model)

set(SRCS Model.cc)
set(INCDIRS Geometry Node Materials Elements Shape Quadrature VoomMath Body Solvers)

add_library(Model ${SRCS})

target_include_directories(Model PUBLIC ${PROJECT_SOURCE_DIR} )
target_include_directories(Model PUBLIC ${PROJECT_SOURCE_DIR}/.. )
target_include_directories(Model PUBLIC ${PROJECT_BINARY_DIR}/..)


foreach(CURRD IN LISTS INCDIRS)
	target_include_directories(Model PUBLIC ${PROJECT_SOURCE_DIR}/../${CURRD} )
endforeach(CURRD)

target_include_directories(Model PUBLIC ${tvmet_INCLUDE_DIRS} )
target_include_directories(Model PUBLIC ${blitz_INCLUDE_DIRS} )
target_include_directories(Model PUBLIC ${EIGEN_INCLUDE_DIRS})