project(Quadrature)

set(SRCS TriangleQuadrature.cc QuadQuadrature.cc TetQuadrature.cc HexQuadrature.cc LineQuadrature.cc BrickQuadrature.cc)

add_library(Quadrature ${SRCS})

target_include_directories(Quadrature PUBLIC ${PROJECT_SOURCE_DIR} )
target_include_directories(Quadrature PUBLIC ${PROJECT_SOURCE_DIR}/.. )
target_include_directories(Quadrature PUBLIC ${PROJECT_BINARY_DIR}/..)

target_include_directories(Quadrature PUBLIC ${tvmet_INCLUDE_DIRS} )
target_include_directories(Quadrature PUBLIC ${blitz_INCLUDE_DIRS} )