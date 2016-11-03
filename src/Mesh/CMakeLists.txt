project(Mesh)

set(SRCS HalfEdgeMesh.cc)

add_library(Mesh ${SRCS})

target_include_directories(Mesh PUBLIC ${PROJECT_SOURCE_DIR} )
target_include_directories(Mesh PUBLIC ${PROJECT_SOURCE_DIR}/.. )
target_include_directories(Mesh PUBLIC ${PROJECT_BINARY_DIR}/..)

target_include_directories(Mesh PUBLIC ${tvmet_INCLUDE_DIRS} )
target_include_directories(Mesh PUBLIC ${blitz_INCLUDE_DIRS} )
target_include_directories(Mesh PUBLIC ${EIGEN_INCLUDE_DIRS})
