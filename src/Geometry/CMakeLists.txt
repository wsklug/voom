project(Geometry)

set(SRCS ShellGeometry.cc UnduloidGeometry.cc)

add_library(Geometry ${SRCS})

target_include_directories(Geometry PUBLIC ${PROJECT_SOURCE_DIR} )
target_include_directories(Geometry PUBLIC ${PROJECT_SOURCE_DIR}/.. )
target_include_directories(Geometry PUBLIC ${PROJECT_BINARY_DIR}/..)
target_include_directories(Geometry PUBLIC ${GSL_INCLUDE_DIR})
target_link_libraries(Geometry ${GSL_LIBRARIES})

target_include_directories(Geometry PUBLIC ${tvmet_INCLUDE_DIRS} )
target_include_directories(Geometry PUBLIC ${blitz_INCLUDE_DIRS} )
target_include_directories(Geometry PUBLIC ${EIGEN_INCLUDE_DIRS})
