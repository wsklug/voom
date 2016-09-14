project(VoomMath)

set(SRCS VoomMath.cc)

add_library(VoomMath ${SRCS})

target_include_directories(VoomMath PUBLIC ${PROJECT_SOURCE_DIR} )
target_include_directories(VoomMath PUBLIC ${PROJECT_SOURCE_DIR}/.. )
target_include_directories(VoomMath PUBLIC ${PROJECT_BINARY_DIR}/..)

target_include_directories(VoomMath PUBLIC ${tvmet_INCLUDE_DIRS} )
target_include_directories(VoomMath PUBLIC ${blitz_INCLUDE_DIRS} )

