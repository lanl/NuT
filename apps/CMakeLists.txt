file(GLOB bh-3_support_SOURCES
    cl-args.cc
 )

message( " PROJECT_SOURCE_DIR = " ${PROJECT_SOURCE_DIR} )

add_library(bh-3_support ${bh-3_support_SOURCES})

include_directories ("${PROJECT_SOURCE_DIR}/lib")
include_directories ("${BOOST_ROOT}/")

add_executable(bh-3 bh-3.cc)

target_link_libraries(bh-3 bh-3_support nut ${Murmeln_LIB})

add_executable(bh-3-cart bh-3-cartesian-3D.cc)

target_link_libraries(bh-3-cart bh-3_support nut ${Murmeln_LIB})

# Build with DBC assertions turned on.
# add_definitions(-DREQUIRE_ON)

