# GCC
set(bh-3_COMPILE_FLAGS "  -O3 -fopenmp")

# set(bh-3_COMPILE_FLAGS " -fopenmp -O3 ")
# set(bh-3_COMPILE_FLAGS "-O3 -msse4.1 -fopenmp -fopt-info-vec-all -march=corei7")

# Intel
# set(bh-3_COMPILE_FLAGS " -openmp -xcore-avx2 -qopt-report5 -profile-loops=all ")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${bh-3_COMPILE_FLAGS}")

file(GLOB bh-3_SOURCES
    cl-args.cc
 )

message( " PROJECT_SOURCE_DIR = " ${PROJECT_SOURCE_DIR} )

add_library(bh-3_support ${bh-3_SOURCES})

include_directories ("${PROJECT_SOURCE_DIR}/lib")
include_directories ("${BOOST_ROOT}/")

add_executable(bh-3 bh-3.cc)

target_link_libraries(bh-3 bh-3_support nut)

# Build with DBC assertions turned on.
# add_definitions(-DREQUIRE_ON)

