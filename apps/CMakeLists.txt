# GCC etc
if(NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    include(FindOpenMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        message(STATUS "open mp c flags:" ${OpenMP_C_FLAGS})
        message(STATUS "open mp CXX flags:" ${OpenMP_CXX_FLAGS})
    else()
        message(STATUS "no omp found!")
    endif()

else()

    message(STATUS "Custom configuring OMP for Clang")
    set(GOMP_LIB_PATH "/Users/tkelley/exports/gcc-5.2.0/lib")
    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${GOMP_LIB_PATH})
    link_directories( ${GOMP_LIB_PATH})

    set(OpenMP_C_FLAGS " -fopenmp")
    set(OpenMP_CXX_FLAGS " -fopenmp")
    set(OpenMP_FOUND TRUE)

#     if(NOT DEFINED ENV ${C_INCLUDE_PATH})
#         C_INCLUDE_PATH=/Users/tkelley/wnld/llvm/omp/libomp_oss/build-clang-3.7.0/src
#     endif()
endif() # Clang


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS ${OpenMP_CXX_FLAGS})


# Avoiding setting flags directly
# set(bh-3_COMPILE_FLAGS "  -O3 -fopenmp")

# set(GOMP_LIB_PATH "/Users/tkelley/exports/gcc-5.2.0/lib")
# set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${GOMP_LIB_PATH})
# link_directories( ${GOMP_LIB_PATH})



# Intel
# set(bh-3_COMPILE_FLAGS " -openmp -xcore-avx2 -qopt-report5 -profile-loops=all ")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${bh-3_COMPILE_FLAGS}")

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

