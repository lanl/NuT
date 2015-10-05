# build google tests for nut

# Borrowed from J. Cohn's test cmake for ENAMR.
# Notes:  CMAKE_SOURCE_DIR is cmake variable set to the
# directory passed as an argument to cmake

# Borrowed from gtest-1.7.0/CMakeLists.txt, which builds gtest
# libraries needed to use gtest as well as example unit tests


# ensure min cmake version, optional - needed when running cmake directly on testG
cmake_minimum_required (VERSION 3.0)


# note:  gtest currently does not have an "install" option
# thus we use an environment variable to find gtest directory

# hold on to that thought!
#if (NOT DEFINED GTEST_ROOT)
#
#	if (DEFINED ENV{GTEST_DIR})
#
#		set(GTEST_ROOT $ENV{GTEST_DIR})
#		message("GTEST_ROOT is ${GTEST_ROOT}")
#
#	else()
#
#	   message("Warning: environment variable GTEST_DIR is not defined.")
#	   message("nut_unittests will not be made")
#	   message("To make nut_unittests, set GTEST_DIR to top level gtest directory.")
#	   return()
#
#	endif()
#
#endif()
#

# Define helper functions and macros used by Google Test.
# include(${GTEST_ROOT}/cmake/internal_utils.cmake)

# config_compiler_and_linker()  # Defined in internal_utils.cmake

# add compile flags for NUT tests if running cmake directly on testG directory
# if (NOT DEFINED NUT_COMPILE_FLAGS)
#
# 	set(NUT_COMPILE_FLAGS "-std=c++11 -g")
# 	# append CMAKE_CXX_FLAGS w/ above variable
# 	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NUT_COMPILE_FLAGS}")
#
# endif()

include_directories(
	${CMAKE_CURRENT_LIST_DIR}/../lib
#    ${RANDOM123_ROOT}/include}
#	${GTEST_ROOT}
#	${GTEST_ROOT}/include
)

# build gtest
# cxx_library(gtest "${cxx_strict}" ${GTEST_ROOT}/src/gtest-all.cc)
# cxx_library(gtest_main "${cxx_strict}" ${GTEST_ROOT}/src/gtest_main.cc)
# target_link_libraries(gtest_main gtest nut)

# build tests
#---------------------------------------------------------
# build nut_unittest1 from sources shown
# TK cargo-culted a glob from MB's approach.
file(GLOB NUT_UNITTESTS
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cc
)

add_executable(nut_unittests ${NUT_UNITTESTS} )

# target_link_libraries(nut_unittests gtest_main)
target_link_libraries(nut_unittests nut)

# Build our tests with DBC assertions turned on.
add_definitions(-DREQUIRE_ON)




