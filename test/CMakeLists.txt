# CMakeLists.txt for unit tests
if (NOT DEFINED GTEST_ROOT)
  if (DEFINED ENV{GTEST_DIR})
    set(GTEST_ROOT $ENV{GTEST_DIR})
    message(STATUS "GTEST_ROOT is ${GTEST_ROOT}")
  else()
     message(WARNING "Warning: environment variable GTEST_DIR is not defined.")
     message(WARNING "NuT Unittests will not be made")
     message(WARNING "To make nut_unittests, set GTEST_DIR to top level gtest directory.")
     return()
  endif()
endif()

# Define helper functions and macros used by Google Test.
include(${GTEST_ROOT}/cmake/internal_utils.cmake)

config_compiler_and_linker()  # Defined in internal_utils.cmake

include_directories(
	${CMAKE_CURRENT_LIST_DIR}/../lib
  ${RANDOM123_ROOT}/include}
	${GTEST_ROOT}
	${GTEST_ROOT}/include
)

# build gtest
cxx_library(gtest "${cxx_strict}" ${GTEST_ROOT}/src/gtest-all.cc)
cxx_library(gtest_main "${cxx_strict}" ${GTEST_ROOT}/src/gtest_main.cc)
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
target_link_libraries(nut_unittests gtest_main gtest nut ${Murmeln_LIB})

# Build our tests with DBC assertions turned on.
add_definitions(-DREQUIRE_ON)

# Add to CTest
add_test(NAME "nut_unittests" COMMAND "nut_unittests")

# End of file
