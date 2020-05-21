# CMakeLists.txt for unit tests
if (NOT DEFINED GTEST_ROOT)
  if (DEFINED ENV{GTEST_ROOT})
    set(GTEST_ROOT $ENV{GTEST_ROOT})
    message(STATUS "GTEST_ROOT is ${GTEST_ROOT}")
  else()
     message(WARNING "Warning: environment variable GTEST_ROOT is not defined.")
     message(WARNING "NuT Unittests will not be made")
     message(WARNING "To make nut_unittests, set GTEST_ROOT to top level gtest directory.")
     return()
  endif()
endif()

set(GTEST_CMAKE_PATH ${GTEST_ROOT}/cmake/cmake)

# Define helper functions and macros used by Google Test.
# include(${GTEST_ROOT}/cmake/internal_utils.cmake)
find_package(GTest REQUIRED
  HINT ${GTEST_CMAKE_PATH})

# config_compiler_and_linker()  # Defined in internal_utils.cmake

include_directories(
	${CMAKE_CURRENT_LIST_DIR}/../lib
  ${RANDOM123_ROOT}/include}
	${GTEST_ROOT}
	${GTEST_ROOT}/include
)

# build tests
#---------------------------------------------------------
# build nut_unittest1 from sources shown
file(GLOB NUT_UNITTESTS
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cc
)

add_executable(nut_unittests ${NUT_UNITTESTS} )

target_include_directories(nut_unittests
  PUBLIC
    ${GTEST_ROOT}/include
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

target_link_libraries(nut_unittests
  GTest::gtest
  GTest::gtest_main
  nut
  #${Murmeln_LIB}
  )

# Build our tests with DBC assertions turned on.
if(NUT_ENABLE_ASSERTIONS)
  target_compile_definitions(nut_unittests
    -DREQUIRE_ON)
endif()

# Add to CTest
add_test(NAME "nut_unittests" COMMAND "nut_unittests")

# End of file
