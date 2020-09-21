# CMakeLists.txt for unit tests
if (NOT DEFINED GTEST_INSTALL_ROOT)
  if (DEFINED ENV{GTEST_ROOT})
    set(GTEST_INSTALL_ROOT $ENV{GTEST_ROOT})
    message(STATUS "GTEST_INSTALL_ROOT is ${GTEST_INSTALL_ROOT}")
  else()
     message(WARNING "Warning: environment variable GTEST_ROOT is not defined.")
     message(WARNING "NuT Unittests will not be made")
     message(WARNING "To make nut_unittests, set GTEST_ROOT to top level gtest directory.")
     return()
  endif()
endif()

set(GTEST_CMAKE_PATH ${GTEST_INSTALL_ROOT}/lib/cmake)
message(STATUS "GTEST_CMAKE_PATH: " ${GTEST_CMAKE_PATH})

# Define helper functions and macros used by Google Test.
# include(${GTEST_ROOT}/cmake/internal_utils.cmake)
find_package(GTest REQUIRED
  CONFIG
  HINTS ${GTEST_CMAKE_PATH}
  )

# build tests
#---------------------------------------------------------
# build nut_unittest1 from sources shown
file(GLOB NUT_UNITTEST_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cc
)

add_executable(nut_unittests ${NUT_UNITTEST_SRC} )

target_include_directories(nut_unittests
  PUBLIC
    ${GTEST_INSTALL_ROOT}/include
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/../lib
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/../lib/meshes
  )

target_link_libraries(nut_unittests
  GTest::gtest
  GTest::gtest_main
  nut
  )

# Build our tests with DBC assertions turned on.
if(NUT_ENABLE_ASSERTIONS)
  target_compile_definitions(nut_unittests
    -DREQUIRE_ON)
endif()

# Add to CTest
add_test(NAME "nut_unittests" COMMAND "nut_unittests")

# End of file
