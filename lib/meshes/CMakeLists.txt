#------------------------------------------------------------------------------#
# Copyright (c) 2020 Triad National Security etc
# All rights reserved.
#------------------------------------------------------------------------------#


string(CONCAT handy_prefix
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/>"
)

message(STATUS "handy handy_prefix:" ${handy_prefix})

# -------------------------------- Unit Testing --------------------------------
option(ENABLE_UNIT_TESTS "Enable unit testing" ON)
option(ENABLE_COLOR_UNIT_TESTS "Enable colorized unit testing output" ON)

if(DEFINED ENABLE_UNIT_TESTS)

  if(NOT TARGET GTest::gtest)

    # set(DEMO_MESHES_TEST_SOURCES CACHE INTERNAL "source files for unit tests" FORCE)
    if (NOT DEFINED GTEST_ROOT)
      if (DEFINED ENV{GTEST_ROOT})
        set(GTEST_ROOT $ENV{GTEST_ROOT})
      else()
         message(WARNING "Warning: environment variable GTEST_DIR is not defined.")
         message(WARNING "nut-meshes-unit-tests will not be built")
         message(WARNING "To make nut_unittests, set GTEST_DIR to top level gtest directory.")
         return()
      endif()
    endif()

    set(GTEST_CMAKE_PATH ${GTEST_ROOT}/lib/cmake)

    message(STATUS "googletest: GTEST_CMAKE_PATH is ${GTEST_CMAKE_PATH}")

    # Define helper functions and macros used by Google Test.
    # include(${GTEST_ROOT}/cmake/internal_utils.cmake)
    find_package(GTest REQUIRED
      HINTS ${GTEST_CMAKE_PATH}
      )

    add_executable(nut-meshes-unit-tests)

  endif()
endif() # ENABLE_UNIT_TESTS


# Now add the nut-meshes library
add_library(nut-meshes)

# ------------------- Now recurse through source directories -------------------
add_subdirectory(base)
add_subdirectory(mesh_common)
add_subdirectory(cartesian)
add_subdirectory(spherical)
add_subdirectory(geometry)
add_subdirectory(mesh_adaptors)

message(STATUS "In lib/meshes, CMAKE_CURRENT_SOURCE_DIR = " ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(nut-meshes
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_features(nut-meshes PRIVATE cxx_std_17)

# TO DO:
# Add install target

# Export target properties

# ---------------------- Now define unit test properties ----------------------
if(DEFINED ENABLE_UNIT_TESTS)

  target_include_directories(nut-meshes-unit-tests
    PUBLIC
      ${GTEST_ROOT}
      ${GTEST_ROOT}/include
    PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

  target_link_libraries(nut-meshes-unit-tests
    nut-meshes
    GTest::gtest
    GTest::gtest_main
    )

  target_compile_features(nut-meshes-unit-tests PRIVATE cxx_std_17)

endif()

# End of CMakeLists.txt
