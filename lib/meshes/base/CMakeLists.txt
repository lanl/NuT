#-----------------------------------------------------------------------------#
# CMakeLists for base directory
#-----------------------------------------------------------------------------#

target_sources(nut-meshes
  INTERFACE
  ${handy_prefix}/base/constants.h
  ${handy_prefix}/base/soft_equiv.h
  ${handy_prefix}/base/test_common.h
  )

# set(D_B_BASE_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/)

target_sources(nut-meshes-unit-tests
  PRIVATE
  ${handy_prefix}/base/test/test_base_alive.cc
  )

# End of file
