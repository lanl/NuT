#-----------------------------------------------------------------------------#
# CMakeLists for mesh_common directory
#-----------------------------------------------------------------------------#
target_sources(nut-meshes
  INTERFACE
  ${handy_prefix}/cartesian/Cartesian_Mesh.h
  ${handy_prefix}/cartesian/Cartesian_Mesh.i.h
  ${handy_prefix}/cartesian/detail/Cartesian_Mesh_detail.h
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/Cartesian_Mesh.cc
  )

set(CARTESIAN_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/Cartesian_Mesh.cc
  )

#-----------------------------------------------------------------------------#
# Unit tests.
#-----------------------------------------------------------------------------#
target_sources(nut-meshes-unit-tests
  PRIVATE
  ${handy_prefix}/cartesian/test/test_Cartesian_Mesh.cc
  )

# End of file
