#-----------------------------------------------------------------------------#
# CMakeLists for mesh_common directory
#-----------------------------------------------------------------------------#
target_sources(nut-meshes
  INTERFACE
  ${handy_prefix}/spherical/Spherical_Mesh.h
  ${handy_prefix}/spherical/Spherical_Mesh.i.h
  PRIVATE
  ${handy_prefix}/spherical/Spherical_Mesh.cc
  )

#-----------------------------------------------------------------------------#
# Unit tests.
#-----------------------------------------------------------------------------#
target_sources(nut-meshes-unit-tests
  PRIVATE
  ${handy_prefix}/spherical/test/test_Spherical_Mesh.cc
  )

# End of file
