
#-----------------------------------------------------------------------------#
# Add header files.Note that these will be exported to the parent scope.
# This is what allows cinch to organize this library.
#-----------------------------------------------------------------------------#

target_sources(nut-meshes
  INTERFACE
  ${handy_prefix}/mesh_adaptors/Cartesian_Mesh_Interface.h
  ${handy_prefix}/mesh_adaptors/Spherical_Mesh_Interface.h
 )

#-----------------------------------------------------------------------------#
# Unit tests.
#-----------------------------------------------------------------------------#

target_sources(nut-meshes-unit-tests
  PRIVATE
  ${handy_prefix}/mesh_adaptors/test/test_Cartesian_Mesh_Interface.cc
  ${handy_prefix}/mesh_adaptors/test/test_Spherical_Mesh_Interface.cc
  )

# End of file
