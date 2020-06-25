#-----------------------------------------------------------------------------#
# CMakeLists for mesh_common directory
#-----------------------------------------------------------------------------#

target_sources(nut-meshes
  INTERFACE
  ${handy_prefix}/mesh_common/Cell.h
  ${handy_prefix}/mesh_common/Cell_Face_Descriptor.h
  ${handy_prefix}/mesh_common/Edge.h
  ${handy_prefix}/mesh_common/Face.h
  ${handy_prefix}/mesh_common/Mesh_Element.h
  ${handy_prefix}/mesh_common/Ray.h
  ${handy_prefix}/mesh_common/Vector.h
  ${handy_prefix}/mesh_common/Vertex.h
  ${handy_prefix}/mesh_common/boundary.h
  )

#-----------------------------------------------------------------------------#
# Unit tests.
#-----------------------------------------------------------------------------#
set(D_B_MESH_COMMON_TEST_SOURCES
  ${handy_prefix}/mesh_common/test/test_Cell.cc
  ${handy_prefix}/mesh_common/test/test_Face.cc
  ${handy_prefix}/mesh_common/test/test_Vector.cc
  ${handy_prefix}/mesh_common/test/test_Edge.cc
  ${handy_prefix}/mesh_common/test/test_Mesh_Element.cc
  ${handy_prefix}/mesh_common/test/test_Vertex.cc)

target_sources(nut-meshes-unit-tests
  PRIVATE
  ${D_B_MESH_COMMON_TEST_SOURCES}
  )

# End of file
