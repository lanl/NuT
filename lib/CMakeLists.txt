# CMakeLists for nut lib.
# Want to accomplish two things:
#   1. build the .cc files into a lib,
#   2. move (and configure if necessary) the header files to the include
#      directory.


set(NUT_HEADERS
  Assert.hh
  Boundary_Cond.hh
  Cartesian_Mesh.h
  Cartesian_Mesh.i.h
  Census.hh
  Density.hh
  Fates.hh
  Log.hh
  Luminosity.hh
  MatState.hh
  Mesh.hh
  Mesh_3D_Cartesian.hh
  Opacity.hh
  Particle.hh
  Planck.hh
  RNG.hh
  Tally.hh
  Temperature.hh
  Vec3D.hh
  Velocity.hh
  apply_event.hh
  constants.hh
  copyright.hh
  decision.hh
  events.hh
  fileio.hh
  lorentz.hh
  mpi_helpers.hh
  partition.hh
  serialize.hh
  sigmas_HBFC.hh
  sim_init.hh
  soft_equiv.hh
  sourcery.hh
  transport.hh
  types.hh
  utilities.hh
  utilities_io.hh
  detail/Cartesian_Mesh_detail.h
  detail/Cell.h
  detail/Cell_Face_Descriptor.h
  detail/Face.h
  detail/Mesh_Element.h
  detail/Ray.h
  detail/Vector.h
  detail/Vertex.h
  detail/boundary.h
)


set(NUT_SRC
  Cartesian_Mesh.cc
  RNG.cc
  copyright.cc
  events.cc
  partition.cc
  types.cc
 )

add_library(nut ${NUT_SRC} )

target_include_directories(nut PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

install(FILES ${NUT_HEADERS} DESTINATION include)
install(TARGETS nut DESTINATION lib)

# End of file
