
#-----------------------------------------------------------------------------#
# Add header files.Note that these will be exported to the parent scope.
# This is what allows cinch to organize this library.
#-----------------------------------------------------------------------------#

target_sources(nut-meshes
  INTERFACE
  ${handy_prefix}/geometry/lorentz.h
 )

#-----------------------------------------------------------------------------#
# Unit tests.
#-----------------------------------------------------------------------------#
target_sources(nut-meshes-unit-tests
  PRIVATE
  ${handy_prefix}/geometry/test/test_lorentz.cc
  )
