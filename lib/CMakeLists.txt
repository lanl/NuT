# CMakeLists for nut lib.
# Want to accomplish two things:
#   1. build the .cc files into a lib,
#   2. move (and configure if necessary) the header files to the include
#      directory.


file(GLOB NUT_HEADERS
   ${CMAKE_CURRENT_SOURCE_DIR}/*.hh
 )


file(GLOB NUT_SRC
   ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
 )

add_library(nut ${NUT_SRC} )

install(FILES ${NUT_HEADERS} DESTINATION include)
install(TARGETS nut DESTINATION lib)

# End of file