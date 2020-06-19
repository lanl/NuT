file(GLOB bh-3_support_SOURCES
    cl-args.cc
 )

message("PROJECT_SOURCE_DIR = " ${PROJECT_SOURCE_DIR} )
message("RANDOM123_INCLUDE_DIR = " ${RANDOM123_INCLUDE_DIR})

add_library(bh-3_support ${bh-3_support_SOURCES})

target_include_directories(bh-3_support
  PUBLIC
    "${PROJECT_SOURCE_DIR}/lib"
    "${PROJECT_SOURCE_DIR}/lib/meshes"
  )

add_executable(bh-3 bh-3.cc)

target_include_directories(bh-3
  PRIVATE
    ${RANDOM123_INCLUDE_DIR}
  )

if(NUT_ENABLE_ASSERTIONS)
  target_compile_definitions(bh-3
    PRIVATE
      REQUIRE_ON=1)
endif()
target_link_libraries(bh-3 bh-3_support nut nut-meshes ${Murmeln_LIB})

add_executable(bh-3-cart bh-3-cartesian-3D.cc)

target_include_directories(bh-3-cart
  PRIVATE
    ${RANDOM123_INCLUDE_DIR}
  )
if(NUT_ENABLE_ASSERTIONS)
  target_compile_definitions(bh-3-cart
    PRIVATE
      REQUIRE_ON=1)
endif()
target_link_libraries(bh-3-cart bh-3_support nut ${Murmeln_LIB})

# Build with DBC assertions turned on.
# add_definitions(-DREQUIRE_ON)

