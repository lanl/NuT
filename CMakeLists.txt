# Top-level cmake file for NuT project

# Quick start:
# Dependencies: gtest 1.7.0
# 1. specify environment variables GTEST_DIR and BOOST_DIR. These
# should point to the root of the respective installations. For example,
#    export GTEST_DIR=/home/me/downloads/google-test/gtest-1.7.0
# 2. Under the root NuT directory, create a directory called build:
#    me@superMachine:~/dev/nut$ mkdir build
#    me@superMachine:~/dev/nut$ cd build
# 3. configure and build, including the Doxygen documentation:
#    me@superMachine:~/dev/nut/build$ cmake -Dbuild_api_doc=1
# -DCMAKE_INSTALL_PREFIX=./nut ..
#    me@superMachine:~/dev/nut/build$ make VERBOSE=on -j 4 2>&1 | tee -a make.out
# 4. run unit tests
#    me@superMachine:~/dev/nut/build$ ./test/nut_unittests
# ... tests all pass!
# 5. install
#    me@superMachine:~/dev/nut/build$ make install
# ... files installed to CMAKE_INSTALL_PREFIX

# To do:
#  1. incorporate build and version information into destination
#  2. incorporate CTest


# To remove all cmake added files and directories (after initial 'cmake .'):
# > make clean-cmake
# Note: This also removes all project build products as they are under dir CMakeFiles.
#       This is an out of src build process (targets are not built in src dirs).
#
#

#------------------------------------------------------------


# ensure min cmake version, optional
cmake_minimum_required (VERSION 3.0)

# name the project
# languages - CXX and C are included by default
project (NuT LANGUAGES CXX)

# set compiler, looking first to the environment, defaulting to
# GNU if nothing better is found
if (DEFINED ENV{CC})
   set(CMAKE_C_COMPILER $ENV{CC})
   set(CMAKE_CXX_COMPILER $ENV{CXX})
   # set(CMAKE_FORTRAN_COMPILER $ENV{F90})
else()
   set(CMAKE_C_COMPILER gcc)
   set(CMAKE_CXX_COMPILER g++)
   # set(CMAKE_FORTRAN_COMPILER gfortran)
endif()


# create variable
set(NuT_COMPILE_FLAGS "-std=c++11 -g -Wall ")

# append CMAKE_CXX_FLAGS w/ above variable
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NuT_COMPILE_FLAGS}")



#--------- Random 123 --------------------------
if (NOT DEFINED RANDOM123_ROOT)

  if (DEFINED ENV{RANDOM123_DIR})

    set(RANDOM123_ROOT $ENV{RANDOM123_DIR})
    message("RANDOM123_ROOT is ${RANDOM123_ROOT}")
    include_directories(
      ${RANDOM123_ROOT}/include
    )

  else()
    # j cohn 21 Jul 2015 - though I thought I had added this earlier
        # RANDOM123_INCLUDES is set by loading an appropriate boost module on darwin, for example
    if (DEFINED ENV{RANDOM123_INCLUDES})

                set(RANDOM123_ROOT $ENV{RANDOM123_INCLUDES})
                message("RANDOM123_ROOT is ${RANDOM123_ROOT}")
                include_directories(
                  ${RANDOM123_ROOT}/include
                )

        else()

                message("Warning: environment variable RANDOM123_DIR is not defined")
                message("nut_unittests will not be made")
                message("To make nut_unittests, set RANDOM123_DIR to top level of Random123 distribution directory (the one with docs, examples, include).")
                return()

    endif()

  endif()

endif()













# add a local directory 'lib' to the include search path
# # (This is for building ¬¬NuT unit tests)
# include_directories(include)
#
# # recurse into the include subdirectory: this is to pick up the installation
# add_subdirectory(include)
#
add_subdirectory(lib)

#------Cmake clean up-------------------------------
# set(NuT_CMAKE_CLEAN_FILES CMakeCache.txt cmake_install.cmake Makefile Doxyfile CMakeFiles)
# set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${NuT_CMAKE_CLEAN_FILES}")
set(CMAKE_CLEAN_CMD COMMAND ${CMAKE_COMMAND} -E remove -f CMakeCache.txt cmake_install.cmake Makefile Doxyfile
    COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles)

if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
   # in SRC build
   # delete cmake created directory(ies)
   set(CMAKE_CLEAN_CMD ${CMAKE_CLEAN_CMD}
     COMMAND ${CMAKE_COMMAND} -E remove_directory doc/Doxygen)
else()
   # out of SRC build
   # delete cmake created directories
   set(CMAKE_CLEAN_CMD ${CMAKE_CLEAN_CMD}
     COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/*)
endif()

add_custom_target(clean-cmake
     ${CMAKE_CLEAN_CMD})


#----------Google Tests------------------------------
# default behaviour is to build gtests
if (NOT DEFINED build_gtests)
   set(build_gtests "yes")
endif()
# test for option "-Dbuild_gtests=yes"  (-D variable is retained from prior runs)
if (${build_gtests} )
   set(NuT_GTEST_DIR ${CMAKE_SOURCE_DIR}/test)

   if (IS_DIRECTORY ${NuT_GTEST_DIR})
      add_subdirectory(${NuT_GTEST_DIR})
   else()
      message("Warning: Test dir ${NuT_GTEST_DIR} not found.")
   endif()
endif()



# Put this down here to catch BOOST configurizing
add_subdirectory(apps)


#----------Doxygen------------------------------
if (NOT DEFINED build_api_doc)
   set(build_api_doc "no")
endif()
if (${build_api_doc} )#
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
     message("Configuring Doxygen")
     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
     # doxygen can't make a sub-dir on the fly, so we're pre-making it here
     #   this should correspond to the OUTPUT_DIRECTORY in the Doxyfile.in
     #   grep OUTPUT_DIR Doxyfile | grep -v ^# | cut -d '=' -f 2
     file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxygen)
     add_custom_target(doc ALL
       COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
       COMMENT "Generating API documentation with Doxygen" VERBATIM
     )
  endif(DOXYGEN_FOUND)

endif() #${build_api_doc}


