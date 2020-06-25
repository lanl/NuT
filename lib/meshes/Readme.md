# nut-meshes

Meshes used by the NuT neutrino Monte Carlo code. The meshes have been factored out of the original code to use as Murmeln Mesh API examples.


## Prerequisites

There are two prerequisites:
* Google test suite, version 1.8.0 or above
* Murmeln.

CMake will look for the GTestConfig.cmake and murmelnConfig.cmake files. If either the environment variables `GTEST_ROOT` or `MURMELN_CMAKE_DIR` are defined, they will be passed as a hint to the respective `find_package` calls.

## Unit Tests
Unit tests will be built by default. Run ccmake and look for `ENABLE_UNIT_TESTS` if for some reason you prefer to not have that happen.

## Use by external projects

External projects would include either the `cartesian/Cartesian_Mesh.h` header to use the 3D Cartesian mesh, or the `spherical/Spherical_Mesh.h` header to use the 1D Spherical mesh.
