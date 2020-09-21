# nut-meshes

Meshes used by the NuT neutrino Monte Carlo code.


## Prerequisites

There are two prerequisites:
* Google test suite, version 1.8.0 or above

CMake will look for the GTestConfig.cmake files. If either the environment variables `GTEST_ROOT` are defined, they will be passed as a hint to the respective `find_package` calls.

## Unit Tests
Unit tests will be built by default. Run ccmake and look for `ENABLE_UNIT_TESTS` if for some reason you prefer to not have that happen.

## Use by external projects

External projects would include either the `cartesian/Cartesian_Mesh.h` header to use the 3D Cartesian mesh, or the `spherical/Spherical_Mesh.h` header to use the 1D Spherical mesh.
