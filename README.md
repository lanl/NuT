# NuT
===

[![Build Status](https://travis-ci.org/lanl/NuT.svg?branch=feat-add-gtest)](https://travis-ci.org/lanl/NuT)

Monte Carlo code for Neutrino Transport

NuT simulates neutrino transport. Currently, it runs in 1D, spherical and 3D Cartesian geometries (the mesh geometry and topology is carefully encapsulated, so adding another mesh requires a few hundred lines of interface implementation). Particles events are tallied in the lab frame, while scattering events (cross sections) are calculated in the co-moving frame. Cross-sections are computed via analytic formulae given in Herant et al [2].

Nut is a C++ analog to the Haskell McPhD code. It currently has serial (master branch) and OpenMP (openmp branch) implementations. NuT is principally aimed at exploring on-node parallelism and performance issues. Like McPhD, NuT's driver application is currently configured to run electron neutrinos and antineutrinos; mu and tau neutrinos are approximated as nu_x. Changing this is a simple edit (give me a star and I'll add it).

This compact app captures many of the computational characteristics and challenges of Monte Carlo transport codes. Random number generation is handled by the Philox class of random number generators[1]. We use the Random123 implementation, [available from D. E. Shaw research](http://www.deshawresearch.com/downloads/download_random123.cgi/ "D. E. Shaw Research").

NuT uses the [cmake build system](http://cmake.org/ "CMake").

Los Alamos National Security, LLC (LANS) owns the copyright to NuT, which it identifies internally as LA-CC-11-087. The license is BSD 3-Clause.

Quick Start
===========
There is one external dependency: Random123, obtain from D. E. Shaw:
   https://www.deshawresearch.com/resources_random123.html
Random123 is implemented entirely in header files, so installation is minimal.

1. Specify environment variables GTEST_DIR, RANDOM123_DIR, CC, and CXX. RANDOM123_DIR
should point to the root of the installation. For example,

   export GTEST_DIR=/path/to/gtest 

   export RANDOM123_DIR=/home/me/downloads/deshaw/Random123-1.08

   export CC=/home/me/downloads/llvm/clang+llvm-3.7.0/bin/clang

   export CXX=/home/me/downloads/llvm/clang+llvm-3.7.0/bin/clang++

2. Under the root NuT directory, create a directory called build:

   me@superMachine:~/dev/nut$ mkdir build

   me@superMachine:~/dev/nut$ cd build

3. Configure and build:

   me@superMachine:~/dev/nut/build$ cmake -DCMAKE_INSTALL_PREFIX=./nut ..

   me@superMachine:~/dev/nut/build$ make VERBOSE=on -j 4 2>&1 | tee -a make.out

4. Run unit tests

   me@superMachine:~/dev/nut/build$ ./test/nut_unittests

... tests all pass!




Repository Structure
====================

lib: headers and source files to build libnut.

test: unit test and integrated test code and data
  * lib: unit tests (uses ancient, simple, home-grown test suite)
  * data: some input files with snapshots from a black hole progenitor (provided by C.Fryer))

apps: bh-3, an application that drives libnut.


The BH-3 application
====================

NuT includes a small application to demonstrate the library, called bh-3.
"bh" is for black hole---the application models the neutrino transport given a snapshot of the material state of a star collapsing to a black hole.

bh-3 configures itself from command line variables, reads in a material state from a file, then transports the indicated set of particles.
It carries out transport for electron and anti-electron neutrinos; everything else is lumped together in a category called &nu;<sub>x</sub>.
A tally is generated for each species.



Quick Introduction to Monte Carlo Neutrino Transport
====================================================

Neutrinos are particles associated with the weak nuclear force.
They are very nearly massless; we treat them as massless; i.e. they move at the speed of light.
Note that neutrinos are very different from neutrons.
Neutrinos are conventionally denoted by the Greek letter &nu; (thus nu transport, or NuT).
Neutrinos come in three flavors, corresponding to the three types of lepton: electron (&nu;<sub>e</sub>), mu (&nu;<sub>&mu;</sub>), and tau (&nu;<sub>&tau;</sub>).
And since neutrinos are fermions, each comes with a corresponding anti-neutrino, usually written with an over bar.

Neutrinos do lots of interesting things; our particular interest is in the role ofneutrinos in  core-collapse supernovae.
Neutrinos play an essential role in carrying energy outward from the collapsing star.
In the heart of a proto-neutron star, the matter is so dense that the neutrinos move diffusively, essentially losing all sense of direction immediately.
Far out in the stellar envelope, neutrinos are essentially streaming, only occasionally being attenuated.
In between these two extremes lies the transport regime, where neutrinos are scattering and exchanging energy and momentum with the star's matter in complex, history dependent ways.

One way to better understand this regime is to simulate it.
NuT is a tool for simulating neutrino transport using Monte Carlo.
In the Monte Carlo method, we partition the neutrino energy into a set of discrete chunks called particles, and track the particles through the material medium with which they interact.
There are a number of interactions that a neutrino can have with matter, such as being scattered by various particles or being absorbed by a nucleus.

NuT has a concept of time stepping and a mesh.
The mesh is just a discretization of space that makes it easy to keep track of what is where.
The NuT mesh encapsulates the geometry and topology of the mesh from the rest of the code.
During the time step, the properties of the material medium are held constant.
At the beginning of the time step, a certain amount of energy is going to be emitted as neutrinos.
Additionally, there may be an initial population, or census, of particles left over from the preceding time step.
The energy to be emitted varies from mesh cell to mesh cell.
NuT divides this energy into particles.
Then it transports each particle using the Monte Carlo method.

The method is simple to state: until the particle is dead, gone, or its time is up,

1. decide the next step,
2. apply the event by advancing the particle to the event site, update its clock, and change the particle state,
3. make notes about what happened (tally).

Some of the outcomes from a Monte Carlo step include the particle being absorbed, the particle escaping from the problem domain, or the particle reaching the end of the time step.
The different interactions are relatively easy to encode.
The famous Monte Carlo bit is that steps are chosen probabilistically:
a probability is calculated for each event happening, and then one uses a (pseudo)random number to pick one of the events.


[1]. "Parallel random numbers: as easy as 1,2,3" J. K. Salmon, Mark A. Moraes, Ron O. Dror, David E. Shaw. Proceeding SC '11 Proceedings of 2011 International Conference for High Performance Computing, Networking, Storage and Analysis. doi 10.1145/2063384.2063405

[2]. Herant, M., Benz, W., Hix, W. R., Fryer, C. L., & Colgate, S. A. 1994, ApJ, 435, 339
