NuT
===

Monte Carlo code for Neutrino Transport

NuT simulates neutrino transport. Currently, it runs in 1D, spherical geometry. Particles events are tallied in the lab frame, while scattering events (cross sections) are calculated in the co-moving frame. Cross-sections are computed via analytic formulae given in Herant et al [2]. We intend to add 3D Cartesian mesh support soon.

Nut is a C++ analog to the Haskell McPhD code. It currently is serial. We intend to take it parallel in the near future. We will incorporate both distributed (MPI) and SMP (OpenMP) parallelism. Like McPhD, NuT's driver application is currently configured to run only electron neutrinos; mu and tau neutrinos are approximated as nu_x. Changing this is a simple edit.

This compact app captures many of the computational characteristics and challenges of Monte Carlo transport codes. Random number generation is handled by the Philox class of random number generators[1]. We use the Random123 implementation, [available from D. E. Shaw research](http://www.deshawresearch.com/downloads/download_random123.cgi/ "D. E. Shaw Research").

NuT uses the [scons build system](http://www.scons.org/ "SCons"). To get started, please see the file "readme_scons" in the scons directory.

Los Alamos National Security, LLC (LANS) owns the copyright to NuT, which it identifies internally as LA-CC-11-087. The license is BSD 3-Clause.


Repository Structure
--------------------

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
They are massless, or very nearly so.
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
At the beginning of the time step, a certain amount of energy is going to be emitted as neutrinos.
The energy varies from mesh cell to mesh cell.
 NuT divides this energy into particles.
Then it transports each particle using the Monte Carlo method.

The method is simple to state: until the particle is dead, gone, or its time is up,
1. decide the next step,
2. apply the event by advancing the particle to the event site, update its clock, and change the particle state,
3. make notes about what happened (tally).

Some of the outcomes from a Monte Carlo step include the particle being absorbed, the particle escaping from the region of interest, or the particle reaching the end of the time step.
The different interactions are relatively easy to encode.
The famous Monte Carlo bit is that steps are chosen probabilistically:
a probability is calculated for each event happening, and then one uses a (pseudo)random number to pick one of the events.
It's a lovely method, invented by Stanislaw Ulam and Nick Metropolis at Los Alamos.


[1]. "Parallel random numbers: as easy as 1,2,3" J. K. Salmon, Mark A. Moraes, Ron O. Dror, David E. Shaw. Proceeding SC '11 Proceedings of 2011 International Conference for High Performance Computing, Networking, Storage and Analysis. doi 10.1145/2063384.2063405

[2]. Herant, M., Benz, W., Hix, W. R., Fryer, C. L., & Colgate, S. A. 1994, ApJ, 435, 339
