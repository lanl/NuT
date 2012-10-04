NuT
===

Monte Carlo code for Neutrino Transport

NuT simulates neutrino transport. Currently, it runs in 1D, spherical geometry. Particles events are tallied in the lab frame, while scattering events (cross sections) are calculated in the co-moving frame. We intend to add 3D Cartesian mesh support soon. 

Nut is a C++ analog to the Haskell McPhD code. It currently is serial. We intend to take it parallel in the near future. We will incorporate both distributed (MPI) and SMP (OpenMP) parallelism.

This compact app captures many of the computational characteristics and challenges of Monte Carlo transport codes. Random numbers generation is handled by the Philox class of random number generators. We use the Random-1-2-3 library, [available from D. E. Shaw research](http://www.deshawresearch.com/downloads/download_random123.cgi/ "D. E. Shaw Research").

NuT uses the [scons build system](http://www.scons.org/ "SCons"). To get started, please see the file "readme_scons" in the scons directory.
