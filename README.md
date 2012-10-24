NuT
===

Monte Carlo code for Neutrino Transport

NuT simulates neutrino transport. Currently, it runs in 1D, spherical geometry. Particles events are tallied in the lab frame, while scattering events (cross sections) are calculated in the co-moving frame. Cross-sections are computed via analytic formulae given in Herant et al [2]. We intend to add 3D Cartesian mesh support soon. 

Nut is a C++ analog to the Haskell McPhD code. It currently is serial. We intend to take it parallel in the near future. We will incorporate both distributed (MPI) and SMP (OpenMP) parallelism. Like McPhD, NuT's driver application is currently configured to run only electron neutrinos; mu and tau neutrinos are approximated as nu_x. Changing this is a simple edit.

This compact app captures many of the computational characteristics and challenges of Monte Carlo transport codes. Random number generation is handled by the Philox class of random number generators[1]. We use the Random123 implementation, [available from D. E. Shaw research](http://www.deshawresearch.com/downloads/download_random123.cgi/ "D. E. Shaw Research").

NuT uses the [scons build system](http://www.scons.org/ "SCons"). To get started, please see the file "readme_scons" in the scons directory.

[1]. "Parallel random numbers: as easy as 1,2,3" J. K. Salmon, Mark A. Moraes, Ron O. Dror, David E. Shaw. Proceeding SC '11 Proceedings of 2011 International Conference for High Performance Computing, Networking, Storage and Analysis. doi 10.1145/2063384.2063405

[2]. Herant, M., Benz, W., Hix, W. R., Fryer, C. L., & Colgate, S. A. 1994, ApJ, 435, 339
