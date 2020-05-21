// cl-args.cc
// T. M. Kelley
// Jul 20, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#include "cl-args.hh"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <getopt.h>

args_t
parseCL(int argc, char ** argv)
{
  args_t args;
  args.n_particles = 0;
  args.inputF = "";
  args.outputF = "tally";
  args.llimit = 0.0;
  args.ulimit = 1e12;
  args.chunkSz = 0;
  args.dt = 1e-7;
  args.alpha = 2.0;
  args.seed = 426356;

  static const char * optString = "n:i:o:l:u:c:t:a:s:";

  const struct option longOpts[] = {
      {"n-particles", required_argument, NULL, 'n'},
      {"input", required_argument, NULL, 'i'},
      {"output", required_argument, NULL, 'o'},
      {"lower-limit", required_argument, NULL, 'l'},
      {"upper-limit", required_argument, NULL, 'u'},
      {"chunk-size", required_argument, NULL, 'c'},
      {"simulation-time", required_argument, NULL, 't'},
      {"alpha", required_argument, NULL, 'a'},
      {"seed", required_argument, NULL, 's'},
      {NULL, no_argument, NULL, 0}};

  int longidx(0);
  int opt = getopt_long(argc, argv, optString, longOpts, &longidx);
  // std::cout << __LINE__ << ": opt = " << opt << std::endl;
  while(opt != -1) {
    if(optarg == NULL) { std::cout << "optarg = NULL" << std::endl; }

    std::stringstream s(optarg);
    // std::cout << "optarg = " << optarg << "; stringstream(optarg): "
    //           << s.str() << "; longidx = " << longidx << std::endl;
    switch(opt) {
      case 'n': s >> args.n_particles; break;
      case 'i': args.inputF = optarg; break;
      case 'o': args.outputF = optarg; break;
      case 'l': s >> args.llimit; break;
      case 'u': s >> args.ulimit; break;
      case 'c': s >> args.chunkSz; break;
      case 't': s >> args.dt; break;
      case 'a': s >> args.alpha; break;
      case 's': s >> args.seed; break;
      default:
        std::stringstream errors("parseCL, unrecognized option: ");
        errors << opt << std::endl;
        std::cerr << __LINE__ << ":" << errors.str() << std::endl;
        // throw(std::runtime_error(errors.str()));
    }
    opt = getopt_long(argc, argv, optString, longOpts, &longidx);
  }  // while(opt!= -1);

  // validation:
  std::stringstream errors;
  bool kill(false);
  if(args.ulimit <= args.llimit) {
    errors << "lower limit (" << args.llimit
           << ") must be less than upper limit (" << args.ulimit << ")"
           << std::endl;
    kill = true;
  }
  if(args.llimit < 0.0) {
    errors << "lower limit must be >= 0, setting to 0" << std::endl;
    args.llimit = 0.0;
  }
  if(args.dt < 0.0) {
    errors << "sim time must be >= 0, you gave " << args.dt << std::endl;
    kill = true;
  }
  if(0 == args.chunkSz) {
    if(0 == args.n_particles) {
      errors << "You must transport some particles: use -n <n> with <n> > 0.";
      kill = true;
    }
    args.chunkSz = args.n_particles;
  }

  std::cout << "Command line options:\n"
            << "n: " << args.n_particles << "\n"
            << "input file: " << args.inputF << "\n"
            << "output file: " << args.outputF << "\n"
            << "llimit: " << args.llimit << "\n"
            << "ulimit: " << args.ulimit << "\n"
            << "chunkSize: " << args.chunkSz << "\n"
            << "sim time: " << args.dt << "\n"
            << "alpha: " << args.alpha << "\n"
            << "seed: " << args.seed << "\n";

  if(errors.str() != "") {
    std::cerr << "CL parsing errors:\n" << errors.str() << std::endl;
    if(kill) { exit(-1); }
  }
  return args;
}  // parseCL

// version
// $Id$

// End of file
