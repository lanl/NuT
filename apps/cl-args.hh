// cl-args.hh
// T. M. Kelley
// Jul 20, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef CL_ARGS_HH
#define CL_ARGS_HH

#include <string>
#include <stdint.h>

struct args_t {
  uint32_t n_particles; /* -n */
  std::string inputF;   /* -i --input */
  std::string outputF;  /* -o --output */
  double llimit;        /* -l */
  double ulimit;        /* -u */
  uint32_t chunkSz;     /* -c */
  double dt;            /* -t */
  double alpha;         /* -a */
  uint32_t seed;        /* -s */
};

// parse CL
args_t
parseCL(int argc, char ** argv);

#endif  // include guard

// version
// $Id$

// End of file
