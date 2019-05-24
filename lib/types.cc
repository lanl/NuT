// types.cc
// T. M. Kelley
// May 23, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#include "types.hh"
#include <string>

namespace nut {

seed_t
species_seed(Species const s)
{
  seed_t n(-1);
  switch(s) {
    case nu_e: n = 0; break;
    case nu_e_bar: n = 1; break;
    case nu_x: n = 6; break;
    case nu_mu: n = 2; break;
    case nu_mu_bar: n = 3; break;
    case nu_tau: n = 4; break;
    case nu_tau_bar: n = 5; break;
  }
  return n;
}

std::string
species_name(Species const s)
{
  std::string n;
  switch(s) {
    case nu_e: n = "nu_e"; break;
    case nu_e_bar: n = "nu_e_bar"; break;
    case nu_x: n = "nu_x"; break;
    case nu_mu: n = "nu_mu"; break;
    case nu_mu_bar: n = "nu_mu_bar"; break;
    case nu_tau: n = "nu_tau"; break;
    case nu_tau_bar: n = "nu_tau_bar"; break;
  }
  return n;
}

Species
species_type(std::string const & n)
{
  Species s;
  if(n == "nu_e")
    s = nu_e;
  else if(n == "nu_e_bar")
    s = nu_e_bar;
  else if(n == "nu_x")
    s = nu_x;
  else if(n == "nu_mu")
    s = nu_mu;
  else if(n == "nu_mu_bar")
    s = nu_mu_bar;
  else if(n == "nu_tau")
    s = nu_tau;
  else
    s = nu_tau_bar;
  return s;
}

}  // namespace nut

// version
// $Id$

// End of file
