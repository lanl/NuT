// RNG.cc
// T. M. Kelley
// Sep 20, 2012
// (c) Copyright 2012 LANSLLC, all rights reserved

#include "RNG.hh"

namespace nut {

std::ostream &
operator<<(std::ostream & os, Philox4x32_RNG const & r)
{
  os << "Philox4x32RNG, ctr = {" << r.m_c.v[0] << " " << r.m_c.v[1] << " "
     << r.m_c.v[2] << " " << r.m_c.v[3] << "}"
     << "; key = {" << r.m_k.v[0] << " " << r.m_k.v[1] << "}"
     << "; bank = {" << r.m_bank[0] << " " << r.m_bank[1] << "}"
     << "; currently on " << r.m_idx << ".";
  return os;
}

}  // namespace nut

// End of file
