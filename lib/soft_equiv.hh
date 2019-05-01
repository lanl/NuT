// soft_equiv.hh
// T. M. Kelley
// Jan 11, 2011
// Header for soft_equiv
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef SOFT_EQUIV_H
#define SOFT_EQUIV_H

#include <cmath>  // std::abs

namespace nut {
template <typename fp_t>
bool
soft_equiv(fp_t const val, fp_t const ref, fp_t const tol = fp_t(1e-12))
{
  bool equiv = (std::abs(val - ref) < std::abs(ref) * tol);
  /* if ref is zero, or close thereto, compare directly. */
  if(!equiv && (ref < 1e-14)) { equiv = std::abs(val) < tol; }
  return equiv;
}  // soft_equiv
}  // namespace nut

#endif

// version
// $Id$

// End of file
