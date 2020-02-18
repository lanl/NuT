// soft_equiv.hh
// T. M. Kelley
// Jan 11, 2011
// Header for soft_equiv
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef SOFT_EQUIV_H
#define SOFT_EQUIV_H

#include <cmath>  // std::abs
#include <iomanip>
#include <iostream>

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

/**\brief Decide whether a value is within relative error of a reference value.
 * \param value: Value to check
 * \param reference: the reference
 * \param rel_error: the relative error.
 * \remark If the reference is nearly zero, check that the value is nearly zero;
 * also, if the difference is subnormal, you get another chance to agree.  */
template <typename T>
bool
soft_equiv_os(T value,
              T reference,
              std::string const & name,
              T rel_error = 1.0e-15,
              std::ostream & o = std::cerr)
{
  bool ok = std::abs(value - reference) < rel_error * std::abs(reference);
  if(std::abs(reference) < 1.0e-15) { ok = std::abs(value) < rel_error; }
  ok = ok || (std::abs(value - reference) < std::numeric_limits<T>::min());
  if(!ok) {
    o << std::setprecision(12) << "soft_equiv of " << name
      << " failed: expected = " << reference << ", actual = " << value
      << ", tol = " << rel_error << "\n";
  }
  return ok;
}  // soft_equiv_os

}  // namespace nut

#endif

// End of file
