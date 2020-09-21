// soft_equiv.h
// T. M. Kelley
// Feb 20, 2019
// (c) Copyright 2019 LANSLLC, all rights reserved

#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace nut_mesh {

/**\brief Decide whether a value is within relative error of a reference value.
 * \param value: Value to check
 * \param reference: the reference
 * \param rel_error: the relative error.
 * \remark If the reference is nearly zero, check that the value is nearly zero;
 * also, if the difference is subnormal, you get another chance to agree.  */
template <typename T>
bool soft_equiv(T value, T reference, T rel_error = 1.0e-15) {
  bool ok = std::abs(value - reference) < rel_error * std::abs(reference);
  if (std::abs(reference) < 1.0e-15) {
    ok = std::abs(value) < rel_error;
  }
  ok = ok || (std::abs(value - reference) < std::numeric_limits<T>::min());
  return ok;
}

/**\brief Decide whether a value is within relative error of a reference value.
 * \param value: Value to check
 * \param reference: the reference
 * \param rel_error: the relative error.
 * \remark If the reference is nearly zero, check that the value is nearly zero;
 * also, if the difference is subnormal, you get another chance to agree.  */
template <typename T>
bool soft_equiv_os(T value, T reference, std::string const &name,
                   T rel_error = 1.0e-15, std::ostream &o = std::cerr) {
  bool ok = std::abs(value - reference) < rel_error * std::abs(reference);
  if (std::abs(reference) < 1.0e-15) {
    ok = std::abs(value) < rel_error;
  }
  ok = ok || (std::abs(value - reference) < std::numeric_limits<T>::min());
  if (!ok) {
    o << std::setprecision(12) << "soft_equiv of " << name
      << " failed: expected = " << reference << ", actual = " << value
      << ", tol = " << rel_error << "\n";
  }
  return ok;
}

}  // namespace nut_mesh

// End of file
