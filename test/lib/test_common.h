// test_common.h
// Jun 24, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once
#include "soft_equiv.hh"
#include <array>

template <typename vector_t>
bool
vec_soft_equiv(vector_t const & v,
               vector_t const & v_exp,
               std::string const & s,
               double tol = 5e-14)
{
  bool passed(true);
  constexpr size_t dim{vector_t::dim};
  for(size_t i = 0; i < dim; ++i) {
    bool c_ok = nut::soft_equiv_os(
        v[i], v_exp[i], s + "Component " + std::to_string(i) + " ", tol);
    EXPECT_TRUE(c_ok);
    passed = passed && c_ok;
  }
  return passed;
}

namespace nut_test {

/**\brief A class that looks like, but is not, a random number generator */
template <size_t S, typename geom_t = double>
struct Buffer_RNG {
  std::array<geom_t, S> a_;

  size_t current_;

  Buffer_RNG(geom_t * begin, geom_t * end) : current_(0)
  {
    std::copy(begin, end, a_.begin());
  }

  geom_t & random()
  {
    geom_t & urd = a_[current_];
    current_++;
    if(current_ > a_.size()) { current_ = 0; }
    return urd;
  }

};  // Buffer_RNG

}  // namespace nut_test

// End of file
