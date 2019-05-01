// T. M. Kelley (c) 2011 LANS LLC

#define USING_HFBC_SIGMAS
#include "Opacity.hh"
#undef USING_HFBC_SIGMAS

#include "Density.hh"
#include "Temperature.hh"
#include "gtest/gtest.h"

using nut::Opacity;

TEST(nut_opacity, instantiation)
{
  typedef double fp_t;
  // bogus random number gen type: ok to test instantiation
  typedef Opacity<fp_t> op_t;
  typedef nut::Density<fp_t> rho_t;
  typedef nut::Temperature<fp_t> temp_t;

  std::vector<fp_t> nil(10, 0);
  rho_t rho(nil, nil, nil, nil, nil, nil);
  temp_t T(nil, nil, nil);
  op_t op(rho, T);
  EXPECT_TRUE(true);
  return;
}  // test_1

// End of file
