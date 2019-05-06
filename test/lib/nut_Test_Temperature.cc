// T. M. Kelley (c) 2011 LANS LLC

#include "Temperature.hh"
#include "gtest/gtest.h"
#include "test_aux.hh"

TEST(nut_Temperature, instantiate)
{
  bool passed(true);

  typedef float fp_t;
  typedef std::vector<fp_t> v_t;

  size_t n_cells(100);

  v_t T_p(n_cells);
  v_t T_e_minus(n_cells);
  v_t T_e_plus(n_cells);

  nut::Temperature<fp_t> temperature(T_p, T_e_minus, T_e_plus);
  EXPECT_TRUE(passed);
  return;
}

// End of file
