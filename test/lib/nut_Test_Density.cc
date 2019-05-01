// T. M. Kelley (c) 2011 LANS LLC

#include "Density.hh"
#include "gtest/gtest.h"
// #include "test_aux.hh"

TEST(Density_tests, instantiate)
{
  typedef float fp_t;
  typedef std::vector<fp_t> v_t;

  size_t n_cells(100);

  v_t rho_nu_e(n_cells);
  v_t rho_nu_e_bar(n_cells);
  v_t rho_nu_x(n_cells);
  v_t rho_nu_x_bar(n_cells);
  v_t rho_p(n_cells);
  v_t rho_n(n_cells);
  v_t rho_e_minus(n_cells);
  v_t rho_e_plus(n_cells);
  v_t rho_A(n_cells);
  v_t y_e(n_cells);
  v_t abar(n_cells);

  nut::Density<fp_t> density(rho_p, rho_e_minus, rho_e_plus, rho_A, y_e, abar);
  EXPECT_TRUE(true);  // seen to fail 1 May 2019 TK
  return;
}  // Density_tests::

// End of file
