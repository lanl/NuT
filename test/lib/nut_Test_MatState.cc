// T. M. Kelley (c) 2011 LANS LLC

#include "MatState.hh"

#include "expect.hh"
#include "fileio.hh"
#include "gtest/gtest.h"
#include <sstream>
#include <string>

std::string line(
    " 1  1.5000E-02  1.7348E+06  2.7438E+11  3.5650E+08"
    "  8.1058E+00  0.0000E+00  5.1030E+02  1.7185E+02  "
    "5.2054E+21  3.1954E+07  2.2382E+02  3.1954E+07  "
    "2.2248E+02  1.1083E+08  4.7589E+02  1.1083E+08  "
    "4.7589E+02  9.5198E+07  4.7589E+02");

using test_aux::expect;

TEST(Cell_Data, 1D_instantiation_and_initialization)
{
  // cf nut_Test_fileio.cc, test_3:
  typedef double fp_t;
  using vector_t = murmeln_mesh::Spherical_1D_Mesh::Vector;
  using Cell_Data_T = nut::Cell_Data<double, vector_t>;
  typedef nut::MatStateRowP<fp_t, vector_t> row_t;
  typedef std::vector<row_t> vecrow;

  std::stringstream instr(line);
  vecrow rows(nut::read_mat_state_file<fp_t, vector_t>(instr));

  std::vector<Cell_Data_T> state(make_cell_data(rows));

  // check mat state
  // sizes of everything:
  size_t const exp_sz(1u);
  bool sizes_ok = expect(state.size(), exp_sz, "density size");
  EXPECT_TRUE(sizes_ok);

  // values

  // density
  bool density_passed =
      expect(state[0].rho_p, 1.6404204182659538e+35, "rho_p") and
      expect(state[0].rho_e_minus, 1.329691982638017e+36, "rho_e_minus") and
      expect(state[0].rho_e_plus, 0.0, "rho_e_plus") and
      expect(state[0].rho_A, 0.0, "rho_A") and
      expect(state[0].y_e, 8.1058, "y_e") and
      expect(state[0].abar, 0.0, "abar");
  EXPECT_TRUE(density_passed);

  // luminosity
  bool luminosity_passed =
      expect(state[0].l_nue, 1.42784e+59, "nue") and
      expect(state[0].l_nueb, 1.42784e+59, "nuebar") and
      expect(state[0].l_nux, 9.5198000000000001e+58, "nux");
  EXPECT_TRUE(luminosity_passed);

  // temperature
  bool temperature_passed =
      expect(state[0].T_p, 43.974286020000001, "T_p") and
      expect(state[0].T_e_minus, 43.974286020000001, "T_e_minus") and
      expect(state[0].T_e_plus, 43.974286020000001, "T_e_plus");
  EXPECT_TRUE(temperature_passed);

  // velocity
  bool velocity_passed = expect(state[0].velocity, {3.565e+08}, "v");
  EXPECT_TRUE(velocity_passed);
  return;
}

// End of file
