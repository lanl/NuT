// T. M. Kelley (c) 2011 LANS LLC

#define USING_HFBC_SIGMAS
#include "Opacity.hh"
#undef USING_HFBC_SIGMAS
#include "gtest/gtest.h"
#include "meshes/mesh_adaptors/Spherical_Mesh_Interface.h"
#include <vector>

using nut::Cell_Data;
using nut::Opacity;

TEST(nut_opacity, instantiation)
{
  typedef double fp_t;
  // bogus random number gen type: ok to test instantiation
  using vector_t = nut_mesh::Spherical_1D_Mesh::Vector;
  using op_t = Opacity<fp_t, vector_t>;
  using cell_data_t = Cell_Data<fp_t, vector_t>;

  std::vector<cell_data_t> nil{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}}};
  op_t op(std::move(nil));
  EXPECT_TRUE(true);
  return;
}  // test_1

// End of file
