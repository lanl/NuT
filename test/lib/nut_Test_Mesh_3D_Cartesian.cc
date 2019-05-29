// T. M. Kelley (c) 2011 LANS LLC

#include "Mesh_3D_Cartesian.hh"
#include "RNG.hh"
#include "expect.hh"
#include "gtest/gtest.h"

#ifdef HAVE_MURMELN

using test_aux::expect;
using Mesh_Iface_T = murmeln::Cartesian_Mesh_Interface;
using Mesh_T = Mesh_Iface_T::Mesh;
using Cell_T = Mesh_Iface_T::Cell;
using Index_T = Mesh_Iface_T::Index_T;
using Geom_T = Mesh_Iface_T::Geom_T;
using Point = Mesh_Iface_T::Point;

// These are specific to Cartesian_Mesh
using Face = Mesh_Iface_T::Face;
Face const f_l_x{Mesh_T::LOW_X};
Face const f_l_y{Mesh_T::LOW_Y};
Face const f_l_z{Mesh_T::LOW_Z};
Face const f_h_x{Mesh_T::HIGH_X};
Face const f_h_y{Mesh_T::HIGH_Y};
Face const f_h_z{Mesh_T::HIGH_Z};

// helper function
namespace std {
std::ostream &
operator<<(std::ostream & s, Cell_T const & c)
{
  s << c.id();
  return s;
}
}  // namespace std

TEST(nut_mesh_3D_Cartesian, instantiate)
{
  Index_T const nx{1};
  Index_T const ny{1};
  Index_T const nz{1};
  // vbd bds;
  // nut::mkReflectBCs<Mesh_T, Cell_T>(bds, nx, ny, nz);
  Mesh_Iface_T mesh(nx, ny, nz);
  EXPECT_TRUE(true);
  return;
}

template <typename Index_T>
bool
test_2_core(Index_T const nx, Index_T const ny, Index_T const nz)
{
  // vbd bds;
  // nut::mkReflectBCs<Mesh_T, Cell_T>(bds, nx, ny, nz);
  Mesh_Iface_T mesh(nx, ny, nz, 1.0, 1.0, 1.0);
  Index_T const n_cs(nx * ny * nz);
  bool passed = n_cs == mesh.num_cells();
  return passed;
}

TEST(nut_mesh_3D_Cartesian, instantiate_correctly_n_cells)
{
  Index_T const nx = 3;
  Index_T const ny = 5;
  Index_T const nz = 7;
  EXPECT_TRUE(test_2_core(nx, ny, nz));
  return;
}  // test_2

// for any 3D cartesian mesh, the volume should be constant
// (and simple!)
template <typename Index_T>
bool
test_3_core(Geom_T const dx,
            Index_T const nx,
            Geom_T const dy,
            Index_T const ny,
            Geom_T const dz,
            Index_T const nz)
{
  bool passed(true);
  Mesh_T mesh(nx, ny, nz, dx, dy, dz);
  Geom_T const exp_vol(dx * dy * dz);
  for(Index_T cidx = 0; cidx < mesh.num_cells(); ++cidx) {
    Cell_T const c{cidx};
    Geom_T const volume = mesh.volume(c);
    bool cok = volume == exp_vol;
    passed = passed && cok;
  }
  return passed;
}

TEST(nut_mesh_3D_Cartesian, instantiate_correctly_volume)
{
  Geom_T const dx(2.0);
  Geom_T const dy(3.0);
  Geom_T const dz(5.0);
  Index_T const nx = 3;
  Index_T const ny = 5;
  Index_T const nz = 7;

  bool passed = test_3_core(dx, nx, dy, ny, dz, nz);
  EXPECT_TRUE(passed);
  return;
}  // test_3

TEST(nut_mesh_3D_Cartesian, is_boundary_3_cell_line)
{
  using test_aux::expect;
  Geom_T const dx(2.0);
  Geom_T const dy(3.0);
  Geom_T const dz(5.0);

  Index_T const nx = 3;
  Index_T const ny = 1;
  Index_T const nz = 1;

  Mesh_Iface_T mesh(nx, ny, nz, dx, dy, dz);

  {
    Cell_T c(0);
    bool passedc0 = expect(mesh.is_boundary(c, f_l_x), true, "c0 LOW_X") &&
                    expect(mesh.is_boundary(c, f_h_x), false, "c0 HIGH_X") &&
                    expect(mesh.is_boundary(c, f_l_y), true, "c0 LOW_Y") &&
                    expect(mesh.is_boundary(c, f_h_y), true, "c0 HIGH_Y") &&
                    expect(mesh.is_boundary(c, f_l_z), true, "c0 LOW_Z") &&
                    expect(mesh.is_boundary(c, f_h_z), true, "c0 HIGH_Z");
    EXPECT_TRUE(passedc0);
  }
  {
    Cell_T c(1);
    bool passedc1 = expect(mesh.is_boundary(c, f_l_x), false, "c1 LOW_X") &&
                    expect(mesh.is_boundary(c, f_h_x), false, "c1 HIGH_X") &&
                    expect(mesh.is_boundary(c, f_l_y), true, "c1 LOW_Y") &&
                    expect(mesh.is_boundary(c, f_h_y), true, "c1 HIGH_Y") &&
                    expect(mesh.is_boundary(c, f_l_z), true, "c1 LOW_Z") &&
                    expect(mesh.is_boundary(c, f_h_z), true, "c1 HIGH_Z");
    EXPECT_TRUE(passedc1);
  }
  {
    Cell_T c(2);
    bool passedc2 = expect(mesh.is_boundary(c, f_l_x), false, "c2 LOW_X") &&
                    expect(mesh.is_boundary(c, f_h_x), true, "c2 HIGH_X") &&
                    expect(mesh.is_boundary(c, f_l_y), true, "c2 LOW_Y") &&
                    expect(mesh.is_boundary(c, f_h_y), true, "c2 HIGH_Y") &&
                    expect(mesh.is_boundary(c, f_l_z), true, "c2 LOW_Z") &&
                    expect(mesh.is_boundary(c, f_h_z), true, "c2 HIGH_Z");
    EXPECT_TRUE(passedc2);
  }
  return;
}  // test_4

TEST(nut_mesh_3D_Cartesian, is_boundary_6_cell_cluster)
{
  using test_aux::expect;
  Geom_T const dx(2.0);
  Geom_T const dy(3.0);
  Geom_T const dz(5.0);

  Index_T const nx = 1;
  Index_T const ny = 3;
  Index_T const nz = 2;

  Mesh_Iface_T mesh(nx, ny, nz, dx, dy, dz);

  {
    Cell_T c(0);  // (0,0,0)
    bool passedc0 = expect(mesh.is_boundary(c, f_l_x), true, "c0 LOW_X") &&
                    expect(mesh.is_boundary(c, f_h_x), true, "c0 HIGH_X") &&
                    expect(mesh.is_boundary(c, f_l_y), true, "c0 LOW_Y") &&
                    expect(mesh.is_boundary(c, f_h_y), false, "c0 HIGH_Y") &&
                    expect(mesh.is_boundary(c, f_l_z), true, "c0 LOW_Z") &&
                    expect(mesh.is_boundary(c, f_h_z), false, "c0 HIGH_Z");
    EXPECT_TRUE(passedc0);
  }
  {
    Cell_T c(1);  // (0,1,0)
    bool passedc1 = expect(mesh.is_boundary(c, f_l_x), true, "c1 LOW_X") &&
                    expect(mesh.is_boundary(c, f_h_x), true, "c1 HIGH_X") &&
                    expect(mesh.is_boundary(c, f_l_y), false, "c1 LOW_Y") &&
                    expect(mesh.is_boundary(c, f_h_y), false, "c1 HIGH_Y") &&
                    expect(mesh.is_boundary(c, f_l_z), true, "c1 LOW_Z") &&
                    expect(mesh.is_boundary(c, f_h_z), false, "c1 HIGH_Z");
    EXPECT_TRUE(passedc1);
  }
  {
    Cell_T c(2);  // (0,2,0)
    bool passedc2 = expect(mesh.is_boundary(c, f_l_x), true, "c2 LOW_X") &&
                    expect(mesh.is_boundary(c, f_h_x), true, "c2 HIGH_X") &&
                    expect(mesh.is_boundary(c, f_l_y), false, "c2 LOW_Y") &&
                    expect(mesh.is_boundary(c, f_h_y), true, "c2 HIGH_Y") &&
                    expect(mesh.is_boundary(c, f_l_z), true, "c2 LOW_Z") &&
                    expect(mesh.is_boundary(c, f_h_z), false, "c2 HIGH_Z");
    EXPECT_TRUE(passedc2);
  }
  {
    Cell_T c(3);  // (0,0,1)
    bool passedc3 = expect(mesh.is_boundary(c, f_l_x), true, "c3 LOW_X") &&
                    expect(mesh.is_boundary(c, f_h_x), true, "c3 HIGH_X") &&
                    expect(mesh.is_boundary(c, f_l_y), true, "c3 LOW_Y") &&
                    expect(mesh.is_boundary(c, f_h_y), false, "c3 HIGH_Y") &&
                    expect(mesh.is_boundary(c, f_l_z), false, "c3 LOW_Z") &&
                    expect(mesh.is_boundary(c, f_h_z), true, "c3 HIGH_Z");
    EXPECT_TRUE(passedc3);
  }
  {
    Cell_T c(4);  // (0,1,1)
    bool passedc4 = expect(mesh.is_boundary(c, f_l_x), true, "c4 LOW_X") &&
                    expect(mesh.is_boundary(c, f_h_x), true, "c4 HIGH_X") &&
                    expect(mesh.is_boundary(c, f_l_y), false, "c4 LOW_Y") &&
                    expect(mesh.is_boundary(c, f_h_y), false, "c4 HIGH_Y") &&
                    expect(mesh.is_boundary(c, f_l_z), false, "c4 LOW_Z") &&
                    expect(mesh.is_boundary(c, f_h_z), true, "c4 HIGH_Z");
    EXPECT_TRUE(passedc4);
  }
  {
    Cell_T c(5);  // (0,2,1)
    bool passedc5 = expect(mesh.is_boundary(c, f_l_x), true, "c5 LOW_X") &&
                    expect(mesh.is_boundary(c, f_h_x), true, "c5 HIGH_X") &&
                    expect(mesh.is_boundary(c, f_l_y), false, "c5 LOW_Y") &&
                    expect(mesh.is_boundary(c, f_h_y), true, "c5 HIGH_Y") &&
                    expect(mesh.is_boundary(c, f_l_z), false, "c5 LOW_Z") &&
                    expect(mesh.is_boundary(c, f_h_z), true, "c5 HIGH_Z");
    EXPECT_TRUE(passedc5);
  }
  return;
}  // test_5

/* This test is different from the original in that the default is now vacuum
 * bdy cond. */
TEST(nut_mesh_3D_Cartesian, cell_across_face_6_cell_cluster)
{
  Geom_T const dx(2.0);
  Geom_T const dy(3.0);
  Geom_T const dz(5.0);

  Index_T const nx = 1;
  Index_T const ny = 3;
  Index_T const nz = 2;

  Mesh_Iface_T mesh(nx, ny, nz, dx, dy, dz);
  Point p{0.0, 0.0, 0.0};

  auto f = [&mesh, &p](Cell_T c, Face const f, Cell_T const e, const char * s) {
    return expect(mesh.cell_across(c, f, p), e, s);
  };

  Cell_T const null_cell = mesh.null_cell();

  {
    Cell_T c(0);  // (0,0,0)
    bool passedc01 = f(c, f_l_x, null_cell, "c0 low x");
    bool passedc02 = f(c, f_h_x, null_cell, "c0 high x");
    bool passedc03 = f(c, f_l_y, null_cell, "c0 low y");
    bool passedc04 = f(c, f_h_y, Cell_T(1), "c0 high y");
    bool passedc05 = f(c, f_l_z, null_cell, "c0 low z");
    bool passedc06 = f(c, f_h_z, Cell_T(3), "c0 high z");
    EXPECT_TRUE(passedc01);
    EXPECT_TRUE(passedc02);
    EXPECT_TRUE(passedc03);
    EXPECT_TRUE(passedc04);
    EXPECT_TRUE(passedc05);
    EXPECT_TRUE(passedc06);
  }
  {
    Cell_T c(1);
    bool passedc11 = f(c, f_l_x, null_cell, "c1 low x");
    bool passedc12 = f(c, f_h_x, null_cell, "c1 high x");
    bool passedc13 = f(c, f_l_y, Cell_T(0), "c1 low y");
    bool passedc14 = f(c, f_h_y, Cell_T(2), "c1 high y");
    bool passedc15 = f(c, f_l_z, null_cell, "c1 low z");
    bool passedc16 = f(c, f_h_z, Cell_T(4), "c1 high z");
    EXPECT_TRUE(passedc11);
    EXPECT_TRUE(passedc12);
    EXPECT_TRUE(passedc13);
    EXPECT_TRUE(passedc14);
    EXPECT_TRUE(passedc15);
    EXPECT_TRUE(passedc16);
  }
  {
    Cell_T c(2);
    bool passedc21 = f(c, f_l_x, null_cell, "c2 low x");
    bool passedc22 = f(c, f_h_x, null_cell, "c2 high x");
    bool passedc23 = f(c, f_l_y, Cell_T(1), "c2 low y");
    bool passedc24 = f(c, f_h_y, null_cell, "c2 high y");
    bool passedc25 = f(c, f_l_z, null_cell, "c2 low z");
    bool passedc26 = f(c, f_h_z, Cell_T(5), "c2 high z");
    EXPECT_TRUE(passedc21);
    EXPECT_TRUE(passedc22);
    EXPECT_TRUE(passedc23);
    EXPECT_TRUE(passedc24);
    EXPECT_TRUE(passedc25);
    EXPECT_TRUE(passedc26);
  }
  {
    Cell_T c(3);
    bool passedc31 = f(c, f_l_x, null_cell, "c3 low x");
    bool passedc32 = f(c, f_h_x, null_cell, "c3 high x");
    bool passedc33 = f(c, f_l_y, null_cell, "c3 low y");
    bool passedc34 = f(c, f_h_y, Cell_T(4), "c3 high y");
    bool passedc35 = f(c, f_l_z, Cell_T(0), "c3 low z");
    bool passedc36 = f(c, f_h_z, null_cell, "c3 high z");
    EXPECT_TRUE(passedc31);
    EXPECT_TRUE(passedc32);
    EXPECT_TRUE(passedc33);
    EXPECT_TRUE(passedc34);
    EXPECT_TRUE(passedc35);
    EXPECT_TRUE(passedc36);
  }
  {
    Cell_T c(4);
    bool passedc41 = f(c, f_l_x, null_cell, "c4 low x");
    bool passedc42 = f(c, f_h_x, null_cell, "c4 high x");
    bool passedc43 = f(c, f_l_y, Cell_T(3), "c4 low y");
    bool passedc44 = f(c, f_h_y, Cell_T(5), "c4 high y");
    bool passedc45 = f(c, f_l_z, Cell_T(1), "c4 low z");
    bool passedc46 = f(c, f_h_z, null_cell, "c4 high z");
    EXPECT_TRUE(passedc41);
    EXPECT_TRUE(passedc42);
    EXPECT_TRUE(passedc43);
    EXPECT_TRUE(passedc44);
    EXPECT_TRUE(passedc45);
    EXPECT_TRUE(passedc46);
  }
  {
    Cell_T c(5);
    bool passedc51 = f(c, f_l_x, null_cell, "c5 low x");
    bool passedc52 = f(c, f_h_x, null_cell, "c5 high x");
    bool passedc53 = f(c, f_l_y, Cell_T(4), "c5 low y");
    bool passedc54 = f(c, f_h_y, null_cell, "c5 high y");
    bool passedc55 = f(c, f_l_z, Cell_T{2}, "c5 low z");
    bool passedc56 = f(c, f_h_z, null_cell, "c5 high z");
    EXPECT_TRUE(passedc51);
    EXPECT_TRUE(passedc52);
    EXPECT_TRUE(passedc53);
    EXPECT_TRUE(passedc54);
    EXPECT_TRUE(passedc55);
    EXPECT_TRUE(passedc56);
  }
  return;
}  // test_6

TEST(nut_mesh_3D_Cartesian, sample_position_buffer_rng)
{
  using test_aux::soft_expect;
  typedef nut::Buffer_RNG<double> rng_t;
  size_t const szRns(5);
  double rns[szRns] = {0.5, 0.5, 0.5, 0.5, 0.5};
  rng_t rng(rns, szRns);

  Geom_T const dx(2.0);
  Geom_T const dy(3.0);
  Geom_T const dz(5.0);
  Index_T const nx = 1;
  Index_T const ny = 3;
  Index_T const nz = 2;
  Mesh_Iface_T mesh(nx, ny, nz, dx, dy, dz);

  Cell_T c0{0};
  Mesh_Iface_T::Vector x = mesh.sample_position(rng, c0);
  Mesh_Iface_T::Vector o = mesh.sample_direction_isotropic(rng);

  bool passed = expect(x[0], 1.0, "cell 0, x") &&
                expect(x[1], 1.5, "cell 0, y") &&
                expect(x[2], 2.5, "cell 0, z") &&
                soft_expect(o[0], -1.0, "cell 0, omega_x") &&
                soft_expect(o[1], 0.0, "cell 0, omega_y") &&
                soft_expect(o[2], 0.0, "cell 0, omega_z");
  EXPECT_TRUE(passed);
  return;
}  // test_7

// TEST(nut_mesh_3D_Cartesian, distance_to_bdy)
// {
//   using nut::vec_t;
//   typedef nut::Buffer_RNG<double> rng_t;
//   size_t const szRns(5);
//   double rns[szRns] = {0.5, 0.5, 0.5, 0.5, 0.5};
//   rng_t rng(rns, szRns);

//   Geom_T const dx(2.0);
//   Geom_T const dy(3.0);
//   Geom_T const dz(5.0);
//   Index_T const nx = 1;
//   Index_T const ny = 3;
//   Index_T const nz = 2;
//   vbd bds;
//   nut::mkReflectBCs<Mesh_Iface_T, Cell_T>(bds, nx, ny, nz);
//   Mesh_Iface_T mesh(nx, ny, nz, dx, dy, dz, bds);

//   vec_t<3> x, o;
//   x.v = {{0.8, 2.9, 4.9}};
//   o.v = {{1.0, 0.0, 0.0}};
//   Mesh_Iface_T::coord_t crd = {x, o};
//   Mesh_Iface_T::d_to_b_t d2b = mesh.distance_to_bdy(crd, Cell_T(0));

//   bool passed = expect(d2b.d, 1.2, "distance") &&
//                 expect(d2b.face, f_h_x, "face");
//   EXPECT_TRUE(passed);
//   return;
// }  // test_8

#endif  // HAVE_MURMELN

// End of file
