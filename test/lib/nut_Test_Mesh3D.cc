// T. M. Kelley (c) 2011 LANS LLC

#include "Mesh3DCar.hh"
#include "expect.hh"
#include "gtest/gtest.h"
#include "types.hh"

using nut::Equal;
using nut::geom_t;
typedef uint64_t cell_t;
typedef nut::Cartesian_3D<cell_t, double> mesh_t;
typedef mesh_t::vec_bdy_descs vbd;
using test_aux::expect;
using test_aux::soft_expect;

using Face = mesh_t::Face;
Face const f_l_x{mesh_t::low_x};
Face const f_l_y{mesh_t::low_y};
Face const f_l_z{mesh_t::low_z};
Face const f_h_x{mesh_t::high_x};
Face const f_h_y{mesh_t::high_y};
Face const f_h_z{mesh_t::high_z};

TEST(nut_mesh_3D, instantiate)
{
  cell_t const nx = 1;
  cell_t const ny = 1;
  cell_t const nz = 1;
  vbd bds;
  nut::mkReflectBCs<mesh_t, cell_t>(bds, nx, ny, nz);
  mesh_t mesh(nx, ny, nz, 1.0, 1.0, 1.0, bds);
  EXPECT_TRUE(true);
  return;
}

template <typename cell_t>
bool
test_2_core(cell_t const nx, cell_t const ny, cell_t const nz)
{
  vbd bds;
  nut::mkReflectBCs<mesh_t, cell_t>(bds, nx, ny, nz);
  mesh_t mesh(nx, ny, nz, 1.0, 1.0, 1.0, bds);
  cell_t const n_cs(nx * ny * nz);
  bool passed = n_cs == mesh.n_cells();
  return passed;
}

TEST(nut_mesh_3D, instantiate_correctly_n_cells)
{
  cell_t const nx = 3;
  cell_t const ny = 5;
  cell_t const nz = 7;
  EXPECT_TRUE(test_2_core(nx, ny, nz));
  return;
}  // test_2

// for any 3D cartesian mesh, the volume should be constant
// (and simple!)
template <typename cell_t>
bool
test_3_core(geom_t const dx,
            cell_t const nx,
            geom_t const dy,
            cell_t const ny,
            geom_t const dz,
            cell_t const nz)
{
  bool passed(true);
  vbd bds;
  nut::mkReflectBCs<mesh_t, cell_t>(bds, nx, ny, nz);
  mesh_t mesh(nx, ny, nz, dx, dy, dz, bds);
  geom_t const exp_vol(dx * dy * dz);
  for(cell_t c = 0; c < mesh.n_cells(); ++c) {
    geom_t const volume = mesh.volume(c);
    bool cok = volume == exp_vol;
    passed = passed && cok;
  }
  return passed;
}

TEST(nut_mesh_3D, instantiate_correctly_volume)
{
  geom_t const dx(2.0);
  geom_t const dy(3.0);
  geom_t const dz(5.0);
  cell_t const nx = 3;
  cell_t const ny = 5;
  cell_t const nz = 7;

  bool passed = test_3_core(dx, nx, dy, ny, dz, nz);
  EXPECT_TRUE(passed);
  return;
}  // test_3

TEST(nut_mesh_3D, isBoundary_3_cell_line)
{
  using test_aux::expect;
  geom_t const dx(2.0);
  geom_t const dy(3.0);
  geom_t const dz(5.0);

  cell_t const nx = 3;
  cell_t const ny = 1;
  cell_t const nz = 1;

  vbd bds;
  nut::mkReflectBCs<mesh_t, cell_t>(bds, nx, ny, nz);
  mesh_t mesh(nx, ny, nz, dx, dy, dz, bds);

  {
    cell_t c(0);
    bool passedc0 =
        expect(mesh.isBoundary(c, mesh_t::low_x), true, "c0 low_x") &&
        expect(mesh.isBoundary(c, mesh_t::high_x), false, "c0 high_x") &&
        expect(mesh.isBoundary(c, mesh_t::low_y), true, "c0 low_y") &&
        expect(mesh.isBoundary(c, mesh_t::high_y), true, "c0 high_y") &&
        expect(mesh.isBoundary(c, mesh_t::low_z), true, "c0 low_z") &&
        expect(mesh.isBoundary(c, mesh_t::high_z), true, "c0 high_z");
    EXPECT_TRUE(passedc0);
  }
  {
    cell_t c(1);
    bool passedc1 =
        expect(mesh.isBoundary(c, mesh_t::low_x), false, "c1 low_x") &&
        expect(mesh.isBoundary(c, mesh_t::high_x), false, "c1 high_x") &&
        expect(mesh.isBoundary(c, mesh_t::low_y), true, "c1 low_y") &&
        expect(mesh.isBoundary(c, mesh_t::high_y), true, "c1 high_y") &&
        expect(mesh.isBoundary(c, mesh_t::low_z), true, "c1 low_z") &&
        expect(mesh.isBoundary(c, mesh_t::high_z), true, "c1 high_z");
    EXPECT_TRUE(passedc1);
  }
  {
    cell_t c(2);
    bool passedc2 =
        expect(mesh.isBoundary(c, mesh_t::low_x), false, "c2 low_x") &&
        expect(mesh.isBoundary(c, mesh_t::high_x), true, "c2 high_x") &&
        expect(mesh.isBoundary(c, mesh_t::low_y), true, "c2 low_y") &&
        expect(mesh.isBoundary(c, mesh_t::high_y), true, "c2 high_y") &&
        expect(mesh.isBoundary(c, mesh_t::low_z), true, "c2 low_z") &&
        expect(mesh.isBoundary(c, mesh_t::high_z), true, "c2 high_z");
    EXPECT_TRUE(passedc2);
  }
  return;
}  // test_4

TEST(nut_mesh_3D, isBoundary_6_cell_cluster)
{
  using nut::ijk_t;
  using test_aux::expect;
  geom_t const dx(2.0);
  geom_t const dy(3.0);
  geom_t const dz(5.0);

  cell_t const nx = 1;
  cell_t const ny = 3;
  cell_t const nz = 2;

  vbd bds;
  nut::mkReflectBCs<mesh_t, cell_t>(bds, nx, ny, nz);
  mesh_t mesh(nx, ny, nz, dx, dy, dz, bds);

  {
    cell_t c(0);  // (0,0,0)
    bool passedc0 =
        expect(mesh.isBoundary(c, mesh_t::low_x), true, "c0 low_x") &&
        expect(mesh.isBoundary(c, mesh_t::high_x), true, "c0 high_x") &&
        expect(mesh.isBoundary(c, mesh_t::low_y), true, "c0 low_y") &&
        expect(mesh.isBoundary(c, mesh_t::high_y), false, "c0 high_y") &&
        expect(mesh.isBoundary(c, mesh_t::low_z), true, "c0 low_z") &&
        expect(mesh.isBoundary(c, mesh_t::high_z), false, "c0 high_z");
    EXPECT_TRUE(passedc0);
  }
  {
    cell_t c(1);  // (0,1,0)
    bool passedc1 =
        expect(mesh.isBoundary(c, mesh_t::low_x), true, "c1 low_x") &&
        expect(mesh.isBoundary(c, mesh_t::high_x), true, "c1 high_x") &&
        expect(mesh.isBoundary(c, mesh_t::low_y), false, "c1 low_y") &&
        expect(mesh.isBoundary(c, mesh_t::high_y), false, "c1 high_y") &&
        expect(mesh.isBoundary(c, mesh_t::low_z), true, "c1 low_z") &&
        expect(mesh.isBoundary(c, mesh_t::high_z), false, "c1 high_z");
    EXPECT_TRUE(passedc1);
  }
  {
    cell_t c(2);  // (0,2,0)
    bool passedc2 =
        expect(mesh.isBoundary(c, mesh_t::low_x), true, "c2 low_x") &&
        expect(mesh.isBoundary(c, mesh_t::high_x), true, "c2 high_x") &&
        expect(mesh.isBoundary(c, mesh_t::low_y), false, "c2 low_y") &&
        expect(mesh.isBoundary(c, mesh_t::high_y), true, "c2 high_y") &&
        expect(mesh.isBoundary(c, mesh_t::low_z), true, "c2 low_z") &&
        expect(mesh.isBoundary(c, mesh_t::high_z), false, "c2 high_z");
    EXPECT_TRUE(passedc2);
  }
  {
    cell_t c(3);  // (0,0,1)
    bool passedc3 =
        expect(mesh.isBoundary(c, mesh_t::low_x), true, "c3 low_x") &&
        expect(mesh.isBoundary(c, mesh_t::high_x), true, "c3 high_x") &&
        expect(mesh.isBoundary(c, mesh_t::low_y), true, "c3 low_y") &&
        expect(mesh.isBoundary(c, mesh_t::high_y), false, "c3 high_y") &&
        expect(mesh.isBoundary(c, mesh_t::low_z), false, "c3 low_z") &&
        expect(mesh.isBoundary(c, mesh_t::high_z), true, "c3 high_z");
    EXPECT_TRUE(passedc3);
  }
  {
    cell_t c(4);  // (0,1,1)
    bool passedc4 =
        expect(mesh.isBoundary(c, mesh_t::low_x), true, "c4 low_x") &&
        expect(mesh.isBoundary(c, mesh_t::high_x), true, "c4 high_x") &&
        expect(mesh.isBoundary(c, mesh_t::low_y), false, "c4 low_y") &&
        expect(mesh.isBoundary(c, mesh_t::high_y), false, "c4 high_y") &&
        expect(mesh.isBoundary(c, mesh_t::low_z), false, "c4 low_z") &&
        expect(mesh.isBoundary(c, mesh_t::high_z), true, "c4 high_z");
    EXPECT_TRUE(passedc4);
  }
  {
    cell_t c(5);  // (0,2,1)
    bool passedc5 =
        expect(mesh.isBoundary(c, mesh_t::low_x), true, "c5 low_x") &&
        expect(mesh.isBoundary(c, mesh_t::high_x), true, "c5 high_x") &&
        expect(mesh.isBoundary(c, mesh_t::low_y), false, "c5 low_y") &&
        expect(mesh.isBoundary(c, mesh_t::high_y), true, "c5 high_y") &&
        expect(mesh.isBoundary(c, mesh_t::low_z), false, "c5 low_z") &&
        expect(mesh.isBoundary(c, mesh_t::high_z), true, "c5 high_z");
    EXPECT_TRUE(passedc5);
  }
  return;
}  // test_5

TEST(nut_mesh_3D, cell_across_6_cell_cluster)
{
  using nut::ijk_t;
  geom_t const dx(2.0);
  geom_t const dy(3.0);
  geom_t const dz(5.0);

  cell_t const nx = 1;
  cell_t const ny = 3;
  cell_t const nz = 2;

  vbd bds;
  nut::mkReflectBCs<mesh_t, cell_t>(bds, nx, ny, nz);
  mesh_t mesh(nx, ny, nz, dx, dy, dz, bds);
  mesh_t::Vector unused{0.0};
  auto f = [&mesh, &unused](cell_t c, mesh_t::face_t const f, cell_t const e,
                            const char * s) {
    return expect(mesh.cell_across(c, f, unused), e, s);
  };
  cell_t const null_cell = mesh_t::null_cell();
  {
    cell_t c(0);  // (0,0,0)
    bool passedc01 = f(c, f_l_x, null_cell, "c0 low x");
    bool passedc02 = f(c, f_h_x, null_cell, "c0 high x");
    bool passedc03 = f(c, f_l_y, null_cell, "c0 low y");
    bool passedc04 = f(c, f_h_y, cell_t(1), "c0 high y");
    bool passedc05 = f(c, f_l_z, null_cell, "c0 low z");
    bool passedc06 = f(c, f_h_z, cell_t(3), "c0 high z");
    EXPECT_TRUE(passedc01);
    EXPECT_TRUE(passedc02);
    EXPECT_TRUE(passedc03);
    EXPECT_TRUE(passedc04);
    EXPECT_TRUE(passedc05);
    EXPECT_TRUE(passedc06);
  }
  {
    cell_t c(1);
    bool passedc11 = f(c, f_l_x, null_cell, "c1 low x");
    bool passedc12 = f(c, f_h_x, null_cell, "c1 high x");
    bool passedc13 = f(c, f_l_y, cell_t(0), "c1 low y");
    bool passedc14 = f(c, f_h_y, cell_t(2), "c1 high y");
    bool passedc15 = f(c, f_l_z, null_cell, "c1 low z");
    bool passedc16 = f(c, f_h_z, cell_t(4), "c1 high z");
    EXPECT_TRUE(passedc11);
    EXPECT_TRUE(passedc12);
    EXPECT_TRUE(passedc13);
    EXPECT_TRUE(passedc14);
    EXPECT_TRUE(passedc15);
    EXPECT_TRUE(passedc16);
  }
  {
    cell_t c(2);
    bool passedc21 = f(c, f_l_x, null_cell, "c2 low x");
    bool passedc22 = f(c, f_h_x, null_cell, "c2 high x");
    bool passedc23 = f(c, f_l_y, cell_t(1), "c2 low y");
    bool passedc24 = f(c, f_h_y, null_cell, "c2 high y");
    bool passedc25 = f(c, f_l_z, null_cell, "c2 low z");
    bool passedc26 = f(c, f_h_z, cell_t(5), "c2 high z");
    EXPECT_TRUE(passedc21);
    EXPECT_TRUE(passedc22);
    EXPECT_TRUE(passedc23);
    EXPECT_TRUE(passedc24);
    EXPECT_TRUE(passedc25);
    EXPECT_TRUE(passedc26);
  }
  {
    cell_t c(3);
    bool passedc31 = f(c, f_l_x, null_cell, "c3 low x");
    bool passedc32 = f(c, f_h_x, null_cell, "c3 high x");
    bool passedc33 = f(c, f_l_y, null_cell, "c3 low y");
    bool passedc34 = f(c, f_h_y, cell_t(4), "c3 high y");
    bool passedc35 = f(c, f_l_z, cell_t(0), "c3 low z");
    bool passedc36 = f(c, f_h_z, null_cell, "c3 high z");
    EXPECT_TRUE(passedc31);
    EXPECT_TRUE(passedc32);
    EXPECT_TRUE(passedc33);
    EXPECT_TRUE(passedc34);
    EXPECT_TRUE(passedc35);
    EXPECT_TRUE(passedc36);
  }
  {
    cell_t c(4);
    bool passedc41 = f(c, f_l_x, null_cell, "c4 low x");
    bool passedc42 = f(c, f_h_x, null_cell, "c4 high x");
    bool passedc43 = f(c, f_l_y, cell_t(3), "c4 low y");
    bool passedc44 = f(c, f_h_y, cell_t(5), "c4 high y");
    bool passedc45 = f(c, f_l_z, cell_t(1), "c4 low z");
    bool passedc46 = f(c, f_h_z, null_cell, "c4 high z");
    EXPECT_TRUE(passedc41);
    EXPECT_TRUE(passedc42);
    EXPECT_TRUE(passedc43);
    EXPECT_TRUE(passedc44);
    EXPECT_TRUE(passedc45);
    EXPECT_TRUE(passedc46);
  }
  {
    cell_t c(5);
    bool passedc51 = f(c, f_l_x, null_cell, "c5 low x");
    bool passedc52 = f(c, f_h_x, null_cell, "c5 high x");
    bool passedc53 = f(c, f_l_y, cell_t(4), "c5 low y");
    bool passedc54 = f(c, f_h_y, null_cell, "c5 high y");
    bool passedc55 = f(c, f_l_z, cell_t{2}, "c5 low z");
    bool passedc56 = f(c, f_h_z, null_cell, "c5 high z");
    EXPECT_TRUE(passedc51);
    EXPECT_TRUE(passedc52);
    EXPECT_TRUE(passedc53);
    EXPECT_TRUE(passedc54);
    EXPECT_TRUE(passedc55);
    EXPECT_TRUE(passedc56);
  }

  // {
  //   cell_t c(4);
  //   bool passedc4 = f(c, mesh_t::low_x, cell_t(4), "c4 low x") &&
  //                   f(c, mesh_t::high_x, cell_t(4), "c4 high x") &&
  //                   f(c, mesh_t::low_y, cell_t(3), "c4 low y") &&
  //                   f(c, mesh_t::high_y, cell_t(5), "c4 high y") &&
  //                   f(c, mesh_t::low_z, cell_t(1), "c4 low z") &&
  //                   f(c, mesh_t::high_z, cell_t(4), "c4 high z");
  //   EXPECT_TRUE(passedc4);
  // }
  // {
  //   cell_t c(5);
  //   bool passedc5 = f(c, mesh_t::low_x, cell_t(5), "c5 low x") &&
  //                   f(c, mesh_t::high_x, cell_t(5), "c5 high x") &&
  //                   f(c, mesh_t::low_y, cell_t(4), "c5 low y") &&
  //                   f(c, mesh_t::high_y, cell_t(5), "c5 high y") &&
  //                   f(c, mesh_t::low_z, cell_t(2), "c5 low z") &&
  //                   f(c, mesh_t::high_z, cell_t(5), "c5 high z");
  //   EXPECT_TRUE(passedc5);
  // }
  return;
}  // test_6

TEST(nut_mesh_3D, sample_position_buffer_rng)
{
  using nut::vec_t;
  typedef nut::Buffer_RNG<double> rng_t;
  size_t const szRns(5);
  double rns[szRns] = {0.5, 0.5, 0.5, 0.5, 0.5};
  rng_t rng(rns, szRns);

  geom_t const dx(2.0);
  geom_t const dy(3.0);
  geom_t const dz(5.0);
  cell_t const nx = 1;
  cell_t const ny = 3;
  cell_t const nz = 2;
  vbd bds;
  nut::mkReflectBCs<mesh_t, cell_t>(bds, nx, ny, nz);
  mesh_t mesh(nx, ny, nz, dx, dy, dz, bds);

  mesh_t::Vector x = mesh.sample_position(rng, 0);
  mesh_t::Vector o = mesh.sample_direction_isotropic(rng);

  bool passed = expect(x.v[0], 1.0, "cell 0, x") &&
                expect(x.v[1], 1.5, "cell 0, y") &&
                expect(x.v[2], 2.5, "cell 0, z") &&
                soft_expect(o.v[0], -1.0, "cell 0, omega_x") &&
                soft_expect(o.v[1], 0.0, "cell 0, omega_y") &&
                soft_expect(o.v[2], 0.0, "cell 0, omega_z");
  EXPECT_TRUE(passed);
  return;
}  // test_7

TEST(nut_mesh_3D, distance_to_bdy)
{
  using nut::vec_t;
  typedef nut::Buffer_RNG<double> rng_t;
  size_t const szRns(5);
  double rns[szRns] = {0.5, 0.5, 0.5, 0.5, 0.5};
  rng_t rng(rns, szRns);

  geom_t const dx(2.0);
  geom_t const dy(3.0);
  geom_t const dz(5.0);
  cell_t const nx = 1;
  cell_t const ny = 3;
  cell_t const nz = 2;
  vbd bds;
  nut::mkReflectBCs<mesh_t, cell_t>(bds, nx, ny, nz);
  mesh_t mesh(nx, ny, nz, dx, dy, dz, bds);

  vec_t<3> x, o;
  x.v = {{0.8, 2.9, 4.9}};
  o.v = {{1.0, 0.0, 0.0}};
  mesh_t::coord_t crd = {x, o};
  mesh_t::d_to_b_t d2b = mesh.distance_to_bdy(crd, cell_t(0));

  bool passed = expect(d2b.d, 1.2, "distance") &&
                expect(d2b.face, mesh_t::high_x, "face");
  EXPECT_TRUE(passed);
  return;
}  // test_8

// End of file
