// test_Cartesian_Mesh.cc
// Apr 04, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#include "base/soft_equiv.h"
#include "base/test_common.h"
#include "cartesian/Cartesian_Mesh.h"
#include "gtest/gtest.h"
#include <array>

using nut_mesh::geom_t;
using nut_mesh::soft_equiv_os;
using namespace nut_mesh;
using cell_handle_t = Cartesian_Mesh::cell_handle_t;
using Face_Name = Cartesian_Mesh::Face_Name;
using Vector = Cartesian_Mesh::Vector;
using Point = Cartesian_Mesh::Point;
using geometron = Cartesian_Mesh::Geom_State;

nut_mesh::Cartesian_Face const f_l_x{Cartesian_Mesh::LOW_X};
nut_mesh::Cartesian_Face const f_l_y{Cartesian_Mesh::LOW_Y};
nut_mesh::Cartesian_Face const f_l_z{Cartesian_Mesh::LOW_Z};
nut_mesh::Cartesian_Face const f_h_x{Cartesian_Mesh::HIGH_X};
nut_mesh::Cartesian_Face const f_h_y{Cartesian_Mesh::HIGH_Y};
nut_mesh::Cartesian_Face const f_h_z{Cartesian_Mesh::HIGH_Z};

Face_Name const fn_l_x{Cartesian_Mesh::LOW_X};
Face_Name const fn_l_y{Cartesian_Mesh::LOW_Y};
Face_Name const fn_l_z{Cartesian_Mesh::LOW_Z};
Face_Name const fn_h_x{Cartesian_Mesh::HIGH_X};
Face_Name const fn_h_y{Cartesian_Mesh::HIGH_Y};
Face_Name const fn_h_z{Cartesian_Mesh::HIGH_Z};

// A 1x1x1 Cartesian_Mesh
Cartesian_Mesh one_cell_cmi() {
  index_t const nx = 1, ny = 1, nz = 1;
  geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
  geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
  Cartesian_Mesh m(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
  return m;
}

// A 2x2x2 Cartesian_Mesh
Cartesian_Mesh eight_cell_cmi() {
  index_t const nx = 2, ny = 2, nz = 2;
  geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
  geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
  Cartesian_Mesh m(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
  return m;
}

// A 3x5x7 Cartesian_Mesh
namespace {
Cartesian_Mesh the_other_cmi() {
  index_t const nx = 3, ny = 5, nz = 7;
  geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
  geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
  Cartesian_Mesh m(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
  return m;
}
} // namespace

TEST(default_mesh_Cartesian_mesh, boundary_face_iterators) {
  using i_t = Cartesian_Mesh::boundary_face_iterator;
  Cartesian_Mesh m{one_cell_cmi()};
  {
    i_t i_l_x_b = m.boundary_faces_begin(Cartesian_Mesh::LOW_X);
    i_t const i_l_x_e = m.boundary_faces_end(Cartesian_Mesh::LOW_X);
    EXPECT_EQ(i_l_x_b.ix, 0u);
    EXPECT_EQ(i_l_x_b.xmax, 1u);
    EXPECT_EQ(i_l_x_b.iy, 0u);
    EXPECT_EQ(i_l_x_b.ymax, 1u);
    EXPECT_EQ(i_l_x_b.iz, 0u);
    EXPECT_EQ(i_l_x_b.zmax, 1u);
    auto big_diff{std::distance(i_l_x_b, i_l_x_e)};
    EXPECT_EQ(big_diff, 1);
    EXPECT_TRUE(i_l_x_b); // iterator converts to true when not at end
    for (; i_l_x_b != i_l_x_e; ++i_l_x_b) {
      Cartesian_Face f{*i_l_x_b};
      EXPECT_EQ(f.as_id(), 0u); // there's only one, so...
    }
    EXPECT_FALSE(i_l_x_b); // iterator converts to false when at end
  }
  {
    // low-y
    i_t i_l_y_b = m.boundary_faces_begin(Cartesian_Mesh::LOW_Y);
    i_t const i_l_y_e = m.boundary_faces_end(Cartesian_Mesh::LOW_Y);
    auto big_diff{std::distance(i_l_y_b, i_l_y_e)};
    EXPECT_EQ(big_diff, 1);
    EXPECT_TRUE(i_l_y_b); // iterator converts to true when not at end
    for (; i_l_y_b != i_l_y_e; ++i_l_y_b) {
      Cartesian_Face f{*i_l_y_b};
      EXPECT_EQ(f.as_id(), 2u); // there's only one, so...
    }
    EXPECT_FALSE(i_l_y_b); // iterator converts to false when at end
  }
  {
    // low-z
    i_t i_l_z_b = m.boundary_faces_begin(Cartesian_Mesh::LOW_Z);
    i_t const i_l_z_e = m.boundary_faces_end(Cartesian_Mesh::LOW_Z);
    auto big_diff{std::distance(i_l_z_b, i_l_z_e)};
    EXPECT_EQ(big_diff, 1);
    EXPECT_TRUE(i_l_z_b); // iterator converts to true when not at end
    for (; i_l_z_b != i_l_z_e; ++i_l_z_b) {
      Cartesian_Face f{*i_l_z_b};
      EXPECT_EQ(f.as_id(), 4u); // there's only one, so...
    }
    EXPECT_FALSE(i_l_z_b); // iterator converts to false when at end
  }
  {
    // high-x
    i_t i_h_x_b = m.boundary_faces_begin(Cartesian_Mesh::HIGH_X);
    i_t const i_h_x_e = m.boundary_faces_end(Cartesian_Mesh::HIGH_X);
    auto big_diff{std::distance(i_h_x_b, i_h_x_e)};
    EXPECT_EQ(big_diff, 1);
    EXPECT_TRUE(i_h_x_b); // iterator converts to true when not at end
    for (; i_h_x_b != i_h_x_e; ++i_h_x_b) {
      Cartesian_Face f{*i_h_x_b};
      EXPECT_EQ(f.as_id(), 1u); // there's only one, so...
    }
    EXPECT_FALSE(i_h_x_b); // iterator converts to false when at end
  }
  {
    // high-y
    i_t i_h_y_b = m.boundary_faces_begin(Cartesian_Mesh::HIGH_Y);
    i_t const i_h_y_e = m.boundary_faces_end(Cartesian_Mesh::HIGH_Y);
    auto big_diff{std::distance(i_h_y_b, i_h_y_e)};
    EXPECT_EQ(big_diff, 1);
    EXPECT_TRUE(i_h_y_b); // iterator converts to true when not at end
    for (; i_h_y_b != i_h_y_e; ++i_h_y_b) {
      Cartesian_Face f{*i_h_y_b};
      EXPECT_EQ(f.as_id(), 3u); // there's only one, so...
    }
    EXPECT_FALSE(i_h_y_b); // iterator converts to false when at end
  }
  {
    // high-z
    i_t i_h_z_b = m.boundary_faces_begin(Cartesian_Mesh::HIGH_Z);
    i_t const i_h_z_e = m.boundary_faces_end(Cartesian_Mesh::HIGH_Z);
    auto big_diff{std::distance(i_h_z_b, i_h_z_e)};
    EXPECT_EQ(big_diff, 1);
    EXPECT_TRUE(i_h_z_b); // iterator converts to true when not at end
    for (; i_h_z_b != i_h_z_e; ++i_h_z_b) {
      Cartesian_Face f{*i_h_z_b};
      EXPECT_EQ(f.as_id(), 5u); // there's only one, so...
    }
    EXPECT_FALSE(i_h_z_b); // iterator converts to false when at end
  }
  return;
} // TEST(default_mesh_Cartesian_mesh, boundary_face_iterators){

TEST(default_mesh_Cartesian_mesh, boundary_faces) {
  using r_t = std::vector<Cartesian_Mesh::face_handle_t>;
  // one cell mesh
  {
    Cartesian_Mesh m{one_cell_cmi()};
    auto v_l_x = m.boundary_faces(Cartesian_Mesh::LOW_X);
    EXPECT_EQ(v_l_x.size(), 1u);
    EXPECT_EQ(v_l_x[0].as_id(), 0u);
    auto v_h_x = m.boundary_faces(Cartesian_Mesh::HIGH_X);
    EXPECT_EQ(v_h_x.size(), 1u);
    EXPECT_EQ(v_h_x[0].as_id(), 1u);
    auto v_l_y = m.boundary_faces(Cartesian_Mesh::LOW_Y);
    EXPECT_EQ(v_l_y.size(), 1u);
    EXPECT_EQ(v_l_y[0].as_id(), 2u);
    auto v_h_y = m.boundary_faces(Cartesian_Mesh::HIGH_Y);
    EXPECT_EQ(v_h_y.size(), 1u);
    EXPECT_EQ(v_h_y[0].as_id(), 3u);
    auto v_l_z = m.boundary_faces(Cartesian_Mesh::LOW_Z);
    EXPECT_EQ(v_l_z.size(), 1u);
    EXPECT_EQ(v_l_z[0].as_id(), 4u);
    auto v_h_z = m.boundary_faces(Cartesian_Mesh::HIGH_Z);
    EXPECT_EQ(v_h_z.size(), 1u);
    EXPECT_EQ(v_h_z[0].as_id(), 5u);
  }
  // 2 x 2 x 2 mesh
  {
    Cartesian_Mesh m{eight_cell_cmi()};
    r_t v_l_x = m.boundary_faces(Cartesian_Mesh::LOW_X);
    r_t v_l_x_exp{0u, 3u, 6u, 9u};
    EXPECT_EQ(v_l_x.size(), v_l_x_exp.size());
    for (size_t i = 0; i < v_l_x_exp.size(); ++i) {
      EXPECT_EQ(v_l_x[i], v_l_x_exp[i]);
    }
    r_t v_l_y = m.boundary_faces(Cartesian_Mesh::LOW_Y);
    r_t v_l_y_exp{12u, 13u, 18u, 19u};
    EXPECT_EQ(v_l_y.size(), v_l_y_exp.size());
    for (size_t i = 0; i < v_l_y_exp.size(); ++i) {
      EXPECT_EQ(v_l_y[i], v_l_y_exp[i]);
    }
    r_t v_l_z = m.boundary_faces(Cartesian_Mesh::LOW_Z);
    r_t v_l_z_exp{24u, 25u, 26u, 27u};
    EXPECT_EQ(v_l_z.size(), v_l_z_exp.size());
    for (size_t i = 0; i < v_l_z_exp.size(); ++i) {
      EXPECT_EQ(v_l_z[i].as_id(), v_l_z_exp[i].as_id());
    }
    r_t v_h_x = m.boundary_faces(Cartesian_Mesh::HIGH_X);
    r_t v_h_x_exp{2u, 5u, 8u, 11u};
    EXPECT_EQ(v_h_x.size(), v_h_x_exp.size());
    for (size_t i = 0; i < v_h_x_exp.size(); ++i) {
      EXPECT_EQ(v_h_x[i], v_h_x_exp[i]);
    }
    r_t v_h_y = m.boundary_faces(Cartesian_Mesh::HIGH_Y);
    r_t v_h_y_exp{16u, 17u, 22u, 23u};
    EXPECT_EQ(v_h_y.size(), v_h_y_exp.size());
    for (size_t i = 0; i < v_h_y_exp.size(); ++i) {
      EXPECT_EQ(v_h_y[i], v_h_y_exp[i]);
    }
    r_t v_h_z = m.boundary_faces(Cartesian_Mesh::HIGH_Z);
    r_t v_h_z_exp{32u, 33u, 34u, 35u};
    EXPECT_EQ(v_h_z.size(), v_h_z_exp.size());
    for (size_t i = 0; i < v_h_z_exp.size(); ++i) {
      EXPECT_EQ(v_h_z[i], v_h_z_exp[i]);
    }
  }
  // EXPECT_TRUE(false);  // Failed 11/22/2019
  return;
} // TEST(default_mesh_Cartesian_mesh, boundary_faces){

TEST(default_mesh_Cartesian_mesh, linear_to_cartesian) {
  using nut_mesh::linear_to_cartesian;
  {
    index_t nx = 1;
    index_t ny = 1;
    index_t c0 = 0;
    auto [ix, iy, iz] = linear_to_cartesian(c0, nx, ny);
    EXPECT_EQ(ix, 0);
    EXPECT_EQ(iy, 0);
    EXPECT_EQ(iz, 0);
  }
  {
    index_t nx = 2;
    index_t ny = 2;
    {
      index_t c0 = 0;
      auto [ix, iy, iz] = linear_to_cartesian(c0, nx, ny);
      EXPECT_EQ(ix, 0);
      EXPECT_EQ(iy, 0);
      EXPECT_EQ(iz, 0);
    }
    {
      index_t c1 = 1;
      auto [ix, iy, iz] = linear_to_cartesian(c1, nx, ny);
      EXPECT_EQ(ix, 1);
      EXPECT_EQ(iy, 0);
      EXPECT_EQ(iz, 0);
    }
    {
      index_t c7 = 7;
      auto [ix, iy, iz] = linear_to_cartesian(c7, nx, ny);
      EXPECT_EQ(ix, 1);
      EXPECT_EQ(iy, 1);
      EXPECT_EQ(iz, 1);
    }
  }
  return;
}

TEST(default_mesh_Cartesian_mesh, cell_to_face) {
  {
    Cartesian_Mesh m{the_other_cmi()};
    cell_handle_t c0{0};
    /* First, test face indices of lowest X, Y, and Z faces*/
    {
      Face_Name lx{Cartesian_Mesh::LOW_X};
      Cartesian_Face cf = m.cell_to_face(c0, lx);
      EXPECT_EQ(0, cf.as_id());
    }
    // There are 140 faces perpendicular to X: first Y face at 140
    {
      Face_Name ly{Cartesian_Mesh::LOW_Y};
      Cartesian_Face cf = m.cell_to_face(c0, ly);
      EXPECT_EQ(140, cf.as_id());
    }
    // There are 126 faces perpendicular to Y, so first Z face is at 140+126
    {
      Face_Name lz{Cartesian_Mesh::LOW_Z};
      Cartesian_Face cf = m.cell_to_face(c0, lz);
      EXPECT_EQ(266, cf.as_id());
    }
    // Next, confirm that all X faces have index lower than any Y, Z face,
    // etc.
    index_t max_x{0}, max_y{0};
    index_t min_y{100000}, min_z{100000};
    for (index_t i = 0; i < m.num_cells(); ++i) {
      cell_handle_t c{i};
      Cartesian_Face flx = m.cell_to_face(c, fn_l_x);
      Cartesian_Face fly = m.cell_to_face(c, fn_l_y);
      Cartesian_Face flz = m.cell_to_face(c, fn_l_z);
      Cartesian_Face fhx = m.cell_to_face(c, fn_h_x);
      Cartesian_Face fhy = m.cell_to_face(c, fn_h_y);
      Cartesian_Face fhz = m.cell_to_face(c, fn_h_z);
      if (flx.as_id() > max_x)
        max_x = flx.as_id();
      if (fhx.as_id() > max_x)
        max_x = fhx.as_id();
      if (fly.as_id() > max_y)
        max_y = fly.as_id();
      if (fhy.as_id() > max_y)
        max_y = fhy.as_id();
      if (fly.as_id() < min_y)
        min_y = fly.as_id();
      if (fhy.as_id() < min_y)
        min_y = fhy.as_id();
      if (flz.as_id() < min_z)
        min_z = flz.as_id();
      if (fhz.as_id() < min_z)
        min_z = fhz.as_id();
    }
    EXPECT_TRUE(max_x < min_y);
    EXPECT_TRUE(max_x < min_z);
    EXPECT_TRUE(max_y < min_z);
    EXPECT_TRUE(max_y > min_y);

    /* In this part of the test, check the faces for every cell. Recall that
     * mesh.get_faces returns faces in the order {low-x, low-y, low-z,
     high-x,
     * high-y, high-z}. The low- faces will be mapped as in the cell, the
     high
     * faces as the low face in the cell that's DIM + 1.  */
    index_t n_cells{m.num_cells()};
    index_t const nx = 3, ny = 5;
    for (index_t i = 0; i < n_cells; ++i) {
      cell_handle_t c{i};
      auto faces = m.get_faces(c);
      index_t ix, iy, iz; // Intel doesn't do structured init til I18
      std::tie(ix, iy, iz) = linear_to_cartesian(c.as_id(), nx, ny);
      cell_handle_t plus_x_cell{cartesian_to_linear(ix + 1, iy, iz, nx, ny)};
      cell_handle_t plus_y_cell{cartesian_to_linear(ix, iy + 1, iz, nx, ny)};
      cell_handle_t plus_z_cell{cartesian_to_linear(ix, iy, iz + 1, nx, ny)};
      Cartesian_Face f0 = m.cell_to_face(c, fn_l_x); // low-x
      Cartesian_Face f1 = m.cell_to_face(c, fn_l_y); // low-y
      Cartesian_Face f2 = m.cell_to_face(c, fn_l_z); // low-z
      Cartesian_Face f3 = m.cell_to_face(c, fn_h_x); // high-x
      Cartesian_Face f4 = m.cell_to_face(c, fn_h_y); // high-y
      Cartesian_Face f5 = m.cell_to_face(c, fn_h_z); // high-z
      EXPECT_EQ(f0, faces[0]);
      EXPECT_EQ(f1, faces[1]);
      EXPECT_EQ(f2, faces[2]);
      EXPECT_EQ(f3, faces[3]);
      EXPECT_EQ(f4, faces[4]);
      EXPECT_EQ(f5, faces[5]);
    }
  }
  /* Want to exercise a couple special cases*/
  {
    Cartesian_Mesh m{eight_cell_cmi()};
    // We should be able to enumerate all the faces?
    size_t const nx{2u};
    size_t const ny{2u};
    {
      cell_handle_t c0{0};
      auto faces1 = m.get_faces(c0);
      // name those faces!
      EXPECT_EQ(faces1[0].as_id(), 0u);
      EXPECT_EQ(faces1[1].as_id(), 12u);
      EXPECT_EQ(faces1[2].as_id(), 24u);
      EXPECT_EQ(faces1[3].as_id(), 1u);
      EXPECT_EQ(faces1[4].as_id(), 14u);
      EXPECT_EQ(faces1[5].as_id(), 28u);
      Cartesian_Face f0 = m.cell_to_face(c0, fn_l_x); // low-x
      Cartesian_Face f1 = m.cell_to_face(c0, fn_l_y); // low-x
      Cartesian_Face f2 = m.cell_to_face(c0, fn_l_z); // low-x
      /* the high x face for cell ix is the low-x face for ix + 1, etc. */
      index_t ix, iy, iz; // Intel doesn't do structured init til I18
      std::tie(ix, iy, iz) = linear_to_cartesian(c0.as_id(), nx, ny);
      cell_handle_t plus_x_cell{cartesian_to_linear(ix + 1, iy, iz, nx, ny)};
      cell_handle_t plus_y_cell{cartesian_to_linear(ix, iy + 1, iz, nx, ny)};
      cell_handle_t plus_z_cell{cartesian_to_linear(ix, iy, iz + 1, nx, ny)};
      Cartesian_Face f3 = m.cell_to_face(plus_x_cell, fn_l_x); // low-x
      Cartesian_Face f4 = m.cell_to_face(plus_y_cell, fn_l_y); // low-x
      Cartesian_Face f5 = m.cell_to_face(plus_z_cell, fn_l_z); // low-x
      EXPECT_EQ(f0, faces1[0]);
      EXPECT_EQ(f1, faces1[1]);
      EXPECT_EQ(f2, faces1[2]);
      EXPECT_EQ(f3, faces1[3]);
      EXPECT_EQ(f4, faces1[4]);
      EXPECT_EQ(f5, faces1[5]);
    }
    {
      cell_handle_t c1{1};
      auto faces1 = m.get_faces(c1);
      // name those faces!
      EXPECT_EQ(faces1[0].as_id(), 1u);
      EXPECT_EQ(faces1[1].as_id(), 13u);
      EXPECT_EQ(faces1[2].as_id(), 25u);
      EXPECT_EQ(faces1[3].as_id(), 2u);
      EXPECT_EQ(faces1[4].as_id(), 15u);
      EXPECT_EQ(faces1[5].as_id(), 29u);
      Cartesian_Face f0 = m.cell_to_face(c1, fn_l_x); // low-x
      Cartesian_Face f1 = m.cell_to_face(c1, fn_l_y); // low-y
      Cartesian_Face f2 = m.cell_to_face(c1, fn_l_z); // low-z
      /* the high x face for cell ix is the low-x face for ix + 1, etc. */
      index_t ix, iy, iz; // Intel doesn't do structured init til I18
      std::tie(ix, iy, iz) = linear_to_cartesian(c1.as_id(), nx, ny);
      cell_handle_t plus_x_cell{cartesian_to_linear(ix + 1, iy, iz, nx, ny)};
      cell_handle_t plus_y_cell{cartesian_to_linear(ix, iy + 1, iz, nx, ny)};
      cell_handle_t plus_z_cell{cartesian_to_linear(ix, iy, iz + 1, nx, ny)};
      Cartesian_Face f3 = m.cell_to_face(c1, fn_h_x); // high-x
      Cartesian_Face f4 = m.cell_to_face(c1, fn_h_y); // high-y
      Cartesian_Face f5 = m.cell_to_face(c1, fn_h_z); // high-z
      EXPECT_EQ(f0.as_id(), faces1[0].as_id());
      EXPECT_EQ(f1.as_id(), faces1[1].as_id());
      EXPECT_EQ(f2.as_id(), faces1[2].as_id());
      EXPECT_EQ(f3.as_id(), faces1[3].as_id());
      EXPECT_EQ(f4.as_id(), faces1[4].as_id());
      EXPECT_EQ(f5.as_id(), faces1[5].as_id());
    }
    {
      cell_handle_t c2{2};
      auto faces1 = m.get_faces(c2);
      // name those faces!
      EXPECT_EQ(faces1[0].as_id(), 3u);
      EXPECT_EQ(faces1[1].as_id(), 14u);
      EXPECT_EQ(faces1[2].as_id(), 26u);
      EXPECT_EQ(faces1[3].as_id(), 4u);
      EXPECT_EQ(faces1[4].as_id(), 16u);
      EXPECT_EQ(faces1[5].as_id(), 30u);
      Cartesian_Face f0 = m.cell_to_face(c2, fn_l_x); // low-x
      Cartesian_Face f1 = m.cell_to_face(c2, fn_l_y); // low-y
      Cartesian_Face f2 = m.cell_to_face(c2, fn_l_z); // low-z
      /* the high x face for cell ix is the low-x face for ix + 1, etc. */
      index_t ix, iy, iz; // Intel doesn't do structured init til I18
      std::tie(ix, iy, iz) = linear_to_cartesian(c2.as_id(), nx, ny);
      cell_handle_t plus_x_cell{cartesian_to_linear(ix + 1, iy, iz, nx, ny)};
      cell_handle_t plus_y_cell{cartesian_to_linear(ix, iy + 1, iz, nx, ny)};
      cell_handle_t plus_z_cell{cartesian_to_linear(ix, iy, iz + 1, nx, ny)};
      Cartesian_Face f3 = m.cell_to_face(c2, fn_h_x); // high-x
      Cartesian_Face f4 = m.cell_to_face(c2, fn_h_y); // high-y
      Cartesian_Face f5 = m.cell_to_face(c2, fn_h_z); // high-z
      EXPECT_EQ(f0.as_id(), faces1[0].as_id());
      EXPECT_EQ(f1.as_id(), faces1[1].as_id());
      EXPECT_EQ(f2.as_id(), faces1[2].as_id());
      EXPECT_EQ(f3.as_id(), faces1[3].as_id());
      EXPECT_EQ(f4.as_id(), faces1[4].as_id());
      EXPECT_EQ(f5.as_id(), faces1[5].as_id());
    }
    {
      cell_handle_t c7{7};
      auto faces1 = m.get_faces(c7);
      // name those faces!
      EXPECT_EQ(faces1[0].as_id(), 10u);
      EXPECT_EQ(faces1[1].as_id(), 21u);
      EXPECT_EQ(faces1[2].as_id(), 31u);
      EXPECT_EQ(faces1[3].as_id(), 11u);
      EXPECT_EQ(faces1[4].as_id(), 23u);
      EXPECT_EQ(faces1[5].as_id(), 35u);
      Cartesian_Face f0 = m.cell_to_face(c7, fn_l_x); // low-x
      Cartesian_Face f1 = m.cell_to_face(c7, fn_l_y); // low-y
      Cartesian_Face f2 = m.cell_to_face(c7, fn_l_z); // low-z
      /* the high x face for cell ix is the low-x face for ix + 1, etc. */
      index_t ix, iy, iz; // Intel doesn't do structured init til I18
      std::tie(ix, iy, iz) = linear_to_cartesian(c7.as_id(), nx, ny);
      cell_handle_t plus_x_cell{cartesian_to_linear(ix + 1, iy, iz, nx, ny)};
      cell_handle_t plus_y_cell{cartesian_to_linear(ix, iy + 1, iz, nx, ny)};
      cell_handle_t plus_z_cell{cartesian_to_linear(ix, iy, iz + 1, nx, ny)};
      Cartesian_Face f3 = m.cell_to_face(c7, fn_h_x); // high-x
      Cartesian_Face f4 = m.cell_to_face(c7, fn_h_y); // high-y
      Cartesian_Face f5 = m.cell_to_face(c7, fn_h_z); // high-z
      EXPECT_EQ(f0.as_id(), faces1[0].as_id());
      EXPECT_EQ(f1.as_id(), faces1[1].as_id());
      EXPECT_EQ(f2.as_id(), faces1[2].as_id());
      EXPECT_EQ(f3.as_id(), faces1[3].as_id());
      EXPECT_EQ(f4.as_id(), faces1[4].as_id());
      EXPECT_EQ(f5.as_id(), faces1[5].as_id());
    }
  }
  return;
} // cell_to_face

TEST(default_mesh_Cartesian_mesh, face_to_cell) {
  {
    Cartesian_Mesh m{one_cell_cmi()};

    index_t const nx = 1, ny = 1;
    for(index_t i = 0; i < 1; ++i) {
      cell_handle_t c{i};
      auto faces = m.get_faces(c);
      index_t ix, iy, iz; // Intel doesn't do structured init til I18
      std::tie(ix, iy, iz) = linear_to_cartesian(c.as_id(), nx, ny);
      cell_handle_t c3_exp{cartesian_to_linear(ix + 1, iy, iz, nx, ny)};
      cell_handle_t c4_exp{cartesian_to_linear(ix, iy + 1, iz, nx, ny)};
      cell_handle_t c5_exp{cartesian_to_linear(ix, iy, iz + 1, nx, ny)};

      cell_handle_t c0 = m.face_to_cell(faces[0]);
      // EXPECT_EQ(faces[0].as_id(), 0u);
      EXPECT_EQ(c0.as_id(), c.as_id());
      cell_handle_t c1 = m.face_to_cell(faces[1]);
      EXPECT_EQ(c1.as_id(), c.as_id());
      cell_handle_t c2 = m.face_to_cell(faces[2]);
      EXPECT_EQ(c2.as_id(), c.as_id());
      cell_handle_t c3 = m.face_to_cell(faces[3]);
      EXPECT_EQ(c3.as_id(), c3_exp.as_id());
      cell_handle_t c4 = m.face_to_cell(faces[4]);
      EXPECT_EQ(c4.as_id(), c4_exp.as_id());
      cell_handle_t c5 = m.face_to_cell(faces[5]);
      EXPECT_EQ(c5.as_id(), c5_exp.as_id());
    } // for(index_t i = 0; i < 1; ++i)
  }   // scope: one cell cartesian mesh iface
  {
    Cartesian_Mesh m{eight_cell_cmi()};

    index_t const nx = 2, ny = 2;
    // for (index_t i = 0; i < n_cells; ++i) {
    for (index_t i = 0; i < 1; ++i) {
      cell_handle_t c{i};
      auto faces = m.get_faces(c);
      index_t ix, iy, iz; // Intel doesn't do structured init til I18
      std::tie(ix, iy, iz) = linear_to_cartesian(c.as_id(), nx, ny);
      cell_handle_t c3_exp{cartesian_to_linear(ix + 1, iy, iz, nx, ny)};
      cell_handle_t c4_exp{cartesian_to_linear(ix, iy + 1, iz, nx, ny)};
      cell_handle_t c5_exp{cartesian_to_linear(ix, iy, iz + 1, nx, ny)};

      cell_handle_t c0 = m.face_to_cell(faces[0]);
      // EXPECT_EQ(faces[0].as_id(), 0u);
      EXPECT_EQ(c0.as_id(), c.as_id());
      cell_handle_t c1 = m.face_to_cell(faces[1]);
      EXPECT_EQ(c1.as_id(), c.as_id());
      cell_handle_t c2 = m.face_to_cell(faces[2]);
      EXPECT_EQ(c2.as_id(), c.as_id());
      cell_handle_t c3 = m.face_to_cell(faces[3]);
      EXPECT_EQ(c3.as_id(), c3_exp.as_id());
      cell_handle_t c4 = m.face_to_cell(faces[4]);
      EXPECT_EQ(c4.as_id(), c4_exp.as_id());
      cell_handle_t c5 = m.face_to_cell(faces[5]);
      EXPECT_EQ(c5.as_id(), c5_exp.as_id());
    } // for(index_t i = 0; i < 1; ++i)
  }   // scope: eight cell cartesian mesh iface
  {
    Cartesian_Mesh m{the_other_cmi()};

    index_t const nx = 3, ny = 5;
    // for (index_t i = 0; i < n_cells; ++i) {
    for (index_t i = 0; i < 1; ++i) {
      cell_handle_t c{i};
      auto faces = m.get_faces(c);
      index_t ix, iy, iz; // Intel doesn't do structured init til I18
      std::tie(ix, iy, iz) = linear_to_cartesian(c.as_id(), nx, ny);
      cell_handle_t c3_exp{cartesian_to_linear(ix + 1, iy, iz, nx, ny)};
      cell_handle_t c4_exp{cartesian_to_linear(ix, iy + 1, iz, nx, ny)};
      cell_handle_t c5_exp{cartesian_to_linear(ix, iy, iz + 1, nx, ny)};

      cell_handle_t c0 = m.face_to_cell(faces[0]);
      // EXPECT_EQ(faces[0].as_id(), 0u);
      EXPECT_EQ(c0.as_id(), c.as_id());
      cell_handle_t c1 = m.face_to_cell(faces[1]);
      EXPECT_EQ(c1.as_id(), c.as_id());
      cell_handle_t c2 = m.face_to_cell(faces[2]);
      EXPECT_EQ(c2.as_id(), c.as_id());
      cell_handle_t c3 = m.face_to_cell(faces[3]);
      EXPECT_EQ(c3.as_id(), c3_exp.as_id());
      cell_handle_t c4 = m.face_to_cell(faces[4]);
      EXPECT_EQ(c4.as_id(), c4_exp.as_id());
      cell_handle_t c5 = m.face_to_cell(faces[5]);
      EXPECT_EQ(c5.as_id(), c5_exp.as_id());
    } // for(index_t i = 0; i < 1; ++i)
  }   // scope: 105 cell cartesian mesh iface
  return;
} // face_to_cell

/* This is just a smaller version of the original instantiate test (below)
 * ... less to read in the error reporting :) */
TEST(default_mesh_Cartesian_mesh, instantiate_tiny) {
  index_t const nx = 2, ny = 2, nz = 2;
  index_t const n_yz_faces(num_yz_faces(nx, ny, nz));
  index_t const n_xz_faces(num_xz_faces(nx, ny, nz));
  Cartesian_Mesh m(nx, ny, nz);
  index_t exp_n_cells = nx * ny * nz;
  auto &cells = m.cells();
  EXPECT_EQ(cells.size(), exp_n_cells);

  // offsets in face arrays
  index_t const vertex_offset_x = 1;
  index_t const vertex_offset_y = nx + 1;
  index_t const vertex_offset_z = (nx + 1) * (ny + 1);

  auto valid_cell = [&cells](cell_handle_t const &c) {
    return cells.find(c) != cells.end();
  };

  for (index_t i = 0; i < 1; ++i) {
    // We expect each cell from [0,n_cells) to be in the mesh and no others
    cell_handle_t c(i);
    bool c_ok(valid_cell(c));
    EXPECT_TRUE(c_ok);
    // if (!c_ok) {
    //   printf("%s:%i cell %lu not ok\n", __FUNCTION__, __LINE__, c.as_id());
    // }
    // Each cell should have the six neighbor faces you'd expect. This is
    // reimplementing the logic in Mesh.h ll. 157-66, and 177-80. If those
    // change, this needs to change.
    auto &faces = m.get_faces(c);
    EXPECT_EQ(faces.size(), 6u);
    // The first face should have the same value as the cell index
    EXPECT_EQ(c.as_id(), faces[0].as_id());
    // The next face should be the XY plane above, it should be nx * ny along
    EXPECT_EQ(c.as_id() + n_yz_faces, faces[1].as_id());
    // Ney, the YZ planes
    // The first YZ face should have the same value as the cell index
    EXPECT_EQ(c.as_id() + n_yz_faces + n_xz_faces, faces[2].as_id());
    // The higher YZ face should just be one further on
    index_t ix, iy, iz;
    std::tie(ix, iy, iz) = linear_to_cartesian(c.as_id(), nx, ny);
    index_t plus_x_cell{cartesian_to_linear(ix + 1, iy, iz, nx, ny)};
    EXPECT_EQ(plus_x_cell, faces[3].as_id());
    // And the XZ faces
    // The first should be the cell id plus the offset
    index_t plus_y_cell{cartesian_to_linear(ix, iy + 1, iz, nx, ny)};
    EXPECT_EQ(plus_y_cell + n_yz_faces, faces[4].as_id());
    // The first should be the cell id plus the offset
    index_t plus_z_cell{cartesian_to_linear(ix, iy, iz + 1, nx, ny)};
    EXPECT_EQ(plus_z_cell + n_yz_faces + n_xz_faces, faces[5].as_id());

    // same thing for vertices
    auto &vs = m.get_vertices(c);
    EXPECT_EQ(vs.size(), 8u);
    // vertices just have the same numbering as the cell that "owns" each
    // vertex. That is, they aren't 'striped' the way faces are.
    // Need to re-cast the cell id into the expanded vertex index set
    index_t vid = cartesian_to_linear(ix, iy, iz, (nx + 1), (ny + 1));
    EXPECT_EQ(vid, vs[0].as_id());
    EXPECT_EQ(vid + vertex_offset_x, vs[1].as_id());
    EXPECT_EQ(vid + vertex_offset_y, vs[2].as_id());
    EXPECT_EQ(vid + vertex_offset_x + vertex_offset_y, vs[3].as_id());
    EXPECT_EQ(vid + vertex_offset_z, vs[4].as_id());
    EXPECT_EQ(vid + vertex_offset_x + vertex_offset_z, vs[5].as_id());
    EXPECT_EQ(vid + vertex_offset_y + vertex_offset_z, vs[6].as_id());
    EXPECT_EQ(vid + vertex_offset_x + vertex_offset_y + vertex_offset_z,
              vs[7].as_id());
  }
  return;
} // TEST(default_mesh_Cartesian_mesh, instantiate_tiny) {

TEST(default_mesh_Cartesian_mesh, instantiate) {
  {
    index_t const nx = 5, ny = 3, nz = 2;
    index_t const n_yz_faces(num_yz_faces(nx, ny, nz));
    index_t const n_xz_faces(num_xz_faces(nx, ny, nz));
    index_t exp_n_cells = nx * ny * nz;
    Cartesian_Mesh m(nx, ny, nz);
    auto &cells = m.cells();
    EXPECT_EQ(cells.size(), exp_n_cells);

    // offsets in face arrays
    index_t const vertex_offset_x = 1;
    index_t const vertex_offset_y = nx + 1;
    index_t const vertex_offset_z = (nx + 1) * (ny + 1);

    auto valid_cell = [&cells](cell_handle_t const &c) {
      return cells.find(c) != cells.end();
    };

    for (index_t i = 0; i < exp_n_cells; ++i) {
      // We expect each cell from [0,n_cells) to be in the mesh and no others
      cell_handle_t c(i);
      bool c_ok(valid_cell(c));
      EXPECT_TRUE(c_ok);
      // if (!c_ok) {
      //   printf("%s:%i cell %lu not ok\n", __FUNCTION__, __LINE__, c.as_id());
      // }
      // Assert that each cell should have the six neighbor faces you'd expect.
      auto &faces = m.get_faces(c);
      EXPECT_EQ(faces.size(), 6u);
      // Faces are numbered as if they belong to a mesh that is one greater
      // than the number of cells in that dimension.
      index_t ix, iy, iz;
      std::tie(ix, iy, iz) = linear_to_cartesian(c.as_id(), nx, ny);
      index_t c0_eff = cartesian_to_linear(ix, iy, iz, (nx + 1), ny);
      bool f0_ok{c0_eff == faces[0].as_id()};
      if (!f0_ok) {
        index_t ix, iy, iz;
        std::tie(ix, iy, iz) = linear_to_cartesian(c.as_id(), nx, ny);
        // printf("%s:%i ix=%lu, iy=%lu, iz=%lu, cell_id = %llu, f0 = %llu\n",
        //        __FUNCTION__, __LINE__, ix, iy, iz, c.as_id(),
        //        faces[0].as_id());
      }
      EXPECT_TRUE(f0_ok);
      // The next face should be the XZ plane above, it should be
      index_t c1_eff = cartesian_to_linear(ix, iy, iz, nx, (ny + 1));
      EXPECT_EQ(c1_eff + n_yz_faces, faces[1].as_id());
      index_t c2_eff = cartesian_to_linear(ix, iy, iz, nx, ny);
      EXPECT_EQ(c2_eff + n_yz_faces + n_xz_faces, faces[2].as_id());
      index_t c3_eff = cartesian_to_linear(ix + 1, iy, iz, (nx + 1), ny);
      EXPECT_EQ(c3_eff, faces[3].as_id());
      index_t c4_eff = cartesian_to_linear(ix, iy + 1, iz, nx, (ny + 1));
      EXPECT_EQ(c4_eff + n_yz_faces, faces[4].as_id());
      index_t c5_eff = cartesian_to_linear(ix, iy, iz + 1, nx, ny);
      EXPECT_EQ(c5_eff + n_yz_faces + n_xz_faces, faces[5].as_id());
      // similar thing for vertices
      auto &vs = m.get_vertices(c);
      EXPECT_EQ(vs.size(), 8u);
      // vertices just have the same numbering as the cell that "owns" each
      // vertex. That is, they aren't 'striped' the way faces are.
      // Need to re-cast the cell id into the expanded vertex index set
      index_t vid = cartesian_to_linear(ix, iy, iz, (nx + 1), (ny + 1));
      EXPECT_EQ(vid, vs[0].as_id());
      EXPECT_EQ(vid + vertex_offset_x, vs[1].as_id());
      EXPECT_EQ(vid + vertex_offset_y, vs[2].as_id());
      EXPECT_EQ(vid + vertex_offset_x + vertex_offset_y, vs[3].as_id());
      EXPECT_EQ(vid + vertex_offset_z, vs[4].as_id());
      EXPECT_EQ(vid + vertex_offset_x + vertex_offset_z, vs[5].as_id());
      EXPECT_EQ(vid + vertex_offset_y + vertex_offset_z, vs[6].as_id());
      EXPECT_EQ(vid + vertex_offset_x + vertex_offset_y + vertex_offset_z,
                vs[7].as_id());
    } // for c in cells
  }   // scope 5x3x2 mesh
  return;
} // TEST(default_mesh_Cartesian_mesh,instantiate){

TEST(default_mesh_Cartesian_mesh, num_faces) {
  index_t nx = 5, ny = 3, nz = 2;
  index_t exp_n_xy_faces = 45;
  index_t n_xy_faces = num_xy_faces(nx, ny, nz);
  EXPECT_EQ(n_xy_faces, exp_n_xy_faces);
  index_t exp_n_yz_faces = 36;
  index_t n_yz_faces = num_yz_faces(nx, ny, nz);
  EXPECT_EQ(n_yz_faces, exp_n_yz_faces);
  index_t exp_n_xz_faces = 40;
  index_t n_xz_faces = num_xz_faces(nx, ny, nz);
  EXPECT_EQ(n_xz_faces, exp_n_xz_faces);
}

TEST(default_mesh_Cartesian_mesh, cartesian_to_linear) {
  /* Case 1: nx = 5, ny = 3, nz = 2 */
  {
    index_t nx = 5, ny = 3;
    index_t ix = 0, iy = 0, iz = 0;
    index_t exp_idx = 0;
    index_t idx = cartesian_to_linear(ix, iy, iz, nx, ny);
    EXPECT_EQ(idx, exp_idx);
    index_t exp_ix = 0, exp_iy = 0, exp_iz = 0;
    std::tie(exp_ix, exp_iy, exp_iz) = linear_to_cartesian(idx, nx, ny);
    EXPECT_EQ(ix, exp_ix);
    EXPECT_EQ(iy, exp_iy);
    EXPECT_EQ(iz, exp_iz);
  } // scope
  {
    index_t nx = 5, ny = 3;
    index_t ix = 1, iy = 0, iz = 0;
    index_t exp_idx = 1;
    index_t idx = cartesian_to_linear(ix, iy, iz, nx, ny);
    EXPECT_EQ(idx, exp_idx);
    index_t exp_ix = 0, exp_iy = 0, exp_iz = 0;
    std::tie(exp_ix, exp_iy, exp_iz) = linear_to_cartesian(idx, nx, ny);
    EXPECT_EQ(ix, exp_ix);
    EXPECT_EQ(iy, exp_iy);
    EXPECT_EQ(iz, exp_iz);
  } // scope
  {
    index_t nx = 5, ny = 3;
    index_t ix = 4, iy = 1, iz = 1;
    index_t exp_idx = 24;
    index_t idx = cartesian_to_linear(ix, iy, iz, nx, ny);
    EXPECT_EQ(idx, exp_idx);
    index_t exp_ix = 0, exp_iy = 0, exp_iz = 0;
    std::tie(exp_ix, exp_iy, exp_iz) = linear_to_cartesian(idx, nx, ny);
    EXPECT_EQ(ix, exp_ix);
    EXPECT_EQ(iy, exp_iy);
    EXPECT_EQ(iz, exp_iz);
  } // scope
  {
    index_t nx = 5, ny = 3;
    index_t ix = 4, iy = 2, iz = 1;
    index_t exp_idx = 29;
    index_t idx = cartesian_to_linear(ix, iy, iz, nx, ny);
    EXPECT_EQ(idx, exp_idx);
    index_t exp_ix = 0, exp_iy = 0, exp_iz = 0;
    std::tie(exp_ix, exp_iy, exp_iz) = linear_to_cartesian(idx, nx, ny);
    EXPECT_EQ(ix, exp_ix);
    EXPECT_EQ(iy, exp_iy);
    EXPECT_EQ(iz, exp_iz);
  } // scope
} // TEST(default_mesh,cartesian_to_linear){

TEST(default_mesh_Cartesian_mesh, get_extents) {
  {
    index_t const nx = 5, ny = 3, nz = 2;
    Cartesian_Mesh m(nx, ny, nz);
    // check x extents
    for (index_t i = 0; i < nx; ++i) {
      index_t idx = cartesian_to_linear(i, 0, 0, nx, ny);
      auto [x_lo, x_hi] = m.get_x_extents(cell_handle_t(idx));
      geom_t x_lo_exp = static_cast<double>(i);
      geom_t x_hi_exp = static_cast<double>(i + 1);
      bool x_lo_ok = soft_equiv_os(x_lo, x_lo_exp, "x_lo");
      bool x_hi_ok = soft_equiv_os(x_hi, x_hi_exp, "x_hi");
      EXPECT_TRUE(x_lo_ok);
      EXPECT_TRUE(x_hi_ok);
    }
    for (index_t i = 0; i < ny; ++i) {
      index_t idx = cartesian_to_linear(0, i, 0, nx, ny);
      auto [y_lo, y_hi] = m.get_y_extents(cell_handle_t(idx));
      geom_t y_lo_exp = static_cast<double>(i);
      geom_t y_hi_exp = static_cast<double>(i + 1);
      bool y_lo_ok = soft_equiv_os(y_lo, y_lo_exp, "y_lo");
      bool y_hi_ok = soft_equiv_os(y_hi, y_hi_exp, "y_hi");
      EXPECT_TRUE(y_lo_ok);
      EXPECT_TRUE(y_hi_ok);
    }
    for (index_t i = 0; i < nz; ++i) {
      index_t idx = cartesian_to_linear(0, 0, i, nx, ny);
      auto [z_lo, z_hi] = m.get_z_extents(cell_handle_t(idx));
      geom_t z_lo_exp = static_cast<double>(i);
      geom_t z_hi_exp = static_cast<double>(i + 1);
      bool z_lo_ok = soft_equiv_os(z_lo, z_lo_exp, "z_lo");
      bool z_hi_ok = soft_equiv_os(z_hi, z_hi_exp, "z_hi");
      EXPECT_TRUE(z_lo_ok);
      EXPECT_TRUE(z_hi_ok);
    }
  }
  // same test, only nor using default mesh sizes
  {
    index_t const nx = 5, ny = 3, nz = 2;
    geom_t dx = 2.0, dy = 3.0, dz = 4.0;
    Cartesian_Mesh m(nx, ny, nz, dx, dy, dz);
    // check x extents
    for (index_t i = 0; i < nx; ++i) {
      index_t idx = cartesian_to_linear(i, 0, 0, nx, ny);
      auto [x_lo, x_hi] = m.get_x_extents(cell_handle_t(idx));
      geom_t x_lo_exp = 2.0 * static_cast<double>(i);
      geom_t x_hi_exp = 2.0 * static_cast<double>(i + 1);
      bool x_lo_ok = soft_equiv_os(x_lo, x_lo_exp, "x_lo");
      bool x_hi_ok = soft_equiv_os(x_hi, x_hi_exp, "x_hi");
      EXPECT_TRUE(x_lo_ok);
      EXPECT_TRUE(x_hi_ok);
    }
    for (index_t i = 0; i < ny; ++i) {
      index_t idx = cartesian_to_linear(0, i, 0, nx, ny);
      auto [y_lo, y_hi] = m.get_y_extents(cell_handle_t(idx));
      geom_t y_lo_exp = 3.0 * static_cast<double>(i);
      geom_t y_hi_exp = 3.0 * static_cast<double>(i + 1);
      bool y_lo_ok = soft_equiv_os(y_lo, y_lo_exp, "y_lo");
      bool y_hi_ok = soft_equiv_os(y_hi, y_hi_exp, "y_hi");
      EXPECT_TRUE(y_lo_ok);
      EXPECT_TRUE(y_hi_ok);
    }
    for (index_t i = 0; i < nz; ++i) {
      index_t idx = cartesian_to_linear(0, 0, i, nx, ny);
      auto [z_lo, z_hi] = m.get_z_extents(cell_handle_t(idx));
      geom_t z_lo_exp = 4.0 * static_cast<double>(i);
      geom_t z_hi_exp = 4.0 * static_cast<double>(i + 1);
      bool z_lo_ok = soft_equiv_os(z_lo, z_lo_exp, "z_lo");
      bool z_hi_ok = soft_equiv_os(z_hi, z_hi_exp, "z_hi");
      EXPECT_TRUE(z_lo_ok);
      EXPECT_TRUE(z_hi_ok);
    }
  }
  return;
} // TEST(default_mesh,get_extents)

TEST(default_mesh_Cartesian_mesh, in_cell) {
  Vector const dir{1.0, 0.0, 0.0};
  cell_handle_t c{1};
  auto ray = [&](Vector const &v) { return geometron{Ray{v, dir}, c}; };
  {
    index_t const nx = 1, ny = 1, nz = 1;
    Cartesian_Mesh const m(nx, ny, nz);
    cell_handle_t c{0};
    Vector const v1{0.5, 0.5, 0.5};
    bool const v1_ok = m.in_cell(ray(v1), c);
    EXPECT_TRUE(v1_ok);

    // in cell if on a lower boundary
    Vector const v2{0.5, 0.0, 0.5};
    bool const v2_ok = m.in_cell(v2, c);
    EXPECT_TRUE(v2_ok);

    // not in cell if on an upper boundary
    Vector const v3{0.5, 1.0, 0.5};
    bool const v3_ok = !m.in_cell(v3, c);
    EXPECT_TRUE(v3_ok);
  }
  // test with non-zero minima
  {
    index_t const nx = 2, ny = 3, nz = 4;
    geom_t const dx = 1.1, dy = 2.2, dz = 3.3;
    geom_t const xmin = 0.2, ymin = 0.0, zmin = -1.0;
    Cartesian_Mesh const m(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
    Vector const v1{0.15, 0.5, 0.5};
    cell_handle_t const c0{0};
    bool const v1_ok{m.in_cell(v1, c0)};
    EXPECT_FALSE(v1_ok);

    Vector const v2{0.25, 0.5, 0.5};
    bool const v2_ok{m.in_cell(v2, c0)};
    EXPECT_TRUE(v2_ok);

    Vector const v3{0.15, 0.5, 2.300000001};
    cell_handle_t const c6{6};
    bool const v3_ok{m.in_cell(v3, c6)};
    EXPECT_FALSE(v3_ok);
  }
} // TEST(default_mesh_Cartesian_mesh,in_cell)

TEST(default_mesh_Cartesian_mesh, intersection) {
  {
    index_t const nx = 1, ny = 1, nz = 1;
    Cartesian_Mesh const m(nx, ny, nz);
    cell_handle_t const c{0};

    // Headed toward HIGH_X
    Ray const r1{{0.5, 0.5, 0.5}, {1.0, 0.0, 0.0}};
    auto const [f1, d1] = m.intersection(r1, c);
    EXPECT_TRUE(f1 == Cartesian_Face(3));
    EXPECT_TRUE(soft_equiv_os(d1, 0.5, "distance to face (1x)"));

    const geom_t s2 = std::sqrt(2.0);
    // If headed directly toward X-Y edge, prefers X
    Ray const r2{{0.5, 0.5, 0.5}, {1.0 / s2, 1.0 / s2, 0.0}};
    auto const [f2, d2] = m.intersection(r2, c);
    EXPECT_TRUE(f2 == Cartesian_Face(3));
    EXPECT_TRUE(soft_equiv_os(d2, s2 / 2.0, "distance to face (2)"));

    const geom_t s3 = std::sqrt(3.0);
    // If headed directly toward X-Y-Z corner, prefers X
    Ray const r3{{0.5, 0.5, 0.5}, {1.0 / s3, 1.0 / s3, 1.0 / s3}};
    auto const [f3, d3] = m.intersection(r3, c);
    EXPECT_TRUE(f3 == Cartesian_Face(3));
    EXPECT_TRUE(soft_equiv_os(d3, s3 / 2.0, "distance to face (3)"));

    // If headed directly toward Y-Z edge, prefers Y
    Ray const r4{{0.5, 0.5, 0.5}, {0.0, 1.0 / s2, 1.0 / s2}};
    auto const [f4, d4] = m.intersection(r4, c);
    EXPECT_TRUE(f4 == Cartesian_Face(4));
    EXPECT_TRUE(soft_equiv_os(d4, s2 / 2.0, "distance to face (4)"));

    // Headed toward LOW_Z
    Ray const r5{{0.5, 0.5, 0.5}, {0.0, 0.0, -1.0}};
    auto const [f5, d5] = m.intersection(r5, c);
    EXPECT_TRUE(f5 == Cartesian_Face(2));
    EXPECT_TRUE(soft_equiv_os(d5, 0.5, "distance to face (5)"));

    // Headed toward LOW_Y
    Ray const r6{{0.5, 0.5, 0.5}, {0.0, -1.0, 0.0}};
    auto const [f6, d6] = m.intersection(r6, c);
    EXPECT_TRUE(f6 == Cartesian_Face(1));
    EXPECT_TRUE(soft_equiv_os(d6, 0.5, "distance to face (6)"));

    // Headed toward LOW_X
    Ray const r7{{0.5, 0.5, 0.5}, {-1.0, 0.0, 0.0}};
    auto const [f7, d7] = m.intersection(r7, c);
    EXPECT_TRUE(f7 == Cartesian_Face(0));
    EXPECT_TRUE(soft_equiv_os(d7, 0.5, "distance to face (7)"));

    // Headed toward HIGH_Z
    Ray const r8{{0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}};
    auto const [f8, d8] = m.intersection(r8, c);
    EXPECT_TRUE(f8 == Cartesian_Face(5));
    EXPECT_TRUE(soft_equiv_os(d8, 0.5, "distance to face (8)"));

    // Headed toward HIGH_Y
    Ray const r9{{0.5, 0.5, 0.5}, {0.0, 1.0, 0.0}};
    auto const [f9, d9] = m.intersection(r9, c);
    EXPECT_TRUE(f9 == Cartesian_Face(4));
    EXPECT_TRUE(soft_equiv_os(d9, 0.5, "distance to face (9)"));

    // Headed toward LOW_X, LOW_Y corner
    Ray const r10{{0.5, 0.5, 0.5}, {-1.0 / s2, -1.0 / s2, 0.0}};
    auto const [f10, d10] = m.intersection(r10, c);
    EXPECT_TRUE(f10 == Cartesian_Face(0));
    EXPECT_TRUE(soft_equiv_os(d10, s2 / 2.0, "distance to face (10)"));

    // Headed toward LOW_X, LOW_Y corner
    Ray const r11{{0.5, 0.5, 0.5}, {0.0, -1.0 / s2, -1.0 / s2}};
    auto const [f11, d11] = m.intersection(r11, c);
    EXPECT_TRUE(f11 == Cartesian_Face(1));
    EXPECT_TRUE(soft_equiv_os(d11, s2 / 2.0, "distance to face (11)"));
  }
} // TEST(default_mesh_Cartesian_mesh,intersection)

// Here is a list of randomly generated tests for example 2 in the next test
size_t constexpr n_example2_points{20};
extern const std::tuple<index_t, index_t, index_t, Cartesian_Mesh::Point>
    example2_points[n_example2_points];

size_t constexpr n_example3_points_a{300};
extern const std::tuple<index_t, index_t, index_t, Cartesian_Mesh::Point>
    example3_points_a[n_example3_points_a];

TEST(default_mesh_Cartesian_mesh, find_cell) {
  // example one: one cell mesh
  {
    index_t const nx = 1, ny = 1, nz = 1;
    Cartesian_Mesh const m(nx, ny, nz);
    cell_handle_t c0{0};
    Point p1{0.25, 0.25, 0.25};
    auto c1{m.find_cell(p1)};
    EXPECT_TRUE(m.in_cell(p1, c1));
    EXPECT_EQ(c1, cell_handle_t(0));

    Point p2{0.0, 0.0, 0.0};
    auto c2{m.find_cell(p2)};
    EXPECT_TRUE(m.in_cell(p2, c2));
    EXPECT_EQ(c2, cell_handle_t(0));

  } // scope (one-cell mesh)
  // example two: two cell mesh
  {
    index_t const nx = 1, ny = 1, nz = 2;
    geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
    geom_t const xmin = 0.0, ymin = 0.0, zmin = -1.0;
    Cartesian_Mesh const m(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
    cell_handle_t const c0{0}, c1{1};
    Point const p1{0.25, 0.25, 0.25};
    auto const result1{m.find_cell(p1)};
    EXPECT_TRUE(m.in_cell(p1, result1));
    EXPECT_TRUE(m.in_cell(p1, c0));
    EXPECT_EQ(result1, c0);

    Point const p2{0.01, 0.01, 3.9};
    auto const result2{m.find_cell(p2)};
    EXPECT_TRUE(m.in_cell(p2, result2));
    EXPECT_TRUE(m.in_cell(p2, c1));
    EXPECT_EQ(result2, c1);

    for (size_t i = 0; i < n_example2_points; ++i) {
      auto [ix, iy, iz, p] = example2_points[i];
      cell_handle_t expected_cell{m.make_cell(ix, iy, iz)};
      auto const result{m.find_cell(p)};
      EXPECT_TRUE(m.in_cell(p, result));
      EXPECT_TRUE(m.in_cell(p, expected_cell));
      EXPECT_EQ(result, expected_cell);
    }
  } // scope (two-cell mesh)
  // example 3: many cell mesh
  {
    index_t const nx = 2, ny = 3, nz = 5;
    geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
    geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
    Cartesian_Mesh const m(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
    for (size_t i = 0; i < n_example3_points_a; ++i) {
      auto [ix, iy, iz, p] = example3_points_a[i];
      cell_handle_t expected_cell{m.make_cell(ix, iy, iz)};
      auto const result{m.find_cell(p)};
      EXPECT_TRUE(m.in_cell(p, result));
      EXPECT_TRUE(m.in_cell(p, expected_cell));
      EXPECT_EQ(result, expected_cell);
    }
  }
} // TEST(default_mesh_Cartesian_mesh,find_cell)

TEST(default_mesh_Cartesian_mesh, cell_across) {
  {
    index_t const nx = 3, ny = 5, nz = 7;
    geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
    geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
    Cartesian_Mesh const m(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
    Cartesian_Mesh::Point punused{0.0, 0.0, 0.0};
    cell_handle_t const &null_cell = Cartesian_Mesh::null_cell();
    auto run_test = [&](cell_handle_t const &cin, Cartesian_Face const &f,
                        cell_handle_t const &cexp, const char *name) {
      SCOPED_TRACE(name);
      cell_handle_t ca = m.cell_across(cin, f, punused);
      EXPECT_EQ(ca, cexp);
    };
    {
      cell_handle_t c0{0};
      run_test(c0, f_l_x, null_cell, "c0 low x");
      run_test(c0, f_l_y, null_cell, "c0 low y");
      run_test(c0, f_l_z, null_cell, "c0 low z");
      run_test(c0, f_h_x, cell_handle_t{1}, "c0 high x");
      run_test(c0, f_h_y, cell_handle_t{3}, "c0 high y");
      run_test(c0, f_h_z, cell_handle_t{15}, "c0 high z");
    }
    {
      cell_handle_t c4{4};
      run_test(c4, f_l_x, cell_handle_t{3}, "c4 low x");
      run_test(c4, f_l_y, cell_handle_t{1}, "c4 low y");
      run_test(c4, f_l_z, null_cell, "c4 low z");
      run_test(c4, f_h_x, cell_handle_t{5}, "c4 high x");
      run_test(c4, f_h_y, cell_handle_t{7}, "c4 high y");
      run_test(c4, f_h_z, cell_handle_t{19}, "c4 high z");
    }
    {
      cell_handle_t c19{19};
      run_test(c19, f_l_x, cell_handle_t{18}, "c19 low x");
      run_test(c19, f_l_y, cell_handle_t{16}, "c19 low y");
      run_test(c19, f_l_z, cell_handle_t{4}, "c19 low z");
      run_test(c19, f_h_x, cell_handle_t{20}, "c19 high x");
      run_test(c19, f_h_y, cell_handle_t{22}, "c19 high y");
      run_test(c19, f_h_z, cell_handle_t{34}, "c19 high z");
    }
    {
      cell_handle_t c104{104};
      run_test(c104, f_l_x, cell_handle_t{103}, "c104 low x");
      run_test(c104, f_l_y, cell_handle_t{101}, "c104 low y");
      run_test(c104, f_l_z, cell_handle_t{89}, "c104 low z");
      run_test(c104, f_h_x, null_cell, "c104 high x");
      run_test(c104, f_h_y, null_cell, "c104 high y");
      run_test(c104, f_h_z, null_cell, "c104 high z");
    }
  }
  return;
} // TEST(default_mesh_Cartesian_mesh,cell_across)

TEST(default_mesh_Cartesian_mesh, is_boundary_cell_face_3_cell_line) {
  index_t const nx = 1, ny = 3, nz = 1;
  // geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
  // geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
  Cartesian_Mesh const m(nx, ny, nz);

  {
    cell_handle_t c0{0};
    auto const &faces{m.get_faces(c0)};
    // per face, only the high-y is not a boundary face in this cell
    EXPECT_TRUE(m.is_boundary(faces[0]));
    EXPECT_TRUE(m.is_boundary(faces[1]));
    EXPECT_TRUE(m.is_boundary(faces[2]));
    EXPECT_TRUE(m.is_boundary(faces[3]));
    EXPECT_FALSE(m.is_boundary(faces[4]));
    EXPECT_TRUE(m.is_boundary(faces[5]));
  }
  {
    cell_handle_t c1{1};
    auto const &faces{m.get_faces(c1)};
    // per face, only the low-y & high-y are not a boundary face in this cell
    EXPECT_TRUE(m.is_boundary(faces[0]));
    EXPECT_FALSE(m.is_boundary(faces[1]));
    EXPECT_TRUE(m.is_boundary(faces[2]));
    EXPECT_TRUE(m.is_boundary(faces[3]));
    EXPECT_FALSE(m.is_boundary(faces[4]));
    EXPECT_TRUE(m.is_boundary(faces[5]));
  }
  {
    cell_handle_t c2{2};
    auto const &faces{m.get_faces(c2)};
    // per face, only the low-y is not a boundary face in this cell
    EXPECT_TRUE(m.is_boundary(faces[0]));
    EXPECT_FALSE(m.is_boundary(faces[1]));
    EXPECT_TRUE(m.is_boundary(faces[2]));
    EXPECT_TRUE(m.is_boundary(faces[3]));
    EXPECT_TRUE(m.is_boundary(faces[4]));
    EXPECT_TRUE(m.is_boundary(faces[5]));
  }
  return;
}

TEST(default_mesh_Cartesian_mesh, is_boundary_cell_large_mesh) {

  index_t const nx = 3, ny = 5, nz = 7;
  // geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
  // geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
  Cartesian_Mesh const m(nx, ny, nz /*, dx, dy, dz, xmin, ymin, zmin*/);

  // interior cells: is_boundary false
  for (index_t iz = 1; iz < (nz - 1); ++iz) {
    for (index_t iy = 1; iy < (ny - 1); ++iy) {
      for (index_t ix = 1; ix < (nx - 1); ++ix) {
        std::stringstream serr;
        serr << "ix = " << ix << ", iy = " << iy << ", iz = " << iz;
        SCOPED_TRACE(serr.str());
        // index_t cidx = cartesian_to_linear(ix, iy, iz, nx, ny);
      }
    }
  }

  // exterior cells: is_boundary is true
  for (index_t iz = 1; iz < (nz - 1); ++iz) {
    for (index_t iy = 1; iy < (ny - 1); ++iy) {
      {
        std::stringstream serr;
        serr << "ix = " << 0 << ", iy = " << iy << ", iz = " << iz;
        SCOPED_TRACE(serr.str());
      }
      {
        std::stringstream serr;
        serr << "ix = " << (nx - 1) << ", iy = " << iy << ", iz = " << iz;
        SCOPED_TRACE(serr.str());
      }
    }
  }
  for (index_t ix = 1; ix < (nx - 1); ++ix) {
    for (index_t iy = 1; iy < (ny - 1); ++iy) {
      {
        std::stringstream serr;
        serr << "ix = " << ix << ", iy = " << iy << ", iz = " << 0;
        SCOPED_TRACE(serr.str());
      }
      {
        std::stringstream serr;
        serr << "ix = " << ix << ", iy = " << iy << ", iz = " << (nz - 1);
        SCOPED_TRACE(serr.str());
      }
    }
  }
  for (index_t iz = 1; iz < (nz - 1); ++iz) {
    for (index_t ix = 1; ix < (nx - 1); ++ix) {
      {
        std::stringstream serr;
        serr << "ix = " << ix << ", iy = " << 0 << ", iz = " << iz;
        SCOPED_TRACE(serr.str());
      }
      {
        std::stringstream serr;
        serr << "ix = " << ix << ", iy = " << (ny - 1) << ", iz = " << iz;
        SCOPED_TRACE(serr.str());
      }
    }
  }
  return;
} // TEST(default_mesh_Cartesian_mesh,is_boundary)

TEST(default_mesh_Cartesian_mesh, is_boundary_cell_face_large_mesh) {
  index_t const nx = 3, ny = 5, nz = 7;
  geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
  geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
  Cartesian_Mesh const m(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);

  // interior cells: is_boundary false
  for (index_t iz = 1; iz < (nz - 1); ++iz) {
    for (index_t iy = 1; iy < (ny - 1); ++iy) {
      for (index_t ix = 1; ix < (nx - 1); ++ix) {
        std::stringstream serr;
        serr << "ix = " << ix << ", iy = " << iy << ", iz = " << iz;
        SCOPED_TRACE(serr.str());
        index_t cidx = cartesian_to_linear(ix, iy, iz, nx, ny);
        cell_handle_t c{cidx};
        auto const &faces{m.get_faces(c)};
        for (auto &f : faces) {
          EXPECT_FALSE(m.is_boundary(f));
        }
      } // for ix
    }   // for iy
  }     // for iz

  // boundary cells: is_boundary is true for some, false for others
  // low-x & high-x faces
  for (index_t iz = 1; iz < (nz - 1); ++iz) {
    for (index_t iy = 1; iy < (ny - 1); ++iy) {
      {
        std::stringstream serr;
        serr << "ix = " << 0 << ", iy = " << iy << ", iz = " << iz;
        SCOPED_TRACE(serr.str());
        index_t cidx = cartesian_to_linear(0, iy, iz, nx, ny);
        cell_handle_t c{cidx};
        auto const &faces{m.get_faces(c)};
        EXPECT_TRUE(m.is_boundary(faces[0]));
        EXPECT_FALSE(m.is_boundary(faces[1]));
        EXPECT_FALSE(m.is_boundary(faces[2]));
        EXPECT_FALSE(m.is_boundary(faces[3]));
        EXPECT_FALSE(m.is_boundary(faces[4]));
        EXPECT_FALSE(m.is_boundary(faces[5]));
      }
      {
        std::stringstream serr;
        serr << "ix = " << (nx - 1) << ", iy = " << iy << ", iz = " << iz;
        SCOPED_TRACE(serr.str());
        index_t cidx = cartesian_to_linear((nx - 1), iy, iz, nx, ny);
        cell_handle_t c{cidx};
        auto const &faces{m.get_faces(c)};
        EXPECT_FALSE(m.is_boundary(faces[0]));
        EXPECT_FALSE(m.is_boundary(faces[1]));
        EXPECT_FALSE(m.is_boundary(faces[2]));
        EXPECT_TRUE(m.is_boundary(faces[3]));
        EXPECT_FALSE(m.is_boundary(faces[4]));
        EXPECT_FALSE(m.is_boundary(faces[5]));
      }
    } // for iy
  }   // for iz
  // low-z and high-z faces
  for (index_t ix = 1; ix < (nx - 1); ++ix) {
    for (index_t iy = 1; iy < (ny - 1); ++iy) {
      {
        std::stringstream serr;
        serr << "ix = " << ix << ", iy = " << iy << ", iz = " << 0;
        SCOPED_TRACE(serr.str());
        index_t cidx = cartesian_to_linear(ix, iy, 0, nx, ny);
        cell_handle_t c{cidx};
        auto const &faces{m.get_faces(c)};
        EXPECT_FALSE(m.is_boundary(faces[0]));
        EXPECT_FALSE(m.is_boundary(faces[1]));
        EXPECT_TRUE(m.is_boundary(faces[2]));
        EXPECT_FALSE(m.is_boundary(faces[3]));
        EXPECT_FALSE(m.is_boundary(faces[4]));
        EXPECT_FALSE(m.is_boundary(faces[5]));
      }
      {
        std::stringstream serr;
        serr << "ix = " << ix << ", iy = " << iy << ", iz = " << (nz - 1);
        SCOPED_TRACE(serr.str());
        index_t cidx = cartesian_to_linear(ix, iy, (nz - 1), nx, ny);
        cell_handle_t c{cidx};
        auto const &faces{m.get_faces(c)};
        EXPECT_FALSE(m.is_boundary(faces[0]));
        EXPECT_FALSE(m.is_boundary(faces[1]));
        EXPECT_FALSE(m.is_boundary(faces[2]));
        EXPECT_FALSE(m.is_boundary(faces[3]));
        EXPECT_FALSE(m.is_boundary(faces[4]));
        EXPECT_TRUE(m.is_boundary(faces[5]));
      }
    } // for iy
  }   // for ix
  // low-y and high-y faces
  for (index_t iz = 1; iz < (nz - 1); ++iz) {
    for (index_t ix = 1; ix < (nx - 1); ++ix) {
      {
        std::stringstream serr;
        serr << "ix = " << ix << ", iy = " << 0 << ", iz = " << iz;
        SCOPED_TRACE(serr.str());
        index_t cidx = cartesian_to_linear(ix, 0, iz, nx, ny);
        cell_handle_t c{cidx};
        auto const &faces{m.get_faces(c)};
        EXPECT_FALSE(m.is_boundary(faces[0]));
        EXPECT_TRUE(m.is_boundary(faces[1]));
        EXPECT_FALSE(m.is_boundary(faces[2]));
        EXPECT_FALSE(m.is_boundary(faces[3]));
        EXPECT_FALSE(m.is_boundary(faces[4]));
        EXPECT_FALSE(m.is_boundary(faces[5]));
      }
      {
        std::stringstream serr;
        serr << "ix = " << ix << ", iy = " << (ny - 1) << ", iz = " << iz;
        SCOPED_TRACE(serr.str());
        index_t cidx = cartesian_to_linear(ix, (ny - 1), iz, nx, ny);
        cell_handle_t c{cidx};
        auto const &faces{m.get_faces(c)};
        EXPECT_FALSE(m.is_boundary(faces[0]));
        EXPECT_FALSE(m.is_boundary(faces[1]));
        EXPECT_FALSE(m.is_boundary(faces[2]));
        EXPECT_FALSE(m.is_boundary(faces[3]));
        EXPECT_TRUE(m.is_boundary(faces[4]));
        EXPECT_FALSE(m.is_boundary(faces[5]));
      }
    } // for ix
  }   // for iz
  return;
} // TEST(default_mesh_Cartesian_mesh,is_boundary_cell_face_large_mesh)

TEST(default_mesh_Cartesian_mesh, sample_position) {
  geom_t urds[]{0.1, 0.2, 0.3};

  index_t const nx = 3, ny = 5, nz = 7;
  geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
  geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
  Cartesian_Mesh const m(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
  /* Run sample_position in a number of cells. Since the buffer RNG is simple,
   * we can easily calculate what we should get. Notionally, could reuse the
   * same RNG, since it should just roll over every three calls. But, this
   * is simpler. */
  {
    Buffer_RNG<3> rng(&urds[0], &urds[3]);
    cell_handle_t c0{0};
    auto p = m.sample_position(rng, c0);
    bool p0_ok = soft_equiv_os(p[0], 10.25, "x-component");
    bool p1_ok = soft_equiv_os(p[1], 5.72, "y-component");
    bool p2_ok = soft_equiv_os(p[2], 0.41, "z component");
    EXPECT_TRUE(p0_ok);
    EXPECT_TRUE(p1_ok);
    EXPECT_TRUE(p2_ok);
  }
  {
    Buffer_RNG<3> rng(&urds[0], &urds[3]);
    cell_handle_t c1{1};
    auto p = m.sample_position(rng, c1);
    bool p0_ok = soft_equiv_os(p[0], 12.75, "x-component");
    bool p1_ok = soft_equiv_os(p[1], 5.72, "y-component");
    bool p2_ok = soft_equiv_os(p[2], 0.41, "z component");
    EXPECT_TRUE(p0_ok);
    EXPECT_TRUE(p1_ok);
    EXPECT_TRUE(p2_ok);
  }
  {
    Buffer_RNG<3> rng(&urds[0], &urds[3]);
    cell_handle_t c3{3};
    auto p = m.sample_position(rng, c3);
    bool p0_ok = soft_equiv_os(p[0], 10.25, "x-component");
    bool p1_ok = soft_equiv_os(p[1], 9.32, "y-component");
    bool p2_ok = soft_equiv_os(p[2], 0.41, "z component");
    EXPECT_TRUE(p0_ok);
    EXPECT_TRUE(p1_ok);
    EXPECT_TRUE(p2_ok);
  }
  {
    Buffer_RNG<3> rng(&urds[0], &urds[3]);
    cell_handle_t c18{18};
    auto p = m.sample_position(rng, c18);
    bool p0_ok = soft_equiv_os(p[0], 10.25, "x-component");
    bool p1_ok = soft_equiv_os(p[1], 9.32, "y-component");
    bool p2_ok = soft_equiv_os(p[2], 5.11, "z component");
    EXPECT_TRUE(p0_ok);
    EXPECT_TRUE(p1_ok);
    EXPECT_TRUE(p2_ok);
  }
  {
    Buffer_RNG<3> rng(&urds[0], &urds[3]);
    cell_handle_t c104{104};
    auto p = m.sample_position(rng, c104);
    bool p0_ok = soft_equiv_os(p[0], 15.25, "x-component");
    bool p1_ok = soft_equiv_os(p[1], 20.119999999999997, "y-component");
    bool p2_ok = soft_equiv_os(p[2], 28.610000000000003, "z component");
    EXPECT_TRUE(p0_ok);
    EXPECT_TRUE(p1_ok);
    EXPECT_TRUE(p2_ok);
  }
  return;
}

TEST(default_mesh_Cartesian_mesh, sample_direction) {
  geom_t urds[]{0.5, 0.5};

  index_t const nx = 3, ny = 5, nz = 7;
  geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
  geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
  Cartesian_Mesh const m(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
  {
    Buffer_RNG<3> rng(&urds[0], &urds[2]);
    auto omega = m.sample_direction_isotropic(rng);
    /* ctheta = 0
       stheta = 1
       phi = pi, cphi = -1, sphi = 0
     */
    bool omega_x_ok = soft_equiv_os(omega[0], -1.0, "omega_x");
    bool omega_y_ok = soft_equiv_os(omega[1], 0.0, "omega_y");
    bool omega_z_ok = soft_equiv_os(omega[2], 0.0, "omega_z");
    EXPECT_TRUE(omega_x_ok);
    EXPECT_TRUE(omega_y_ok);
    EXPECT_TRUE(omega_z_ok);
  }
  return;
}

TEST(default_mesh_Cartesian_mesh, face_to_face_name) {
  Cartesian_Mesh m{one_cell_cmi()};
  cell_handle_t c0{0u};
  auto faces{m.get_faces(c0)};
  Face_Name fname{m.face_to_face_name(c0, faces[0])};
  EXPECT_EQ(fname, Cartesian_Mesh::LOW_X);
  return;
}

// Randomly generated {face_handle_t,Vector} pairs for test of reflect.
constexpr size_t n_reflection_cases{600};
extern const std::tuple<nut_mesh::Cartesian_Face, Cartesian_Mesh::Vector,
                        Cartesian_Mesh::Vector>
    test_reflection_cases[n_reflection_cases];

TEST(default_mesh_Cartesian_mesh, reflect) {
  index_t const nx = 1, ny = 1, nz = 1;
  geom_t const dx = 2000.0, dy = 2000.0, dz = 2000.0;
  geom_t const xmin = -1000.0, ymin = -1000.0, zmin = -1000.0;
  Cartesian_Mesh m{nx, ny, nz, dx, dy, dz, xmin, ymin, zmin};
  cell_handle_t c0{0u};
  uint32_t n_fails(0);
  for (size_t i = 0; i < n_reflection_cases; ++i) {
    // for(size_t i = n_reflection_cases - 4; i < n_reflection_cases; ++i) {
    // for(size_t i = 0; i < 3; ++i) {
    std::string const label(std::to_string(i) + " ");
    auto [f, v, r_exp] = test_reflection_cases[i];
    Vector r{m.reflect(f, v)};
    bool refl_ok{
        vec_soft_equiv(r, r_exp, "Cartesian_Mesh::reflect case " + label)};
    if (!refl_ok) {
      n_fails++;
      printf("%s:%i That was case %lu, n_fails = %u\n", __FUNCTION__, __LINE__,
             i, n_fails);
    }
    EXPECT_TRUE(refl_ok);
    // The parallel component of the reflected vector should be equal and
    // opposite to the original.
    Vector normal = m.get_normal(c0, f);
    bool comps_ok =
        soft_equiv_os(v.dot(normal), -r.dot(normal), "normals", 1.0e-13);
    if (!comps_ok) {
      n_fails++;
      printf("%s:%i That was case %lu, n_fails = %u\n", __FUNCTION__, __LINE__,
             i, n_fails);
    }
    if (n_fails > 2) {
      break;
    }
    EXPECT_TRUE(comps_ok);
  } // for case in cases
  return;
} // TEST(default_mesh_Cartesian_mesh, reflect)

//          --------- test problem data ---------

nut_mesh::Cartesian_Face face_lx{0};
nut_mesh::Cartesian_Face face_hx{1};
nut_mesh::Cartesian_Face face_ly{2};
nut_mesh::Cartesian_Face face_hy{3};
nut_mesh::Cartesian_Face face_lz{4};
nut_mesh::Cartesian_Face face_hz{5};

const std::tuple<nut_mesh::Cartesian_Face, Cartesian_Mesh::Vector,
                 Cartesian_Mesh::Vector>
    test_reflection_cases[n_reflection_cases] = {
        {face_lx,
         {-209.751039716082, 225.626549503419, -108.744639537191},
         {209.751039716082, 225.626549503419, -108.744639537191}},
        {face_lx,
         {-369.792022654238, 321.068562195155, 947.913905568476},
         {369.792022654238, 321.068562195155, 947.913905568476}},
        {face_lx,
         {-433.411678869332, 11.7866053203579, 982.002116902816},
         {433.411678869332, 11.7866053203579, 982.002116902816}},
        {face_lx,
         {-294.842339136222, -88.4856187499113, -807.476822674811},
         {294.842339136222, -88.4856187499113, -807.476822674811}},
        {face_lx,
         {-266.370799722577, -412.469114488175, -970.264615794086},
         {266.370799722577, -412.469114488175, -970.264615794086}},
        {face_lx,
         {-839.48856756651, -304.201707881154, -587.251086338735},
         {839.48856756651, -304.201707881154, -587.251086338735}},
        {face_lx,
         {-339.924319552552, -924.425464183352, 809.512552296111},
         {339.924319552552, -924.425464183352, 809.512552296111}},
        {face_lx,
         {-413.58120010526, 291.091841128781, 430.65274578935},
         {413.58120010526, 291.091841128781, 430.65274578935}},
        {face_lx,
         {-749.830821928774, -284.372764811067, 347.535731523861},
         {749.830821928774, -284.372764811067, 347.535731523861}},
        {face_lx,
         {-615.686153775672, -447.632280240663, 628.452352463717},
         {615.686153775672, -447.632280240663, 628.452352463717}},
        {face_lx,
         {-916.050929190633, -883.116983711704, 330.279454141775},
         {916.050929190633, -883.116983711704, 330.279454141775}},
        {face_lx,
         {-868.336176628663, -772.57702816734, 541.278838853745},
         {868.336176628663, -772.57702816734, 541.278838853745}},
        {face_lx,
         {-545.008606777656, -135.954502759918, 168.513654103909},
         {545.008606777656, -135.954502759918, 168.513654103909}},
        {face_lx,
         {-365.639403596436, 526.607768493771, 991.320917006623},
         {365.639403596436, 526.607768493771, 991.320917006623}},
        {face_lx,
         {-82.0053905971536, -633.673872303761, 870.906934644928},
         {82.0053905971536, -633.673872303761, 870.906934644928}},
        {face_lx,
         {-755.619684679741, -978.299381762029, -785.362911980641},
         {755.619684679741, -978.299381762029, -785.362911980641}},
        {face_lx,
         {-740.088299062813, -860.264585144716, -813.773312561387},
         {740.088299062813, -860.264585144716, -813.773312561387}},
        {face_lx,
         {-982.826911109105, 173.444473301614, 497.792582746816},
         {982.826911109105, 173.444473301614, 497.792582746816}},
        {face_lx,
         {-616.458675631643, 571.925084418831, -991.744895117235},
         {616.458675631643, 571.925084418831, -991.744895117235}},
        {face_lx,
         {-302.858815001772, -387.750525986516, -89.637128988395},
         {302.858815001772, -387.750525986516, -89.637128988395}},
        {face_lx,
         {-991.198857089184, 482.195959508328, 728.811640017856},
         {991.198857089184, 482.195959508328, 728.811640017856}},
        {face_lx,
         {-286.427859201278, -124.574223516706, 882.295476507225},
         {286.427859201278, -124.574223516706, 882.295476507225}},
        {face_lx,
         {-532.294391053826, -893.897363236512, -858.710989147761},
         {532.294391053826, -893.897363236512, -858.710989147761}},
        {face_lx,
         {-145.063984665691, -330.638229604477, 936.863035697445},
         {145.063984665691, -330.638229604477, 936.863035697445}},
        {face_lx,
         {-723.701451467513, 708.569412169698, -632.774603957642},
         {723.701451467513, 708.569412169698, -632.774603957642}},
        {face_lx,
         {-544.084979773453, 215.243961155508, 81.1910345273127},
         {544.084979773453, 215.243961155508, 81.1910345273127}},
        {face_lx,
         {-733.664164088126, 733.157757777705, -838.424754872637},
         {733.664164088126, 733.157757777705, -838.424754872637}},
        {face_lx,
         {-981.271128470444, 72.5297487805819, 296.997324885753},
         {981.271128470444, 72.5297487805819, 296.997324885753}},
        {face_lx,
         {-117.991125248945, -198.867029243948, 79.7413821786786},
         {117.991125248945, -198.867029243948, 79.7413821786786}},
        {face_lx,
         {-271.824310624947, -711.692829255442, -859.192333886625},
         {271.824310624947, -711.692829255442, -859.192333886625}},
        {face_lx,
         {-849.90897605065, -718.904549530926, -996.014820515414},
         {849.90897605065, -718.904549530926, -996.014820515414}},
        {face_lx,
         {-275.311082665276, -730.043103111127, -294.823440710374},
         {275.311082665276, -730.043103111127, -294.823440710374}},
        {face_lx,
         {-29.2493225471967, -587.019313849478, -102.676825991553},
         {29.2493225471967, -587.019313849478, -102.676825991553}},
        {face_lx,
         {-821.670251016918, -818.04419301563, 16.7387292922381},
         {821.670251016918, -818.04419301563, 16.7387292922381}},
        {face_lx,
         {-341.054989312576, -723.793807824651, 620.248498585467},
         {341.054989312576, -723.793807824651, 620.248498585467}},
        {face_lx,
         {-990.042259780616, 289.134096996485, 282.829817405029},
         {990.042259780616, 289.134096996485, 282.829817405029}},
        {face_lx,
         {-479.83021505272, -283.036378947297, -692.057663400875},
         {479.83021505272, -283.036378947297, -692.057663400875}},
        {face_lx,
         {-965.232751145869, 925.300474254114, -550.983340185287},
         {965.232751145869, 925.300474254114, -550.983340185287}},
        {face_lx,
         {-210.30711850938, -59.94532483079, 623.268996323146},
         {210.30711850938, -59.94532483079, 623.268996323146}},
        {face_lx,
         {-445.08240367813, -967.039335344599, 140.752560139056},
         {445.08240367813, -967.039335344599, 140.752560139056}},
        {face_lx,
         {-344.620238251995, -559.574333846159, -434.175322299135},
         {344.620238251995, -559.574333846159, -434.175322299135}},
        {face_lx,
         {-161.544979080601, -655.886281384313, 918.924358072439},
         {161.544979080601, -655.886281384313, 918.924358072439}},
        {face_lx,
         {-659.09637501303, -550.558077398378, 555.772855078664},
         {659.09637501303, -550.558077398378, 555.772855078664}},
        {face_lx,
         {-584.080088827097, -817.224875876825, 991.408421705497},
         {584.080088827097, -817.224875876825, 991.408421705497}},
        {face_lx,
         {-394.527492543364, 257.798164804549, 753.521867190832},
         {394.527492543364, 257.798164804549, 753.521867190832}},
        {face_lx,
         {-705.121040291238, 288.494782455992, 626.488820540528},
         {705.121040291238, 288.494782455992, 626.488820540528}},
        {face_lx,
         {-496.011021573509, -665.558117616358, 389.258072372487},
         {496.011021573509, -665.558117616358, 389.258072372487}},
        {face_lx,
         {-865.030035205379, -862.665994388272, -526.769196632989},
         {865.030035205379, -862.665994388272, -526.769196632989}},
        {face_lx,
         {-672.794591092986, -898.474931581002, 422.063530880747},
         {672.794591092986, -898.474931581002, 422.063530880747}},
        {face_lx,
         {-85.9692451650353, -860.802851506662, -541.296145631319},
         {85.9692451650353, -860.802851506662, -541.296145631319}},
        {face_lx,
         {-934.880384228854, 127.700614529907, 275.537418760202},
         {934.880384228854, 127.700614529907, 275.537418760202}},
        {face_lx,
         {-221.718560079821, -948.643544985977, 286.367427944383},
         {221.718560079821, -948.643544985977, 286.367427944383}},
        {face_lx,
         {-820.675668871493, -933.321604920981, 630.235139179669},
         {820.675668871493, -933.321604920981, 630.235139179669}},
        {face_lx,
         {-871.051627421902, 25.7178753267217, -33.4993945618189},
         {871.051627421902, 25.7178753267217, -33.4993945618189}},
        {face_lx,
         {-801.093140178381, 741.885323335845, -796.085377563094},
         {801.093140178381, 741.885323335845, -796.085377563094}},
        {face_lx,
         {-580.65281737756, -338.723880264011, -998.677670777857},
         {580.65281737756, -338.723880264011, -998.677670777857}},
        {face_lx,
         {-935.325175382392, -143.799986074188, -485.634727008253},
         {935.325175382392, -143.799986074188, -485.634727008253}},
        {face_lx,
         {-370.771653285604, -957.962868060698, 171.64821295505},
         {370.771653285604, -957.962868060698, 171.64821295505}},
        {face_lx,
         {-467.957990275341, 56.7960019230632, 967.568104241603},
         {467.957990275341, 56.7960019230632, 967.568104241603}},
        {face_lx,
         {-534.063197563785, 937.435043747504, 190.873033218684},
         {534.063197563785, 937.435043747504, 190.873033218684}},
        {face_lx,
         {-595.973392681638, 997.448089704442, -546.902620241847},
         {595.973392681638, 997.448089704442, -546.902620241847}},
        {face_lx,
         {-117.74012416235, -950.335752910494, 977.936315280482},
         {117.74012416235, -950.335752910494, 977.936315280482}},
        {face_lx,
         {-657.360612523305, -609.334953589731, 662.58710020125},
         {657.360612523305, -609.334953589731, 662.58710020125}},
        {face_lx,
         {-78.8744999594492, -658.202587316166, 9.81380704582261},
         {78.8744999594492, -658.202587316166, 9.81380704582261}},
        {face_lx,
         {-646.627770852789, -364.397142609082, -436.461488890956},
         {646.627770852789, -364.397142609082, -436.461488890956}},
        {face_lx,
         {-610.560483423535, 762.27561485641, 552.200309352325},
         {610.560483423535, 762.27561485641, 552.200309352325}},
        {face_lx,
         {-612.375294788626, 193.387346448907, -845.955378730611},
         {612.375294788626, 193.387346448907, -845.955378730611}},
        {face_lx,
         {-564.724317207482, 335.038485673937, -245.265535953844},
         {564.724317207482, 335.038485673937, -245.265535953844}},
        {face_lx,
         {-910.296196719688, 699.542024728447, 97.8533024194962},
         {910.296196719688, 699.542024728447, 97.8533024194962}},
        {face_lx,
         {-430.3502841618, 63.8721038327781, -114.19917533304},
         {430.3502841618, 63.8721038327781, -114.19917533304}},
        {face_lx,
         {-368.252951096588, -665.104618425566, -665.178254669152},
         {368.252951096588, -665.104618425566, -665.178254669152}},
        {face_lx,
         {-642.04414460488, 935.027336815072, -9.30607608441869},
         {642.04414460488, 935.027336815072, -9.30607608441869}},
        {face_lx,
         {-871.233606240533, 752.027019796737, -115.800397726734},
         {871.233606240533, 752.027019796737, -115.800397726734}},
        {face_lx,
         {-727.809324463428, -394.654327368731, 505.13839311672},
         {727.809324463428, -394.654327368731, 505.13839311672}},
        {face_lx,
         {-478.624684077544, -499.970575373512, 578.008173437864},
         {478.624684077544, -499.970575373512, 578.008173437864}},
        {face_lx,
         {-160.946973463212, 88.8610551467186, 750.825335257757},
         {160.946973463212, 88.8610551467186, 750.825335257757}},
        {face_lx,
         {-664.253692220367, -136.767376443464, 428.751956994162},
         {664.253692220367, -136.767376443464, 428.751956994162}},
        {face_lx,
         {-657.383806162362, 437.703560394032, 702.864092876207},
         {657.383806162362, 437.703560394032, 702.864092876207}},
        {face_lx,
         {-250.912765797588, -360.555551367498, -417.935611182772},
         {250.912765797588, -360.555551367498, -417.935611182772}},
        {face_lx,
         {-676.144408903461, -286.185489466873, 679.186992785375},
         {676.144408903461, -286.185489466873, 679.186992785375}},
        {face_lx,
         {-265.018206494389, 740.225762151526, -992.398427924004},
         {265.018206494389, 740.225762151526, -992.398427924004}},
        {face_lx,
         {-258.000817764067, 25.9758167300834, 203.122941452231},
         {258.000817764067, 25.9758167300834, 203.122941452231}},
        {face_lx,
         {-589.160282006294, -832.378818413793, 628.9865269422},
         {589.160282006294, -832.378818413793, 628.9865269422}},
        {face_lx,
         {-186.852077780023, -862.892172191315, 420.208876251524},
         {186.852077780023, -862.892172191315, 420.208876251524}},
        {face_lx,
         {-151.750962535673, 124.92335014943, -946.522719648955},
         {151.750962535673, 124.92335014943, -946.522719648955}},
        {face_lx,
         {-126.925634383238, 463.957589817906, -333.992285311442},
         {126.925634383238, 463.957589817906, -333.992285311442}},
        {face_lx,
         {-630.741905292101, -158.246599300932, -886.992097662228},
         {630.741905292101, -158.246599300932, -886.992097662228}},
        {face_lx,
         {-532.497190078551, 37.6656242169552, -706.978855036085},
         {532.497190078551, 37.6656242169552, -706.978855036085}},
        {face_lx,
         {-975.1630952092, 142.05246474861, 620.574714967058},
         {975.1630952092, 142.05246474861, 620.574714967058}},
        {face_lx,
         {-867.462282435192, 220.05888539482, 249.922483485251},
         {867.462282435192, 220.05888539482, 249.922483485251}},
        {face_lx,
         {-472.025839204428, -345.565674525928, -275.781692852454},
         {472.025839204428, -345.565674525928, -275.781692852454}},
        {face_lx,
         {-438.404156728431, 267.990529484678, 494.134336500076},
         {438.404156728431, 267.990529484678, 494.134336500076}},
        {face_lx,
         {-147.41135366174, 207.293044203463, 811.556578012792},
         {147.41135366174, 207.293044203463, 811.556578012792}},
        {face_lx,
         {-504.644691499489, -224.41111153912, -549.382011197681},
         {504.644691499489, -224.41111153912, -549.382011197681}},
        {face_lx,
         {-421.563555206835, 212.848375732636, 212.719576980951},
         {421.563555206835, 212.848375732636, 212.719576980951}},
        {face_lx,
         {-15.8237769428783, -353.107512574303, -483.916985301113},
         {15.8237769428783, -353.107512574303, -483.916985301113}},
        {face_lx,
         {-916.571317880852, -790.376520403926, 148.458233217024},
         {916.571317880852, -790.376520403926, 148.458233217024}},
        {face_lx,
         {-697.674879858158, 937.200077446572, -607.16963979874},
         {697.674879858158, 937.200077446572, -607.16963979874}},
        {face_lx,
         {-347.283662982524, 249.427842377049, -660.764755470114},
         {347.283662982524, 249.427842377049, -660.764755470114}},
        {face_lx,
         {-786.836651198861, 244.870688402206, 717.066141550437},
         {786.836651198861, 244.870688402206, 717.066141550437}},
        {face_ly,
         {531.883949299715, -266.017155097828, -744.51849050471},
         {531.883949299715, 266.017155097828, -744.51849050471}},
        {face_ly,
         {-527.268669236957, -665.102163575284, -584.178500433729},
         {-527.268669236957, 665.102163575284, -584.178500433729}},
        {face_ly,
         {103.176226774077, -173.223585684261, -941.760604413665},
         {103.176226774077, 173.223585684261, -941.760604413665}},
        {face_ly,
         {-395.924501499104, -358.50511350722, -851.143303547092},
         {-395.924501499104, 358.50511350722, -851.143303547092}},
        {face_ly,
         {330.735636461778, -827.334889012219, 37.3669351096924},
         {330.735636461778, 827.334889012219, 37.3669351096924}},
        {face_ly,
         {990.02640595657, -125.092024187134, -583.671663020905},
         {990.02640595657, 125.092024187134, -583.671663020905}},
        {face_ly,
         {465.08605471716, -761.069515185445, -500.101490625982},
         {465.08605471716, 761.069515185445, -500.101490625982}},
        {face_ly,
         {209.213046917553, -644.107630519621, 293.226961074639},
         {209.213046917553, 644.107630519621, 293.226961074639}},
        {face_ly,
         {374.263588831724, -847.341259106852, 148.62647707355},
         {374.263588831724, 847.341259106852, 148.62647707355}},
        {face_ly,
         {-46.6176732352437, -237.294221076207, 670.143387938975},
         {-46.6176732352437, 237.294221076207, 670.143387938975}},
        {face_ly,
         {-883.437097633293, -766.540822387572, 46.4567036473109},
         {-883.437097633293, 766.540822387572, 46.4567036473109}},
        {face_ly,
         {873.076669014642, -249.054735805428, -311.616884812153},
         {873.076669014642, 249.054735805428, -311.616884812153}},
        {face_ly,
         {349.499335754416, -57.5447556448034, 613.497770476909},
         {349.499335754416, 57.5447556448034, 613.497770476909}},
        {face_ly,
         {13.4152355077922, -904.000140210949, -710.304703972689},
         {13.4152355077922, 904.000140210949, -710.304703972689}},
        {face_ly,
         {476.825346987395, -763.155848557012, 106.914318723484},
         {476.825346987395, 763.155848557012, 106.914318723484}},
        {face_ly,
         {595.238584180804, -57.8148323744699, -308.250076963946},
         {595.238584180804, 57.8148323744699, -308.250076963946}},
        {face_ly,
         {-819.902951414918, -16.828400452036, 281.649800633024},
         {-819.902951414918, 16.828400452036, 281.649800633024}},
        {face_ly,
         {-30.1483885865878, -952.683977578885, -265.024689308035},
         {-30.1483885865878, 952.683977578885, -265.024689308035}},
        {face_ly,
         {-115.854634895039, -750.253866388227, 342.855274122445},
         {-115.854634895039, 750.253866388227, 342.855274122445}},
        {face_ly,
         {-78.3817955695454, -494.299723861835, 16.3191622951235},
         {-78.3817955695454, 494.299723861835, 16.3191622951235}},
        {face_ly,
         {982.178128267054, -817.256147698321, -655.129356785449},
         {982.178128267054, 817.256147698321, -655.129356785449}},
        {face_ly,
         {-307.571241088121, -934.131404877454, 672.482392572087},
         {-307.571241088121, 934.131404877454, 672.482392572087}},
        {face_ly,
         {-251.671766030614, -793.67298360326, 329.561657447575},
         {-251.671766030614, 793.67298360326, 329.561657447575}},
        {face_ly,
         {175.408327352356, -387.495216770701, 470.297016331152},
         {175.408327352356, 387.495216770701, 470.297016331152}},
        {face_ly,
         {629.978541034191, -93.7201800931143, -776.774336271969},
         {629.978541034191, 93.7201800931143, -776.774336271969}},
        {face_ly,
         {-136.895727222468, -148.885288831187, 678.644387178556},
         {-136.895727222468, 148.885288831187, 678.644387178556}},
        {face_ly,
         {955.165502792529, -168.724856488364, -149.635544506869},
         {955.165502792529, 168.724856488364, -149.635544506869}},
        {face_ly,
         {-210.736196426267, -835.992960194079, 141.541907262025},
         {-210.736196426267, 835.992960194079, 141.541907262025}},
        {face_ly,
         {-236.678351872676, -676.816426051103, -188.529985540615},
         {-236.678351872676, 676.816426051103, -188.529985540615}},
        {face_ly,
         {-49.5119283940439, -828.108591863812, 418.687262453994},
         {-49.5119283940439, 828.108591863812, 418.687262453994}},
        {face_ly,
         {-694.782332926021, -519.575915506875, 735.320332010296},
         {-694.782332926021, 519.575915506875, 735.320332010296}},
        {face_ly,
         {-649.765537420595, -546.676662780783, 859.373757600433},
         {-649.765537420595, 546.676662780783, 859.373757600433}},
        {face_ly,
         {-268.997077258267, -192.091864557849, -69.0361474655147},
         {-268.997077258267, 192.091864557849, -69.0361474655147}},
        {face_ly,
         {264.53762739621, -394.050969568656, -521.487870980175},
         {264.53762739621, 394.050969568656, -521.487870980175}},
        {face_ly,
         {-491.72464557124, -249.949401730752, -370.384319264205},
         {-491.72464557124, 249.949401730752, -370.384319264205}},
        {face_ly,
         {-920.487015241907, -950.316753714247, -836.435124621105},
         {-920.487015241907, 950.316753714247, -836.435124621105}},
        {face_ly,
         {-490.983624887822, -186.962299473375, -8.1143872285611},
         {-490.983624887822, 186.962299473375, -8.1143872285611}},
        {face_ly,
         {208.400649564176, -297.381315433441, -862.655659487997},
         {208.400649564176, 297.381315433441, -862.655659487997}},
        {face_ly,
         {317.640836783701, -376.607921505737, -522.18622830775},
         {317.640836783701, 376.607921505737, -522.18622830775}},
        {face_ly,
         {575.380428572348, -647.493928137712, 805.620228537503},
         {575.380428572348, 647.493928137712, 805.620228537503}},
        {face_ly,
         {-326.888555484475, -316.898957992133, 192.025233982586},
         {-326.888555484475, 316.898957992133, 192.025233982586}},
        {face_ly,
         {-422.707342135244, -858.123687785639, -58.21407132797},
         {-422.707342135244, 858.123687785639, -58.21407132797}},
        {face_ly,
         {-931.920796789307, -743.257376269306, -854.092195927493},
         {-931.920796789307, 743.257376269306, -854.092195927493}},
        {face_ly,
         {-680.519907235183, -184.723760432758, -379.569753276762},
         {-680.519907235183, 184.723760432758, -379.569753276762}},
        {face_ly,
         {-420.332039122782, -17.3563687927826, -922.017061591859},
         {-420.332039122782, 17.3563687927826, -922.017061591859}},
        {face_ly,
         {35.301517149715, -168.990312703359, -488.221472970833},
         {35.301517149715, 168.990312703359, -488.221472970833}},
        {face_ly,
         {-845.065213743073, -453.335211181681, -766.983068662858},
         {-845.065213743073, 453.335211181681, -766.983068662858}},
        {face_ly,
         {605.55862771097, -356.1295212629, -851.365761299805},
         {605.55862771097, 356.1295212629, -851.365761299805}},
        {face_ly,
         {147.332717797885, -572.420763173059, 180.529504177459},
         {147.332717797885, 572.420763173059, 180.529504177459}},
        {face_ly,
         {654.568056683135, -175.014801441182, 101.602393687871},
         {654.568056683135, 175.014801441182, 101.602393687871}},
        {face_ly,
         {923.459362393562, -339.049466103626, 19.0116036881295},
         {923.459362393562, 339.049466103626, 19.0116036881295}},
        {face_ly,
         {952.209146102985, -957.613156600256, -189.940112135924},
         {952.209146102985, 957.613156600256, -189.940112135924}},
        {face_ly,
         {-704.902581113246, -316.01376289761, -185.037475786392},
         {-704.902581113246, 316.01376289761, -185.037475786392}},
        {face_ly,
         {833.808803670295, -34.8193676625888, -58.7185745852885},
         {833.808803670295, 34.8193676625888, -58.7185745852885}},
        {face_ly,
         {800.072267418056, -809.304543811676, 929.496913853168},
         {800.072267418056, 809.304543811676, 929.496913853168}},
        {face_ly,
         {859.270559327888, -127.614811397769, 779.967679340686},
         {859.270559327888, 127.614811397769, 779.967679340686}},
        {face_ly,
         {233.556613734479, -849.00425745839, -330.978002257606},
         {233.556613734479, 849.00425745839, -330.978002257606}},
        {face_ly,
         {-855.827194293201, -203.148732101879, -793.396910597819},
         {-855.827194293201, 203.148732101879, -793.396910597819}},
        {face_ly,
         {-251.261517417701, -839.974437108952, -423.276738213018},
         {-251.261517417701, 839.974437108952, -423.276738213018}},
        {face_ly,
         {732.655316994181, -702.21618521761, -858.275715620337},
         {732.655316994181, 702.21618521761, -858.275715620337}},
        {face_ly,
         {-836.154699516963, -26.0610812707364, -303.08869817201},
         {-836.154699516963, 26.0610812707364, -303.08869817201}},
        {face_ly,
         {-97.2539685819688, -229.27132394929, 586.54565937723},
         {-97.2539685819688, 229.27132394929, 586.54565937723}},
        {face_ly,
         {-527.479919050704, -673.418860010487, 392.95256503314},
         {-527.479919050704, 673.418860010487, 392.95256503314}},
        {face_ly,
         {-1.56091006375527, -729.716319096231, 417.252353796346},
         {-1.56091006375527, 729.716319096231, 417.252353796346}},
        {face_ly,
         {908.862634288762, -47.1117455184194, -186.400036711578},
         {908.862634288762, 47.1117455184194, -186.400036711578}},
        {face_ly,
         {-448.940024783518, -860.866994581676, -643.775046391205},
         {-448.940024783518, 860.866994581676, -643.775046391205}},
        {face_ly,
         {398.754534517446, -658.935371576582, -16.0455697735965},
         {398.754534517446, 658.935371576582, -16.0455697735965}},
        {face_ly,
         {-892.017461356709, -793.144598819467, 738.662366546737},
         {-892.017461356709, 793.144598819467, 738.662366546737}},
        {face_ly,
         {356.599254162338, -336.347392829412, 148.692126503203},
         {356.599254162338, 336.347392829412, 148.692126503203}},
        {face_ly,
         {478.894261552134, -143.637281339942, 264.660162261928},
         {478.894261552134, 143.637281339942, 264.660162261928}},
        {face_ly,
         {-256.736912352354, -833.044149102499, -518.562173267322},
         {-256.736912352354, 833.044149102499, -518.562173267322}},
        {face_ly,
         {-894.725134847385, -603.663892978174, 284.4646019397},
         {-894.725134847385, 603.663892978174, 284.4646019397}},
        {face_ly,
         {826.301797910485, -994.53965816299, 867.996049181071},
         {826.301797910485, 994.53965816299, 867.996049181071}},
        {face_ly,
         {-525.645632086223, -494.333573440135, 731.55203107972},
         {-525.645632086223, 494.333573440135, 731.55203107972}},
        {face_ly,
         {-823.091470047534, -572.72686055026, -162.01167565871},
         {-823.091470047534, 572.72686055026, -162.01167565871}},
        {face_ly,
         {-508.930368137079, -890.126891328408, -750.258267469361},
         {-508.930368137079, 890.126891328408, -750.258267469361}},
        {face_ly,
         {-375.497272069215, -444.26478140417, -843.829344834326},
         {-375.497272069215, 444.26478140417, -843.829344834326}},
        {face_ly,
         {-734.150610656854, -738.18527928116, -612.716317043056},
         {-734.150610656854, 738.18527928116, -612.716317043056}},
        {face_ly,
         {-482.556643522326, -970.141445162975, -735.131619927659},
         {-482.556643522326, 970.141445162975, -735.131619927659}},
        {face_ly,
         {338.598552345885, -679.736409148149, -186.043134388298},
         {338.598552345885, 679.736409148149, -186.043134388298}},
        {face_ly,
         {-258.383266955813, -612.017561400374, -650.814065440379},
         {-258.383266955813, 612.017561400374, -650.814065440379}},
        {face_ly,
         {957.109840201644, -979.310592576725, -998.37825624869},
         {957.109840201644, 979.310592576725, -998.37825624869}},
        {face_ly,
         {952.595241963886, -267.761808799044, 497.964654531279},
         {952.595241963886, 267.761808799044, 497.964654531279}},
        {face_ly,
         {-693.67995869011, -51.2052216016291, 760.230094340523},
         {-693.67995869011, 51.2052216016291, 760.230094340523}},
        {face_ly,
         {-66.2623979869131, -971.002355556833, 939.112315756803},
         {-66.2623979869131, 971.002355556833, 939.112315756803}},
        {face_ly,
         {115.358553786405, -327.115181331631, 509.80613167158},
         {115.358553786405, 327.115181331631, 509.80613167158}},
        {face_ly,
         {424.535509520502, -28.1693109950174, 14.9134205043624},
         {424.535509520502, 28.1693109950174, 14.9134205043624}},
        {face_ly,
         {181.609135399759, -38.3470741069496, 755.048167251248},
         {181.609135399759, 38.3470741069496, 755.048167251248}},
        {face_ly,
         {482.179714694972, -855.119641253141, 285.06091287173},
         {482.179714694972, 855.119641253141, 285.06091287173}},
        {face_ly,
         {-263.878235858462, -358.214660563999, 211.954231292926},
         {-263.878235858462, 358.214660563999, 211.954231292926}},
        {face_ly,
         {-819.495488592471, -202.51734917038, 497.631281324389},
         {-819.495488592471, 202.51734917038, 497.631281324389}},
        {face_ly,
         {-961.944976859819, -885.025322758706, 207.802891690943},
         {-961.944976859819, 885.025322758706, 207.802891690943}},
        {face_ly,
         {-28.0976721173824, -800.173550791591, 318.400620757676},
         {-28.0976721173824, 800.173550791591, 318.400620757676}},
        {face_ly,
         {-795.194476796941, -623.226960408269, 601.653163170672},
         {-795.194476796941, 623.226960408269, 601.653163170672}},
        {face_ly,
         {-781.720383249695, -566.59064813211, 711.441080612667},
         {-781.720383249695, 566.59064813211, 711.441080612667}},
        {face_ly,
         {-858.900293214889, -74.9405767491253, 61.1818729902993},
         {-858.900293214889, 74.9405767491253, 61.1818729902993}},
        {face_ly,
         {929.91839514506, -321.13229060805, 832.087343462797},
         {929.91839514506, 321.13229060805, 832.087343462797}},
        {face_ly,
         {-428.070924956677, -309.843674601644, 532.988257623776},
         {-428.070924956677, 309.843674601644, 532.988257623776}},
        {face_ly,
         {-245.345724979067, -206.87201063254, 655.16252756806},
         {-245.345724979067, 206.87201063254, 655.16252756806}},
        {face_ly,
         {-344.609249962145, -492.39694323798, -123.688567264765},
         {-344.609249962145, 492.39694323798, -123.688567264765}},
        {face_lz,
         {-541.169815812339, 903.826701688362, -468.596154947391},
         {-541.169815812339, 903.826701688362, 468.596154947391}},
        {face_lz,
         {-852.251766499304, 363.745307163111, -721.362413098998},
         {-852.251766499304, 363.745307163111, 721.362413098998}},
        {face_lz,
         {-76.1633556388715, -833.423564101907, -247.3498629329},
         {-76.1633556388715, -833.423564101907, 247.3498629329}},
        {face_lz,
         {-337.068059075234, -822.251274263951, -289.811115009995},
         {-337.068059075234, -822.251274263951, 289.811115009995}},
        {face_lz,
         {339.645117243127, 965.699537331483, -978.306067930705},
         {339.645117243127, 965.699537331483, 978.306067930705}},
        {face_lz,
         {-860.108051023429, -16.4926000867199, -68.0675780668562},
         {-860.108051023429, -16.4926000867199, 68.0675780668562}},
        {face_lz,
         {584.452524784821, -330.290686616053, -805.809427102203},
         {584.452524784821, -330.290686616053, 805.809427102203}},
        {face_lz,
         {360.965679343121, 815.181189950489, -997.796083329012},
         {360.965679343121, 815.181189950489, 997.796083329012}},
        {face_lz,
         {469.830222747411, -739.166377235962, -967.485659083973},
         {469.830222747411, -739.166377235962, 967.485659083973}},
        {face_lz,
         {-279.302861908084, -209.997802782391, -27.6113969882363},
         {-279.302861908084, -209.997802782391, 27.6113969882363}},
        {face_lz,
         {-914.536113203683, -666.860799223259, -787.678302117195},
         {-914.536113203683, -666.860799223259, 787.678302117195}},
        {face_lz,
         {-879.426076769678, 690.840394758031, -285.134448192837},
         {-879.426076769678, 690.840394758031, 285.134448192837}},
        {face_lz,
         {-979.367340860899, -599.108246980694, -451.419166199693},
         {-979.367340860899, -599.108246980694, 451.419166199693}},
        {face_lz,
         {-128.306579992169, 430.812569192646, -91.07625021388},
         {-128.306579992169, 430.812569192646, 91.07625021388}},
        {face_lz,
         {-593.817611961379, 621.49608362216, -503.412195598011},
         {-593.817611961379, 621.49608362216, 503.412195598011}},
        {face_lz,
         {754.3906198587, -142.060655010598, -452.777958187479},
         {754.3906198587, -142.060655010598, 452.777958187479}},
        {face_lz,
         {-399.306850581736, 700.225792421181, -871.429075109207},
         {-399.306850581736, 700.225792421181, 871.429075109207}},
        {face_lz,
         {-358.006381221761, 293.707323022757, -431.176769539842},
         {-358.006381221761, 293.707323022757, 431.176769539842}},
        {face_lz,
         {734.619061714059, -851.401347664291, -241.386007593484},
         {734.619061714059, -851.401347664291, 241.386007593484}},
        {face_lz,
         {888.619617575903, -747.613368699807, -320.995161580218},
         {888.619617575903, -747.613368699807, 320.995161580218}},
        {face_lz,
         {485.311632537711, -721.67833478166, -826.927001722639},
         {485.311632537711, -721.67833478166, 826.927001722639}},
        {face_lz,
         {250.407580248042, 72.2794854972476, -798.141683677778},
         {250.407580248042, 72.2794854972476, 798.141683677778}},
        {face_lz,
         {-659.445203411142, -831.079622336741, -565.959179637154},
         {-659.445203411142, -831.079622336741, 565.959179637154}},
        {face_lz,
         {-838.529790499333, 231.893878032759, -325.267023042959},
         {-838.529790499333, 231.893878032759, 325.267023042959}},
        {face_lz,
         {-577.828459060733, 807.227163270576, -860.949304872046},
         {-577.828459060733, 807.227163270576, 860.949304872046}},
        {face_lz,
         {-321.342118298074, 934.926639352004, -820.960617810757},
         {-321.342118298074, 934.926639352004, 820.960617810757}},
        {face_lz,
         {-380.083164894902, 136.680095536208, -575.882594017857},
         {-380.083164894902, 136.680095536208, 575.882594017857}},
        {face_lz,
         {-16.0404721432392, 183.183089523692, -528.576920455454},
         {-16.0404721432392, 183.183089523692, 528.576920455454}},
        {face_lz,
         {-915.312167306952, -744.642208452078, -766.526533301705},
         {-915.312167306952, -744.642208452078, 766.526533301705}},
        {face_lz,
         {269.10417643328, 614.41906577504, -893.006383934264},
         {269.10417643328, 614.41906577504, 893.006383934264}},
        {face_lz,
         {677.987282023019, -145.062078027395, -980.198652080665},
         {677.987282023019, -145.062078027395, 980.198652080665}},
        {face_lz,
         {-186.039608614671, 649.66205813137, -719.150441982479},
         {-186.039608614671, 649.66205813137, 719.150441982479}},
        {face_lz,
         {212.29607605452, -772.38132748991, -62.185049814891},
         {212.29607605452, -772.38132748991, 62.185049814891}},
        {face_lz,
         {518.295878120685, 207.334190773519, -253.639122909441},
         {518.295878120685, 207.334190773519, 253.639122909441}},
        {face_lz,
         {297.504999912684, -497.743799415169, -349.078799262676},
         {297.504999912684, -497.743799415169, 349.078799262676}},
        {face_lz,
         {-339.517085586305, 548.554979662962, -197.740165200124},
         {-339.517085586305, 548.554979662962, 197.740165200124}},
        {face_lz,
         {-554.653573294661, 889.715262256328, -645.240652440544},
         {-554.653573294661, 889.715262256328, 645.240652440544}},
        {face_lz,
         {18.2506477264956, -837.133483759732, -413.756550832985},
         {18.2506477264956, -837.133483759732, 413.756550832985}},
        {face_lz,
         {-764.582991958129, 47.4235813248051, -256.015561889033},
         {-764.582991958129, 47.4235813248051, 256.015561889033}},
        {face_lz,
         {319.466585938606, -250.578146712892, -742.794704843222},
         {319.466585938606, -250.578146712892, 742.794704843222}},
        {face_lz,
         {93.9735186533176, 877.463845348412, -196.158284910898},
         {93.9735186533176, 877.463845348412, 196.158284910898}},
        {face_lz,
         {-633.542116162362, -661.779091010525, -974.582969458424},
         {-633.542116162362, -661.779091010525, 974.582969458424}},
        {face_lz,
         {968.817154969056, -259.901840646367, -401.52677586979},
         {968.817154969056, -259.901840646367, 401.52677586979}},
        {face_lz,
         {-888.635004042528, -170.318071299306, -586.871442555846},
         {-888.635004042528, -170.318071299306, 586.871442555846}},
        {face_lz,
         {-957.006732596171, 813.504912207231, -917.676649238252},
         {-957.006732596171, 813.504912207231, 917.676649238252}},
        {face_lz,
         {-6.94293624457441, 659.314000658029, -481.781418495995},
         {-6.94293624457441, 659.314000658029, 481.781418495995}},
        {face_lz,
         {229.122972377159, -959.560488070688, -834.076629206174},
         {229.122972377159, -959.560488070688, 834.076629206174}},
        {face_lz,
         {284.34448411595, -291.452815230156, -960.145856608428},
         {284.34448411595, -291.452815230156, 960.145856608428}},
        {face_lz,
         {671.354212833919, 150.85352002718, -986.999565942957},
         {671.354212833919, 150.85352002718, 986.999565942957}},
        {face_lz,
         {-820.861438145815, 45.850647575352, -921.739727313232},
         {-820.861438145815, 45.850647575352, 921.739727313232}},
        {face_lz,
         {-141.706417862516, -456.529268614594, -379.742682212245},
         {-141.706417862516, -456.529268614594, 379.742682212245}},
        {face_lz,
         {266.199177182284, 140.160330912304, -364.089268201184},
         {266.199177182284, 140.160330912304, 364.089268201184}},
        {face_lz,
         {118.824165271699, -423.247076449491, -567.668859762129},
         {118.824165271699, -423.247076449491, 567.668859762129}},
        {face_lz,
         {-383.836204622205, 755.418537745317, -962.373529522459},
         {-383.836204622205, 755.418537745317, 962.373529522459}},
        {face_lz,
         {-162.978221729553, -770.209136809513, -963.325472314353},
         {-162.978221729553, -770.209136809513, 963.325472314353}},
        {face_lz,
         {768.753232824892, 252.436853162838, -534.431183869742},
         {768.753232824892, 252.436853162838, 534.431183869742}},
        {face_lz,
         {488.885655482375, 779.751786369525, -577.6732894396},
         {488.885655482375, 779.751786369525, 577.6732894396}},
        {face_lz,
         {228.895415526839, 189.246657810039, -752.254674458506},
         {228.895415526839, 189.246657810039, 752.254674458506}},
        {face_lz,
         {662.026769059421, -52.996970042771, -314.086247771769},
         {662.026769059421, -52.996970042771, 314.086247771769}},
        {face_lz,
         {-827.417226662719, -789.597469195893, -440.650307813129},
         {-827.417226662719, -789.597469195893, 440.650307813129}},
        {face_lz,
         {-923.260241996761, -293.924021799632, -559.559976983937},
         {-923.260241996761, -293.924021799632, 559.559976983937}},
        {face_lz,
         {687.993764601063, -708.600250252383, -289.098797604787},
         {687.993764601063, -708.600250252383, 289.098797604787}},
        {face_lz,
         {-995.685034980414, 432.272715623884, -522.114033664947},
         {-995.685034980414, 432.272715623884, 522.114033664947}},
        {face_lz,
         {-588.064306556396, -806.039362825604, -292.510762381438},
         {-588.064306556396, -806.039362825604, 292.510762381438}},
        {face_lz,
         {-880.163148597178, -90.4295461219203, -556.571303500395},
         {-880.163148597178, -90.4295461219203, 556.571303500395}},
        {face_lz,
         {-428.724080216792, -238.929535366381, -594.762940807212},
         {-428.724080216792, -238.929535366381, 594.762940807212}},
        {face_lz,
         {724.429194864137, 177.841648762991, -732.537249866271},
         {724.429194864137, 177.841648762991, 732.537249866271}},
        {face_lz,
         {-423.693607734345, -742.280633020631, -810.333258503778},
         {-423.693607734345, -742.280633020631, 810.333258503778}},
        {face_lz,
         {489.779766415749, -151.965803537837, -758.796838654634},
         {489.779766415749, -151.965803537837, 758.796838654634}},
        {face_lz,
         {394.403356253045, 990.040617238218, -605.570597518745},
         {394.403356253045, 990.040617238218, 605.570597518745}},
        {face_lz,
         {-102.529123006719, 177.459122300329, -786.169809263608},
         {-102.529123006719, 177.459122300329, 786.169809263608}},
        {face_lz,
         {849.160382116226, -462.141073563248, -28.1932499754157},
         {849.160382116226, -462.141073563248, 28.1932499754157}},
        {face_lz,
         {-30.2223322632008, -585.473944638256, -986.273633539698},
         {-30.2223322632008, -585.473944638256, 986.273633539698}},
        {face_lz,
         {703.481738263394, 150.297108224577, -819.294285240904},
         {703.481738263394, 150.297108224577, 819.294285240904}},
        {face_lz,
         {798.356743931272, 29.0867730719274, -128.1315336698},
         {798.356743931272, 29.0867730719274, 128.1315336698}},
        {face_lz,
         {-819.43936717721, 351.506239838921, -681.553572341051},
         {-819.43936717721, 351.506239838921, 681.553572341051}},
        {face_lz,
         {333.449385954789, 776.005312102075, -457.472181940357},
         {333.449385954789, 776.005312102075, 457.472181940357}},
        {face_lz,
         {296.838060585644, -233.686013900276, -728.522702840544},
         {296.838060585644, -233.686013900276, 728.522702840544}},
        {face_lz,
         {178.410551944048, -12.7657696221463, -40.3299054601439},
         {178.410551944048, -12.7657696221463, 40.3299054601439}},
        {face_lz,
         {-503.120138961942, -586.028289553155, -791.585668969363},
         {-503.120138961942, -586.028289553155, 791.585668969363}},
        {face_lz,
         {-244.43587364087, 872.124028051038, -101.956233818989},
         {-244.43587364087, 872.124028051038, 101.956233818989}},
        {face_lz,
         {2.85142856433185, 744.620973901498, -704.229201580584},
         {2.85142856433185, 744.620973901498, 704.229201580584}},
        {face_lz,
         {-58.8462442576201, 370.152871739551, -342.698066555073},
         {-58.8462442576201, 370.152871739551, 342.698066555073}},
        {face_lz,
         {-512.598216747857, -376.126436637821, -523.674726207461},
         {-512.598216747857, -376.126436637821, 523.674726207461}},
        {face_lz,
         {498.286444846884, -66.4101710431478, -342.399122362578},
         {498.286444846884, -66.4101710431478, 342.399122362578}},
        {face_lz,
         {12.4017273608115, 619.738564738088, -631.242533314011},
         {12.4017273608115, 619.738564738088, 631.242533314011}},
        {face_lz,
         {313.616570718732, 441.279640715378, -10.4227533381668},
         {313.616570718732, 441.279640715378, 10.4227533381668}},
        {face_lz,
         {805.966735665594, -312.196319621743, -572.700472477812},
         {805.966735665594, -312.196319621743, 572.700472477812}},
        {face_lz,
         {649.426864423126, 275.912571638448, -380.391634834924},
         {649.426864423126, 275.912571638448, 380.391634834924}},
        {face_lz,
         {1.98798760729596, -501.056170100848, -220.748132096945},
         {1.98798760729596, -501.056170100848, 220.748132096945}},
        {face_lz,
         {101.617886966098, 292.807260111882, -199.114818375521},
         {101.617886966098, 292.807260111882, 199.114818375521}},
        {face_lz,
         {331.822008635026, 293.498412485702, -24.222318484598},
         {331.822008635026, 293.498412485702, 24.222318484598}},
        {face_lz,
         {-445.505481897747, -257.553323568728, -377.491745134221},
         {-445.505481897747, -257.553323568728, 377.491745134221}},
        {face_lz,
         {920.636135123735, 38.1816693435803, -710.842354957914},
         {920.636135123735, 38.1816693435803, 710.842354957914}},
        {face_lz,
         {-983.06546685655, -183.846030894833, -53.0725989582788},
         {-983.06546685655, -183.846030894833, 53.0725989582788}},
        {face_lz,
         {-240.650218713057, -177.048847052517, -327.085952399198},
         {-240.650218713057, -177.048847052517, 327.085952399198}},
        {face_lz,
         {-600.068538366976, 892.911237503561, -596.034646012221},
         {-600.068538366976, 892.911237503561, 596.034646012221}},
        {face_lz,
         {-767.97742142427, 813.718021606683, -351.561534683004},
         {-767.97742142427, 813.718021606683, 351.561534683004}},
        {face_lz,
         {69.6397987675382, 551.596288816513, -270.41719082155},
         {69.6397987675382, 551.596288816513, 270.41719082155}},
        {face_lz,
         {974.384657253937, -405.079938265978, -185.501098233034},
         {974.384657253937, -405.079938265978, 185.501098233034}},
        {face_hx,
         {59.5091489597517, 203.754150930265, 261.94901298058},
         {-59.5091489597517, 203.754150930265, 261.94901298058}},
        {face_hx,
         {207.772285303202, 851.568836921269, -738.696486185456},
         {-207.772285303202, 851.568836921269, -738.696486185456}},
        {face_hx,
         {165.757201503146, 894.05416894066, 862.462028679541},
         {-165.757201503146, 894.05416894066, 862.462028679541}},
        {face_hx,
         {696.462895475509, -239.605163607494, 432.850787625788},
         {-696.462895475509, -239.605163607494, 432.850787625788}},
        {face_hx,
         {399.148233193929, -873.983421380046, 496.305984896326},
         {-399.148233193929, -873.983421380046, 496.305984896326}},
        {face_hx,
         {659.579538824482, 825.915318986827, 245.242568357006},
         {-659.579538824482, 825.915318986827, 245.242568357006}},
        {face_hx,
         {426.507374987304, 858.218325917909, 845.188757576926},
         {-426.507374987304, 858.218325917909, 845.188757576926}},
        {face_hx,
         {234.837354156449, 153.863019926311, 302.473824429332},
         {-234.837354156449, 153.863019926311, 302.473824429332}},
        {face_hx,
         {299.68612249762, 834.670348950352, -368.64277807377},
         {-299.68612249762, 834.670348950352, -368.64277807377}},
        {face_hx,
         {907.422771921906, 122.60079603137, -978.489512130402},
         {-907.422771921906, 122.60079603137, -978.489512130402}},
        {face_hx,
         {509.529817917221, 936.773687744686, -578.041757908019},
         {-509.529817917221, 936.773687744686, -578.041757908019}},
        {face_hx,
         {43.8417132807208, -890.830127504604, 514.576088401534},
         {-43.8417132807208, -890.830127504604, 514.576088401534}},
        {face_hx,
         {687.656123222235, 788.107891177032, 637.150356270978},
         {-687.656123222235, 788.107891177032, 637.150356270978}},
        {face_hx,
         {958.92646105781, -206.759244578512, -78.5866001478257},
         {-958.92646105781, -206.759244578512, -78.5866001478257}},
        {face_hx,
         {271.081520039892, -453.390728138416, 822.658264436795},
         {-271.081520039892, -453.390728138416, 822.658264436795}},
        {face_hx,
         {633.802275779455, -912.456380816364, 852.3083509826},
         {-633.802275779455, -912.456380816364, 852.3083509826}},
        {face_hx,
         {681.819815846465, 542.291370816351, -663.094303402379},
         {-681.819815846465, 542.291370816351, -663.094303402379}},
        {face_hx,
         {661.666214280487, 322.35249521093, -464.001665925384},
         {-661.666214280487, 322.35249521093, -464.001665925384}},
        {face_hx,
         {509.319370039633, -306.670546896157, 693.563624496919},
         {-509.319370039633, -306.670546896157, 693.563624496919}},
        {face_hx,
         {207.745263695107, -447.161709160069, -150.778552411705},
         {-207.745263695107, -447.161709160069, -150.778552411705}},
        {face_hx,
         {394.203985348858, -507.052429130302, -934.723889258949},
         {-394.203985348858, -507.052429130302, -934.723889258949}},
        {face_hx,
         {953.32256537525, -961.639939474588, 835.88380831924},
         {-953.32256537525, -961.639939474588, 835.88380831924}},
        {face_hx,
         {478.669120408924, -430.490387281109, 843.559245651588},
         {-478.669120408924, -430.490387281109, 843.559245651588}},
        {face_hx,
         {842.321598964266, -6.61023995378264, 937.523778383424},
         {-842.321598964266, -6.61023995378264, 937.523778383424}},
        {face_hx,
         {940.423240013093, 445.371774890387, 565.409033258748},
         {-940.423240013093, 445.371774890387, 565.409033258748}},
        {face_hx,
         {28.7288886691372, -248.538245175585, 844.022979899088},
         {-28.7288886691372, -248.538245175585, 844.022979899088}},
        {face_hx,
         {510.606777115009, 35.2123141031916, -277.021111458959},
         {-510.606777115009, 35.2123141031916, -277.021111458959}},
        {face_hx,
         {779.795614728176, -697.378187755218, 416.328172912972},
         {-779.795614728176, -697.378187755218, 416.328172912972}},
        {face_hx,
         {827.352252013493, 792.149823093867, -84.2794989275553},
         {-827.352252013493, 792.149823093867, -84.2794989275553}},
        {face_hx,
         {327.936537843117, 313.118701195391, 422.602621812999},
         {-327.936537843117, 313.118701195391, 422.602621812999}},
        {face_hx,
         {854.465908388501, 368.240210455139, -390.739504847724},
         {-854.465908388501, 368.240210455139, -390.739504847724}},
        {face_hx,
         {506.389042074794, 439.229498574179, 92.9905438200153},
         {-506.389042074794, 439.229498574179, 92.9905438200153}},
        {face_hx,
         {868.327787086929, -947.479702789549, 460.87090693548},
         {-868.327787086929, -947.479702789549, 460.87090693548}},
        {face_hx,
         {58.525216932399, 580.811493741098, -453.631738333274},
         {-58.525216932399, 580.811493741098, -453.631738333274}},
        {face_hx,
         {677.237420389217, 515.398891532732, -437.749399173519},
         {-677.237420389217, 515.398891532732, -437.749399173519}},
        {face_hx,
         {844.280726748051, 785.858930994331, 445.73141480835},
         {-844.280726748051, 785.858930994331, 445.73141480835}},
        {face_hx,
         {685.407354125322, 916.333096787576, 319.473173995239},
         {-685.407354125322, 916.333096787576, 319.473173995239}},
        {face_hx,
         {560.903531694373, 808.50559207891, 706.257061494863},
         {-560.903531694373, 808.50559207891, 706.257061494863}},
        {face_hx,
         {994.348897073096, 81.9142452449864, 976.122624270678},
         {-994.348897073096, 81.9142452449864, 976.122624270678}},
        {face_hx,
         {512.254962943644, -258.057238660894, -697.294855518488},
         {-512.254962943644, -258.057238660894, -697.294855518488}},
        {face_hx,
         {732.132019063143, 539.455937923885, 835.51263499569},
         {-732.132019063143, 539.455937923885, 835.51263499569}},
        {face_hx,
         {566.292973312935, 137.383529327263, 273.61662765632},
         {-566.292973312935, 137.383529327263, 273.61662765632}},
        {face_hx,
         {715.112239529223, -74.9933531859119, 540.322475688404},
         {-715.112239529223, -74.9933531859119, 540.322475688404}},
        {face_hx,
         {679.862860699611, 118.334733243783, 279.428418540024},
         {-679.862860699611, 118.334733243783, 279.428418540024}},
        {face_hx,
         {855.950193329858, -700.53945247694, 911.692674011259},
         {-855.950193329858, -700.53945247694, 911.692674011259}},
        {face_hx,
         {652.498761349615, 477.902881033522, -593.017909813015},
         {-652.498761349615, 477.902881033522, -593.017909813015}},
        {face_hx,
         {942.889663715555, 770.542747459513, 426.530548042059},
         {-942.889663715555, 770.542747459513, 426.530548042059}},
        {face_hx,
         {304.480082967823, 783.493447287228, 559.929266934161},
         {-304.480082967823, 783.493447287228, 559.929266934161}},
        {face_hx,
         {194.057010866543, -226.499170357452, -136.350137901409},
         {-194.057010866543, -226.499170357452, -136.350137901409}},
        {face_hx,
         {77.4187549997205, -974.42395407569, -462.462881275607},
         {-77.4187549997205, -974.42395407569, -462.462881275607}},
        {face_hx,
         {672.536328809018, -443.956635520849, 823.071509669268},
         {-672.536328809018, -443.956635520849, 823.071509669268}},
        {face_hx,
         {664.454501580804, -395.386112535346, -385.418078454251},
         {-664.454501580804, -395.386112535346, -385.418078454251}},
        {face_hx,
         {679.715804805717, -781.663944989566, 861.083962373494},
         {-679.715804805717, -781.663944989566, 861.083962373494}},
        {face_hx,
         {164.128870792429, 457.454926026458, -931.375551152472},
         {-164.128870792429, 457.454926026458, -931.375551152472}},
        {face_hx,
         {85.4294881877509, -49.0319248917504, 882.576273410043},
         {-85.4294881877509, -49.0319248917504, 882.576273410043}},
        {face_hx,
         {687.951184016585, -914.434657606541, 257.350606392037},
         {-687.951184016585, -914.434657606541, 257.350606392037}},
        {face_hx,
         {780.906803433797, 150.876156794295, 776.122527551226},
         {-780.906803433797, 150.876156794295, 776.122527551226}},
        {face_hx,
         {507.041064536743, -598.534759989983, 998.523464696281},
         {-507.041064536743, -598.534759989983, 998.523464696281}},
        {face_hx,
         {802.291448717949, -310.069968378694, 610.178767938645},
         {-802.291448717949, -310.069968378694, 610.178767938645}},
        {face_hx,
         {93.4553850385591, -451.705758077839, -240.535232680593},
         {-93.4553850385591, -451.705758077839, -240.535232680593}},
        {face_hx,
         {969.60125355073, 67.8146028305523, 875.319580074908},
         {-969.60125355073, 67.8146028305523, 875.319580074908}},
        {face_hx,
         {147.878921349108, 167.321974640171, -224.797486847384},
         {-147.878921349108, 167.321974640171, -224.797486847384}},
        {face_hx,
         {368.460871231548, -62.9551429145549, -651.80498927839},
         {-368.460871231548, -62.9551429145549, -651.80498927839}},
        {face_hx,
         {23.3799059722364, -637.953943611043, 274.791113878331},
         {-23.3799059722364, -637.953943611043, 274.791113878331}},
        {face_hx,
         {398.112222883406, 256.4629104262, 69.1430839711202},
         {-398.112222883406, 256.4629104262, 69.1430839711202}},
        {face_hx,
         {220.082356484796, -343.841518544893, 775.631555299899},
         {-220.082356484796, -343.841518544893, 775.631555299899}},
        {face_hx,
         {903.892456132789, -440.95172232256, -368.812269551347},
         {-903.892456132789, -440.95172232256, -368.812269551347}},
        {face_hx,
         {854.179960057391, 677.480362345861, -432.230754656468},
         {-854.179960057391, 677.480362345861, -432.230754656468}},
        {face_hx,
         {819.312738301499, 827.045117379108, -210.533999383376},
         {-819.312738301499, 827.045117379108, -210.533999383376}},
        {face_hx,
         {325.42979429321, -440.636015150621, 604.936814559965},
         {-325.42979429321, -440.636015150621, 604.936814559965}},
        {face_hx,
         {805.826915751136, -168.680211039852, 20.1321557009333},
         {-805.826915751136, -168.680211039852, 20.1321557009333}},
        {face_hx,
         {385.092210698838, 897.825512285122, -277.720150186102},
         {-385.092210698838, 897.825512285122, -277.720150186102}},
        {face_hx,
         {721.126783700654, -2.56638476886792, -295.975002426072},
         {-721.126783700654, -2.56638476886792, -295.975002426072}},
        {face_hx,
         {412.360948309504, -345.252993089406, -913.983035753015},
         {-412.360948309504, -345.252993089406, -913.983035753015}},
        {face_hx,
         {412.341101050916, 875.144100903529, 778.947810091549},
         {-412.341101050916, 875.144100903529, 778.947810091549}},
        {face_hx,
         {588.196246054859, -97.6108340422156, 748.645052146399},
         {-588.196246054859, -97.6108340422156, 748.645052146399}},
        {face_hx,
         {422.956479238439, 868.21594197111, -929.390487257382},
         {-422.956479238439, 868.21594197111, -929.390487257382}},
        {face_hx,
         {965.222596285673, -452.877151523075, 132.987362018778},
         {-965.222596285673, -452.877151523075, 132.987362018778}},
        {face_hx,
         {915.217621939752, 471.124500776444, 27.524597488904},
         {-915.217621939752, 471.124500776444, 27.524597488904}},
        {face_hx,
         {539.293562332896, -39.1129703928796, 401.434103172254},
         {-539.293562332896, -39.1129703928796, 401.434103172254}},
        {face_hx,
         {91.8880026376946, -532.623223440535, -539.080308825135},
         {-91.8880026376946, -532.623223440535, -539.080308825135}},
        {face_hx,
         {387.515088809229, -908.977335951784, -760.401643550559},
         {-387.515088809229, -908.977335951784, -760.401643550559}},
        {face_hx,
         {252.077349019163, 244.471999500518, -80.2684556003305},
         {-252.077349019163, 244.471999500518, -80.2684556003305}},
        {face_hx,
         {965.897458532425, -374.995299516564, 763.024641251348},
         {-965.897458532425, -374.995299516564, 763.024641251348}},
        {face_hx,
         {875.80136399273, -566.866414754094, -92.1168933113295},
         {-875.80136399273, -566.866414754094, -92.1168933113295}},
        {face_hx,
         {372.055691670096, -480.710968246511, 341.853659938526},
         {-372.055691670096, -480.710968246511, 341.853659938526}},
        {face_hx,
         {194.360880838675, -878.810106932839, 54.5758539665094},
         {-194.360880838675, -878.810106932839, 54.5758539665094}},
        {face_hx,
         {342.24360282735, -925.041472137937, -664.923369732842},
         {-342.24360282735, -925.041472137937, -664.923369732842}},
        {face_hx,
         {569.248351386235, -854.182374169255, -487.613452954384},
         {-569.248351386235, -854.182374169255, -487.613452954384}},
        {face_hx,
         {968.162695805561, 500.655774280711, 609.468150715802},
         {-968.162695805561, 500.655774280711, 609.468150715802}},
        {face_hx,
         {91.818215036747, 486.570784389461, 53.5818609645298},
         {-91.818215036747, 486.570784389461, 53.5818609645298}},
        {face_hx,
         {566.792840275417, 996.325394653806, 741.536378259008},
         {-566.792840275417, 996.325394653806, 741.536378259008}},
        {face_hx,
         {354.013961300394, 193.742742664374, -901.226962577714},
         {-354.013961300394, 193.742742664374, -901.226962577714}},
        {face_hx,
         {300.922024862231, 113.791431365364, 907.301044318104},
         {-300.922024862231, 113.791431365364, 907.301044318104}},
        {face_hx,
         {392.328682818677, -485.5303114738, -805.884058868823},
         {-392.328682818677, -485.5303114738, -805.884058868823}},
        {face_hx,
         {83.5965628079653, -344.995271427667, 563.068747667007},
         {-83.5965628079653, -344.995271427667, 563.068747667007}},
        {face_hx,
         {719.526886726092, 18.1778471403618, -78.2194638229466},
         {-719.526886726092, 18.1778471403618, -78.2194638229466}},
        {face_hx,
         {74.5608867150118, -326.064266393733, 530.418959515008},
         {-74.5608867150118, -326.064266393733, 530.418959515008}},
        {face_hx,
         {321.264781362128, 234.913214906179, 983.262297334954},
         {-321.264781362128, 234.913214906179, 983.262297334954}},
        {face_hx,
         {329.793894617079, 582.570708371208, 434.677279304304},
         {-329.793894617079, 582.570708371208, 434.677279304304}},
        {face_hy,
         {911.581064518776, 392.759470947783, 4.70433212259104},
         {911.581064518776, -392.759470947783, 4.70433212259104}},
        {face_hy,
         {-625.993390697472, 201.117485184776, -531.671131610641},
         {-625.993390697472, -201.117485184776, -531.671131610641}},
        {face_hy,
         {-568.269417469938, 203.00927120329, -533.88321795581},
         {-568.269417469938, -203.00927120329, -533.88321795581}},
        {face_hy,
         {-325.139807552282, 823.706905731944, 257.196948397389},
         {-325.139807552282, -823.706905731944, 257.196948397389}},
        {face_hy,
         {866.892413271401, 640.788690455485, -176.54000487229},
         {866.892413271401, -640.788690455485, -176.54000487229}},
        {face_hy,
         {-419.309226307149, 498.014883794, -764.036207929972},
         {-419.309226307149, -498.014883794, -764.036207929972}},
        {face_hy,
         {585.236025714257, 783.120110625855, 743.539666556885},
         {585.236025714257, -783.120110625855, 743.539666556885}},
        {face_hy,
         {524.359685879621, 692.384946428762, 578.054745357942},
         {524.359685879621, -692.384946428762, 578.054745357942}},
        {face_hy,
         {897.244726555787, 875.872731479632, 775.61096368399},
         {897.244726555787, -875.872731479632, 775.61096368399}},
        {face_hy,
         {-300.907481496071, 480.625631370619, -301.165637847569},
         {-300.907481496071, -480.625631370619, -301.165637847569}},
        {face_hy,
         {388.166574222673, 933.051708122826, 732.428777397689},
         {388.166574222673, -933.051708122826, 732.428777397689}},
        {face_hy,
         {-942.536579525323, 583.922617587038, -60.2204333836025},
         {-942.536579525323, -583.922617587038, -60.2204333836025}},
        {face_hy,
         {-298.520892306565, 617.646078242717, 222.712547129166},
         {-298.520892306565, -617.646078242717, 222.712547129166}},
        {face_hy,
         {281.692762551692, 838.127088175353, 415.209806372826},
         {281.692762551692, -838.127088175353, 415.209806372826}},
        {face_hy,
         {-892.725428313172, 114.201076939534, 374.782708689109},
         {-892.725428313172, -114.201076939534, 374.782708689109}},
        {face_hy,
         {143.683651937059, 955.915548420687, -846.215168724028},
         {143.683651937059, -955.915548420687, -846.215168724028}},
        {face_hy,
         {106.156621713193, 33.3412090749507, 980.196692262411},
         {106.156621713193, -33.3412090749507, 980.196692262411}},
        {face_hy,
         {-459.242923691718, 247.776523714609, 700.744204376099},
         {-459.242923691718, -247.776523714609, 700.744204376099}},
        {face_hy,
         {-456.949623038379, 118.247707848002, -637.102082984234},
         {-456.949623038379, -118.247707848002, -637.102082984234}},
        {face_hy,
         {-579.059225274804, 273.773120891202, 175.416426878204},
         {-579.059225274804, -273.773120891202, 175.416426878204}},
        {face_hy,
         {872.020646049572, 257.334224771093, -292.257269750882},
         {872.020646049572, -257.334224771093, -292.257269750882}},
        {face_hy,
         {-293.897949590494, 855.89873467486, -811.673221686376},
         {-293.897949590494, -855.89873467486, -811.673221686376}},
        {face_hy,
         {-776.349120404947, 992.420787225825, -42.843122579483},
         {-776.349120404947, -992.420787225825, -42.843122579483}},
        {face_hy,
         {145.126855015214, 27.5903459940014, -692.553024455991},
         {145.126855015214, -27.5903459940014, -692.553024455991}},
        {face_hy,
         {715.444597266735, 292.276527135942, 248.282168459894},
         {715.444597266735, -292.276527135942, 248.282168459894}},
        {face_hy,
         {-638.976865408622, 707.084381651792, 51.8148037407591},
         {-638.976865408622, -707.084381651792, 51.8148037407591}},
        {face_hy,
         {437.204444504621, 773.064226798618, -442.568324140597},
         {437.204444504621, -773.064226798618, -442.568324140597}},
        {face_hy,
         {-287.576723697824, 342.990471827284, 718.650670190588},
         {-287.576723697824, -342.990471827284, 718.650670190588}},
        {face_hy,
         {-115.291180675498, 990.830809434905, -419.866802840993},
         {-115.291180675498, -990.830809434905, -419.866802840993}},
        {face_hy,
         {-185.343414401515, 986.463430963636, 224.859815437623},
         {-185.343414401515, -986.463430963636, 224.859815437623}},
        {face_hy,
         {278.598148086486, 326.152870026196, 228.722381993676},
         {278.598148086486, -326.152870026196, 228.722381993676}},
        {face_hy,
         {-986.870341763274, 78.4623263560739, 982.015869498526},
         {-986.870341763274, -78.4623263560739, 982.015869498526}},
        {face_hy,
         {-650.381607813343, 351.420513801289, 501.628576694939},
         {-650.381607813343, -351.420513801289, 501.628576694939}},
        {face_hy,
         {284.127567891223, 317.243616065089, -529.66901444014},
         {284.127567891223, -317.243616065089, -529.66901444014}},
        {face_hy,
         {-863.41464943173, 709.276002255996, 800.878028954046},
         {-863.41464943173, -709.276002255996, 800.878028954046}},
        {face_hy,
         {1.92581318412704, 610.251779199989, -791.673276822107},
         {1.92581318412704, -610.251779199989, -791.673276822107}},
        {face_hy,
         {745.510178690193, 705.529957578483, -552.216306794638},
         {745.510178690193, -705.529957578483, -552.216306794638}},
        {face_hy,
         {-101.128299886379, 610.448520778059, 285.223092860301},
         {-101.128299886379, -610.448520778059, 285.223092860301}},
        {face_hy,
         {-219.543812436387, 260.837200207668, 353.503513004539},
         {-219.543812436387, -260.837200207668, 353.503513004539}},
        {face_hy,
         {-882.893315124743, 706.121960580093, 137.673458745508},
         {-882.893315124743, -706.121960580093, 137.673458745508}},
        {face_hy,
         {408.343964045841, 271.57583659758, -489.993770018786},
         {408.343964045841, -271.57583659758, -489.993770018786}},
        {face_hy,
         {689.505155428717, 179.388314573675, 809.331340155338},
         {689.505155428717, -179.388314573675, 809.331340155338}},
        {face_hy,
         {-406.707044505717, 212.083379150134, 325.870376084849},
         {-406.707044505717, -212.083379150134, 325.870376084849}},
        {face_hy,
         {708.70692195399, 791.209698512825, -368.514593229001},
         {708.70692195399, -791.209698512825, -368.514593229001}},
        {face_hy,
         {26.3117307140446, 690.187378720799, -865.860301592169},
         {26.3117307140446, -690.187378720799, -865.860301592169}},
        {face_hy,
         {-95.6137498020457, 708.745448336093, 464.857584182731},
         {-95.6137498020457, -708.745448336093, 464.857584182731}},
        {face_hy,
         {-506.435090675513, 402.156789250986, -89.0236058783376},
         {-506.435090675513, -402.156789250986, -89.0236058783376}},
        {face_hy,
         {-828.645837538478, 9.12296654329657, -805.690705748172},
         {-828.645837538478, -9.12296654329657, -805.690705748172}},
        {face_hy,
         {684.667188962493, 946.04096886651, 481.448455418788},
         {684.667188962493, -946.04096886651, 481.448455418788}},
        {face_hy,
         {-262.440115364344, 538.825567161514, 598.113341731745},
         {-262.440115364344, -538.825567161514, 598.113341731745}},
        {face_hy,
         {-856.454105732829, 440.426973581876, -794.940991103271},
         {-856.454105732829, -440.426973581876, -794.940991103271}},
        {face_hy,
         {-960.538093901378, 919.369667338673, -606.57473194454},
         {-960.538093901378, -919.369667338673, -606.57473194454}},
        {face_hy,
         {-925.790619817641, 415.680417665599, -698.952157646932},
         {-925.790619817641, -415.680417665599, -698.952157646932}},
        {face_hy,
         {553.417887798515, 535.483305369734, 668.248698048153},
         {553.417887798515, -535.483305369734, 668.248698048153}},
        {face_hy,
         {968.561205014015, 804.237245198788, 879.050858818334},
         {968.561205014015, -804.237245198788, 879.050858818334}},
        {face_hy,
         {624.475893384827, 268.773330313915, -480.512171923734},
         {624.475893384827, -268.773330313915, -480.512171923734}},
        {face_hy,
         {746.980971115403, 911.211381616316, 26.5453701241422},
         {746.980971115403, -911.211381616316, 26.5453701241422}},
        {face_hy,
         {-984.230638670387, 266.035249392794, -901.446100857861},
         {-984.230638670387, -266.035249392794, -901.446100857861}},
        {face_hy,
         {510.308538298956, 272.268926403569, 765.376263287156},
         {510.308538298956, -272.268926403569, 765.376263287156}},
        {face_hy,
         {485.084793012974, 65.8975573310122, 598.25007840137},
         {485.084793012974, -65.8975573310122, 598.25007840137}},
        {face_hy,
         {-600.639340680485, 128.792746071921, -982.452509615742},
         {-600.639340680485, -128.792746071921, -982.452509615742}},
        {face_hy,
         {255.416934848683, 585.895094978547, 387.666906354863},
         {255.416934848683, -585.895094978547, 387.666906354863}},
        {face_hy,
         {-724.629107611176, 4.53322412483703, -997.786943914479},
         {-724.629107611176, -4.53322412483703, -997.786943914479}},
        {face_hy,
         {-696.255149056623, 683.046253093787, 607.287694419274},
         {-696.255149056623, -683.046253093787, 607.287694419274}},
        {face_hy,
         {270.301560822305, 230.871677003524, -568.370467704564},
         {270.301560822305, -230.871677003524, -568.370467704564}},
        {face_hy,
         {778.785600346282, 598.378013514975, 215.46051164451},
         {778.785600346282, -598.378013514975, 215.46051164451}},
        {face_hy,
         {564.768842306258, 399.980834787384, 986.573234921948},
         {564.768842306258, -399.980834787384, 986.573234921948}},
        {face_hy,
         {-855.842599319878, 566.13969135284, -800.149446260001},
         {-855.842599319878, -566.13969135284, -800.149446260001}},
        {face_hy,
         {253.17034236467, 602.393451068556, 146.940226508168},
         {253.17034236467, -602.393451068556, 146.940226508168}},
        {face_hy,
         {93.2719450413038, 412.580742947157, 226.888806131402},
         {93.2719450413038, -412.580742947157, 226.888806131402}},
        {face_hy,
         {-471.398450792949, 644.219832539138, 534.529581091979},
         {-471.398450792949, -644.219832539138, 534.529581091979}},
        {face_hy,
         {177.908349534864, 400.742708070754, 260.848549741512},
         {177.908349534864, -400.742708070754, 260.848549741512}},
        {face_hy,
         {-325.66309272529, 597.285187264145, -927.018157609554},
         {-325.66309272529, -597.285187264145, -927.018157609554}},
        {face_hy,
         {383.055401703869, 252.467515518198, 925.664174495809},
         {383.055401703869, -252.467515518198, 925.664174495809}},
        {face_hy,
         {-427.106884641442, 285.172720590832, 277.518686921026},
         {-427.106884641442, -285.172720590832, 277.518686921026}},
        {face_hy,
         {737.135188403215, 125.966355092768, 98.793383301379},
         {737.135188403215, -125.966355092768, 98.793383301379}},
        {face_hy,
         {654.171568060202, 543.79938404239, -713.96548960344},
         {654.171568060202, -543.79938404239, -713.96548960344}},
        {face_hy,
         {764.307299171368, 553.108637173103, -754.578098450153},
         {764.307299171368, -553.108637173103, -754.578098450153}},
        {face_hy,
         {463.550760650644, 727.176493492028, 736.519435620767},
         {463.550760650644, -727.176493492028, 736.519435620767}},
        {face_hy,
         {-276.032646033136, 465.892084296532, -465.763770731111},
         {-276.032646033136, -465.892084296532, -465.763770731111}},
        {face_hy,
         {-75.3146392465369, 928.157931898239, 682.945974857608},
         {-75.3146392465369, -928.157931898239, 682.945974857608}},
        {face_hy,
         {115.944831195386, 260.292167792249, 992.487435682268},
         {115.944831195386, -260.292167792249, 992.487435682268}},
        {face_hy,
         {-298.174088557662, 936.522868550669, 961.711773841151},
         {-298.174088557662, -936.522868550669, 961.711773841151}},
        {face_hy,
         {418.780940336478, 198.969494455739, 168.989835256553},
         {418.780940336478, -198.969494455739, 168.989835256553}},
        {face_hy,
         {873.243212358123, 683.35489071528, 351.128799057643},
         {873.243212358123, -683.35489071528, 351.128799057643}},
        {face_hy,
         {-237.886738665195, 593.295294767434, 266.230331890287},
         {-237.886738665195, -593.295294767434, 266.230331890287}},
        {face_hy,
         {-542.811230875982, 706.83364345619, 235.987045430755},
         {-542.811230875982, -706.83364345619, 235.987045430755}},
        {face_hy,
         {-630.985335665582, 400.83799501979, 104.101929945014},
         {-630.985335665582, -400.83799501979, 104.101929945014}},
        {face_hy,
         {-942.550308844471, 293.183740431981, 642.906823423133},
         {-942.550308844471, -293.183740431981, 642.906823423133}},
        {face_hy,
         {-895.925867003155, 438.929299573133, 767.21358646262},
         {-895.925867003155, -438.929299573133, 767.21358646262}},
        {face_hy,
         {815.537936650053, 622.320886319116, 160.211603276537},
         {815.537936650053, -622.320886319116, 160.211603276537}},
        {face_hy,
         {-823.550305008003, 366.105943680593, -989.422843137351},
         {-823.550305008003, -366.105943680593, -989.422843137351}},
        {face_hy,
         {-222.750166567013, 593.699347567781, -203.831671938436},
         {-222.750166567013, -593.699347567781, -203.831671938436}},
        {face_hy,
         {627.756002859593, 670.695915630895, 304.033826750127},
         {627.756002859593, -670.695915630895, 304.033826750127}},
        {face_hy,
         {-885.9625659989, 953.553592916887, -643.134031302816},
         {-885.9625659989, -953.553592916887, -643.134031302816}},
        {face_hy,
         {-32.3757359813158, 361.531494985536, 665.978804786432},
         {-32.3757359813158, -361.531494985536, 665.978804786432}},
        {face_hy,
         {-418.755376746301, 700.742051365307, -273.006772306289},
         {-418.755376746301, -700.742051365307, -273.006772306289}},
        {face_hy,
         {741.475824150541, 4.91159275270547, 323.047920050814},
         {741.475824150541, -4.91159275270547, 323.047920050814}},
        {face_hy,
         {296.213349298208, 744.779666127922, -703.103488243974},
         {296.213349298208, -744.779666127922, -703.103488243974}},
        {face_hy,
         {-228.098795684622, 23.1530689230372, -616.50407096331},
         {-228.098795684622, -23.1530689230372, -616.50407096331}},
        {face_hz,
         {-142.589354985652, 948.907449593867, 749.240783812209},
         {-142.589354985652, 948.907449593867, -749.240783812209}},
        {face_hz,
         {-760.009176273358, -554.181553713331, 680.631313649808},
         {-760.009176273358, -554.181553713331, -680.631313649808}},
        {face_hz,
         {814.747074671149, -875.106353260133, 470.790300543098},
         {814.747074671149, -875.106353260133, -470.790300543098}},
        {face_hz,
         {-916.967011839543, 339.673582280558, 961.598390552482},
         {-916.967011839543, 339.673582280558, -961.598390552482}},
        {face_hz,
         {760.497332539585, -55.9064273183953, 779.705518223494},
         {760.497332539585, -55.9064273183953, -779.705518223494}},
        {face_hz,
         {354.132306348135, -104.865397647635, 933.169825673511},
         {354.132306348135, -104.865397647635, -933.169825673511}},
        {face_hz,
         {439.13208969668, -473.140073441969, 387.750954209216},
         {439.13208969668, -473.140073441969, -387.750954209216}},
        {face_hz,
         {404.904209262001, -255.629483867829, 149.762482157846},
         {404.904209262001, -255.629483867829, -149.762482157846}},
        {face_hz,
         {686.757950514328, -461.236200164356, 634.084734306998},
         {686.757950514328, -461.236200164356, -634.084734306998}},
        {face_hz,
         {-909.584231648876, -936.653885536158, 610.507854705144},
         {-909.584231648876, -936.653885536158, -610.507854705144}},
        {face_hz,
         {126.679255928817, -54.7870417704639, 277.008692424708},
         {126.679255928817, -54.7870417704639, -277.008692424708}},
        {face_hz,
         {472.413627665617, -213.407056130251, 422.182764205493},
         {472.413627665617, -213.407056130251, -422.182764205493}},
        {face_hz,
         {737.068951908428, -292.485638134624, 168.180195950911},
         {737.068951908428, -292.485638134624, -168.180195950911}},
        {face_hz,
         {172.165497544355, -162.378892697337, 198.816824711202},
         {172.165497544355, -162.378892697337, -198.816824711202}},
        {face_hz,
         {-800.094872744446, 141.734813520617, 462.373662361684},
         {-800.094872744446, 141.734813520617, -462.373662361684}},
        {face_hz,
         {455.18908425831, 614.8360024657, 127.821495857755},
         {455.18908425831, 614.8360024657, -127.821495857755}},
        {face_hz,
         {371.419837588136, -485.385586042513, 607.458328525789},
         {371.419837588136, -485.385586042513, -607.458328525789}},
        {face_hz,
         {-22.3822449705121, 831.811350452785, 87.6456325682616},
         {-22.3822449705121, 831.811350452785, -87.6456325682616}},
        {face_hz,
         {-78.4860086106942, 551.733720125018, 93.0876890656209},
         {-78.4860086106942, 551.733720125018, -93.0876890656209}},
        {face_hz,
         {62.4660877276178, 914.401277989567, 201.589070040786},
         {62.4660877276178, 914.401277989567, -201.589070040786}},
        {face_hz,
         {408.147264538527, -698.471667880203, 736.626003420064},
         {408.147264538527, -698.471667880203, -736.626003420064}},
        {face_hz,
         {343.525890901625, 451.626369133957, 441.532582113977},
         {343.525890901625, 451.626369133957, -441.532582113977}},
        {face_hz,
         {662.54149695501, 97.0401368615217, 777.843725515985},
         {662.54149695501, 97.0401368615217, -777.843725515985}},
        {face_hz,
         {-859.171444259512, 49.9281816451862, 272.501383433104},
         {-859.171444259512, 49.9281816451862, -272.501383433104}},
        {face_hz,
         {-844.276776027989, -106.987806944381, 338.321635037784},
         {-844.276776027989, -106.987806944381, -338.321635037784}},
        {face_hz,
         {-600.137532975097, -142.855693389705, 997.981850230843},
         {-600.137532975097, -142.855693389705, -997.981850230843}},
        {face_hz,
         {-620.115582187482, 422.265358531967, 652.202171653},
         {-620.115582187482, 422.265358531967, -652.202171653}},
        {face_hz,
         {-571.709303901497, -143.900342522645, 154.75105122493},
         {-571.709303901497, -143.900342522645, -154.75105122493}},
        {face_hz,
         {-467.40507404322, -29.2373803266469, 207.949907289907},
         {-467.40507404322, -29.2373803266469, -207.949907289907}},
        {face_hz,
         {840.410328419826, -789.060105743441, 517.895707334375},
         {840.410328419826, -789.060105743441, -517.895707334375}},
        {face_hz,
         {-905.838534179775, 329.583445965076, 953.747153313673},
         {-905.838534179775, 329.583445965076, -953.747153313673}},
        {face_hz,
         {-870.888039887825, 210.708232662118, 519.650825778921},
         {-870.888039887825, 210.708232662118, -519.650825778921}},
        {face_hz,
         {-647.697809152584, 160.026894415307, 903.259449955151},
         {-647.697809152584, 160.026894415307, -903.259449955151}},
        {face_hz,
         {32.6494899090749, 544.255744126703, 812.116993264069},
         {32.6494899090749, 544.255744126703, -812.116993264069}},
        {face_hz,
         {807.663586609188, -694.600411131477, 689.609384303542},
         {807.663586609188, -694.600411131477, -689.609384303542}},
        {face_hz,
         {-30.3545017038823, 536.577971405066, 574.114425524989},
         {-30.3545017038823, 536.577971405066, -574.114425524989}},
        {face_hz,
         {-577.443112005014, 254.190987031736, 510.557219721959},
         {-577.443112005014, 254.190987031736, -510.557219721959}},
        {face_hz,
         {68.9010786184731, -379.223291754145, 625.986095974577},
         {68.9010786184731, -379.223291754145, -625.986095974577}},
        {face_hz,
         {877.379044814946, 938.633276172036, 452.667388144136},
         {877.379044814946, 938.633276172036, -452.667388144136}},
        {face_hz,
         {267.455202765665, -668.956354712748, 545.663085813109},
         {267.455202765665, -668.956354712748, -545.663085813109}},
        {face_hz,
         {-195.097369778303, -895.174590127321, 381.53837414486},
         {-195.097369778303, -895.174590127321, -381.53837414486}},
        {face_hz,
         {744.674774321085, -446.008236466472, 216.073220321442},
         {744.674774321085, -446.008236466472, -216.073220321442}},
        {face_hz,
         {-606.300314346873, -291.975538571468, 121.844208807344},
         {-606.300314346873, -291.975538571468, -121.844208807344}},
        {face_hz,
         {267.897450158793, -20.1092056655511, 394.003723460101},
         {267.897450158793, -20.1092056655511, -394.003723460101}},
        {face_hz,
         {-272.605336086465, 402.410906603408, 624.466803618236},
         {-272.605336086465, 402.410906603408, -624.466803618236}},
        {face_hz,
         {20.6714118196333, 446.793660484574, 709.482511541251},
         {20.6714118196333, 446.793660484574, -709.482511541251}},
        {face_hz,
         {498.195831510515, -453.699956215676, 824.493367540551},
         {498.195831510515, -453.699956215676, -824.493367540551}},
        {face_hz,
         {-168.508116460822, -441.759325927119, 375.395632075489},
         {-168.508116460822, -441.759325927119, -375.395632075489}},
        {face_hz,
         {170.877092159612, -15.2204937787196, 223.324653256473},
         {170.877092159612, -15.2204937787196, -223.324653256473}},
        {face_hz,
         {-264.217783287931, -79.2984316058887, 553.998764974831},
         {-264.217783287931, -79.2984316058887, -553.998764974831}},
        {face_hz,
         {165.63667025951, 638.714389418844, 634.061317754254},
         {165.63667025951, 638.714389418844, -634.061317754254}},
        {face_hz,
         {693.391225811594, 298.891785115736, 793.354630794739},
         {693.391225811594, 298.891785115736, -793.354630794739}},
        {face_hz,
         {866.164802436236, -373.614339065919, 42.3720678559662},
         {866.164802436236, -373.614339065919, -42.3720678559662}},
        {face_hz,
         {-795.460345661273, 251.68969094958, 810.119960493716},
         {-795.460345661273, 251.68969094958, -810.119960493716}},
        {face_hz,
         {844.166701072715, -725.314156176755, 686.294130170422},
         {844.166701072715, -725.314156176755, -686.294130170422}},
        {face_hz,
         {-754.359033656237, -569.688426900686, 192.336260356288},
         {-754.359033656237, -569.688426900686, -192.336260356288}},
        {face_hz,
         {853.487521577529, 872.353616260495, 617.538159046226},
         {853.487521577529, 872.353616260495, -617.538159046226}},
        {face_hz,
         {524.577796540693, -667.186214534342, 609.714638163296},
         {524.577796540693, -667.186214534342, -609.714638163296}},
        {face_hz,
         {285.94936074684, 662.126287006228, 373.739593790438},
         {285.94936074684, 662.126287006228, -373.739593790438}},
        {face_hz,
         {451.530034733986, 690.994850791892, 922.221510187649},
         {451.530034733986, 690.994850791892, -922.221510187649}},
        {face_hz,
         {-936.89704219695, 302.866930075926, 172.550312861579},
         {-936.89704219695, 302.866930075926, -172.550312861579}},
        {face_hz,
         {736.603944119448, -535.715856188027, 979.948868534966},
         {736.603944119448, -535.715856188027, -979.948868534966}},
        {face_hz,
         {-718.413533271344, 464.273641676125, 616.235325152142},
         {-718.413533271344, 464.273641676125, -616.235325152142}},
        {face_hz,
         {-40.5826480528149, 494.299115367226, 434.684710946885},
         {-40.5826480528149, 494.299115367226, -434.684710946885}},
        {face_hz,
         {-580.442522485083, -802.728040952062, 224.654527326504},
         {-580.442522485083, -802.728040952062, -224.654527326504}},
        {face_hz,
         {-412.228898120122, 614.445239778892, 254.877094875849},
         {-412.228898120122, 614.445239778892, -254.877094875849}},
        {face_hz,
         {206.709860346523, 435.974662384199, 84.5498680765945},
         {206.709860346523, 435.974662384199, -84.5498680765945}},
        {face_hz,
         {116.772499159327, -883.182002837284, 853.896209656159},
         {116.772499159327, -883.182002837284, -853.896209656159}},
        {face_hz,
         {156.071434144787, 49.0308045770639, 60.9368056086028},
         {156.071434144787, 49.0308045770639, -60.9368056086028}},
        {face_hz,
         {345.810902253838, -593.371641562292, 597.055764392194},
         {345.810902253838, -593.371641562292, -597.055764392194}},
        {face_hz,
         {668.422702786601, -662.615433184497, 195.650802928418},
         {668.422702786601, -662.615433184497, -195.650802928418}},
        {face_hz,
         {-34.3133688375015, -672.835053282155, 606.380279662103},
         {-34.3133688375015, -672.835053282155, -606.380279662103}},
        {face_hz,
         {-790.964071650288, 799.602739486125, 607.594583830452},
         {-790.964071650288, 799.602739486125, -607.594583830452}},
        {face_hz,
         {-246.815017625636, -133.953818046777, 125.670945028522},
         {-246.815017625636, -133.953818046777, -125.670945028522}},
        {face_hz,
         {712.767098891, 757.21522637108, 205.926901892941},
         {712.767098891, 757.21522637108, -205.926901892941}},
        {face_hz,
         {-701.63823218383, 55.8409570603485, 165.250151049232},
         {-701.63823218383, 55.8409570603485, -165.250151049232}},
        {face_hz,
         {-691.596299807776, 208.556184702319, 789.36439211148},
         {-691.596299807776, 208.556184702319, -789.36439211148}},
        {face_hz,
         {-628.029126437042, 780.48989312333, 102.734510102473},
         {-628.029126437042, 780.48989312333, -102.734510102473}},
        {face_hz,
         {-581.880290273543, -860.479599322929, 113.498078422394},
         {-581.880290273543, -860.479599322929, -113.498078422394}},
        {face_hz,
         {-254.751992659712, -311.957901149625, 418.019875800534},
         {-254.751992659712, -311.957901149625, -418.019875800534}},
        {face_hz,
         {-990.445504562064, -592.165173454454, 243.51117307416},
         {-990.445504562064, -592.165173454454, -243.51117307416}},
        {face_hz,
         {-491.773968189849, -372.684598499453, 398.768795501589},
         {-491.773968189849, -372.684598499453, -398.768795501589}},
        {face_hz,
         {-828.268223708646, -347.817259682897, 78.3405271142451},
         {-828.268223708646, -347.817259682897, -78.3405271142451}},
        {face_hz,
         {78.148609313087, 956.488060111675, 200.274500886163},
         {78.148609313087, 956.488060111675, -200.274500886163}},
        {face_hz,
         {70.7328376485525, 770.983165773446, 661.568548901613},
         {70.7328376485525, 770.983165773446, -661.568548901613}},
        {face_hz,
         {-60.3022771587976, -654.099042661386, 377.135647892083},
         {-60.3022771587976, -654.099042661386, -377.135647892083}},
        {face_hz,
         {777.75198327303, 614.190058669481, 296.935645554182},
         {777.75198327303, 614.190058669481, -296.935645554182}},
        {face_hz,
         {-333.756281835664, -793.360267977809, 592.298435626961},
         {-333.756281835664, -793.360267977809, -592.298435626961}},
        {face_hz,
         {-402.871792684228, -970.17334096721, 362.703271341261},
         {-402.871792684228, -970.17334096721, -362.703271341261}},
        {face_hz,
         {-744.515514963589, 610.589586830212, 84.2834989079233},
         {-744.515514963589, 610.589586830212, -84.2834989079233}},
        {face_hz,
         {-451.852800866536, -494.74055518993, 945.026420719275},
         {-451.852800866536, -494.74055518993, -945.026420719275}},
        {face_hz,
         {-288.802269198705, 803.375016339348, 625.418357152602},
         {-288.802269198705, 803.375016339348, -625.418357152602}},
        {face_hz,
         {-352.236791915057, 628.531030612468, 370.590665302225},
         {-352.236791915057, 628.531030612468, -370.590665302225}},
        {face_hz,
         {705.713357980406, -521.078236200635, 924.588282851862},
         {705.713357980406, -521.078236200635, -924.588282851862}},
        {face_hz,
         {472.458162717989, -647.691025474027, 738.43178705582},
         {472.458162717989, -647.691025474027, -738.43178705582}},
        {face_hz,
         {-854.6145021607, 988.334088336939, 888.046457797077},
         {-854.6145021607, 988.334088336939, -888.046457797077}},
        {face_hz,
         {787.009254550609, 693.527487724782, 293.428064185744},
         {787.009254550609, 693.527487724782, -293.428064185744}},
        {face_hz,
         {-862.477984440612, -519.437640055255, 861.266603450098},
         {-862.477984440612, -519.437640055255, -861.266603450098}},
        {face_hz,
         {-868.690905438034, -527.766692074094, 361.073605607096},
         {-868.690905438034, -527.766692074094, -361.073605607096}},
        {face_hz,
         {-280.293398458194, 738.110968523316, 210.800591086721},
         {-280.293398458194, 738.110968523316, -210.800591086721}}};

const std::tuple<index_t, index_t, index_t, Cartesian_Mesh::Point>
    example2_points[n_example2_points] = {
        {0, 0, 0, {0.299378, 1.51744, 0.435515}},
        {0, 0, 0, {0.734359, 3.55144, -0.785705}},
        {0, 0, 0, {2.01539, 2.93677, 2.88501}},
        {0, 0, 0, {0.912936, 3.49049, 0.866316}},
        {0, 0, 0, {1.42459, 3.34635, 1.94259}},
        {0, 0, 0, {2.37027, 0.242776, -0.176895}},
        {0, 0, 0, {0.364258, 1.06642, 3.0381}},
        {0, 0, 0, {1.81347, 1.4174, -0.481571}},
        {0, 0, 0, {0.285775, 1.75558, 3.29675}},
        {0, 0, 0, {1.44854, 1.46518, -0.168813}},
        {0, 0, 1, {0.579817, 2.15485, 3.88679}},
        {0, 0, 1, {0.207307, 0.101534, 6.19151}},
        {0, 0, 1, {0.0117318, 0.299886, 7.6527}},
        {0, 0, 1, {0.360264, 0.0688, 4.59735}},
        {0, 0, 1, {0.393615, 0.228056, 6.64515}},
        {0, 0, 1, {0.663279, 1.91505, 6.73251}},
        {0, 0, 1, {1.48841, 0.220323, 4.92339}},
        {0, 0, 1, {0.265479, 2.14881, 7.63632}},
        {0, 0, 1, {0.49181, 2.04206, 6.11238}},
        {0, 0, 1, {0.361884, 1.0122, 5.22487}}};

const std::tuple<index_t, index_t, index_t, Cartesian_Mesh::Point>
    example3_points_a[n_example3_points_a] = {
        {0, 0, 0, {12.4499, 7.88437, 1.72463}},
        {0, 0, 0, {11.0627, 7.8388, 3.57411}},
        {0, 0, 0, {11.5748, 8.55784, -0.0905416}},
        {0, 0, 0, {12.2735, 7.80448, 1.40329}},
        {0, 0, 0, {10.3159, 8.58263, 1.55553}},
        {0, 0, 0, {11.9099, 6.69346, 2.21014}},
        {0, 0, 0, {10.666, 5.03022, 2.2075}},
        {0, 0, 0, {10.2761, 5.42246, 1.68062}},
        {0, 0, 0, {11.2509, 5.425, 2.45527}},
        {0, 0, 0, {11.5471, 5.90287, 1.90455}},
        {0, 0, 1, {10.4437, 7.50375, 8.23605}},
        {0, 0, 1, {11.4423, 5.48882, 6.47841}},
        {0, 0, 1, {10.628, 5.30595, 4.45729}},
        {0, 0, 1, {11.28, 7.70344, 6.29462}},
        {0, 0, 1, {10.286, 5.16273, 5.21729}},
        {0, 0, 1, {10.0665, 6.05057, 5.36717}},
        {0, 0, 1, {11.3366, 7.55334, 7.91843}},
        {0, 0, 1, {12.2812, 6.6113, 4.19231}},
        {0, 0, 1, {10.4575, 5.06655, 4.94785}},
        {0, 0, 1, {11.2151, 6.04506, 6.94928}},
        {0, 0, 2, {12.2835, 7.01425, 11.35}},
        {0, 0, 2, {10.783, 5.09215, 10.454}},
        {0, 0, 2, {11.2439, 6.06581, 11.4384}},
        {0, 0, 2, {10.0058, 7.32029, 10.7383}},
        {0, 0, 2, {11.2618, 8.56166, 11.1393}},
        {0, 0, 2, {11.8293, 7.81586, 9.05286}},
        {0, 0, 2, {11.7788, 6.61674, 11.4452}},
        {0, 0, 2, {10.8677, 7.84811, 12.9441}},
        {0, 0, 2, {11.4383, 5.22317, 11.255}},
        {0, 0, 2, {11.5693, 7.0514, 10.0191}},
        {0, 0, 3, {10.1561, 6.95995, 14.4225}},
        {0, 0, 3, {11.5574, 8.55155, 13.7393}},
        {0, 0, 3, {11.3701, 6.8473, 14.8947}},
        {0, 0, 3, {10.0133, 8.05175, 13.2126}},
        {0, 0, 3, {11.8719, 6.07934, 17.3197}},
        {0, 0, 3, {11.9379, 5.45446, 17.1313}},
        {0, 0, 3, {10.8162, 5.76048, 15.9064}},
        {0, 0, 3, {10.633, 5.94175, 13.7076}},
        {0, 0, 3, {12.0021, 7.56385, 17.2005}},
        {0, 0, 3, {10.3407, 5.93518, 17.7746}},
        {0, 0, 4, {11.3611, 6.91503, 18.7874}},
        {0, 0, 4, {11.2001, 5.37041, 17.9199}},
        {0, 0, 4, {11.7552, 5.6045, 20.3501}},
        {0, 0, 4, {10.0327, 5.54495, 19.6107}},
        {0, 0, 4, {11.6428, 8.02966, 18.6277}},
        {0, 0, 4, {10.5419, 8.15004, 19.2976}},
        {0, 0, 4, {10.8887, 7.94712, 19.8476}},
        {0, 0, 4, {10.6212, 7.11531, 21.1931}},
        {0, 0, 4, {10.4841, 5.02567, 21.5158}},
        {0, 0, 4, {11.5779, 7.87962, 18.8017}},
        {0, 1, 0, {10.7424, 10.0067, -0.0250019}},
        {0, 1, 0, {10.7407, 10.9369, -0.800228}},
        {0, 1, 0, {12.3801, 9.62258, 1.34608}},
        {0, 1, 0, {11.5139, 10.474, -0.130728}},
        {0, 1, 0, {12.2535, 10.1628, -0.131074}},
        {0, 1, 0, {10.8009, 9.42416, 0.730755}},
        {0, 1, 0, {12.4137, 9.06453, 2.55955}},
        {0, 1, 0, {11.4027, 10.3329, 0.0275045}},
        {0, 1, 0, {11.1144, 8.76315, 0.567393}},
        {0, 1, 0, {10.0867, 11.6125, 2.61873}},
        {0, 1, 1, {10.2439, 8.91675, 4.90114}},
        {0, 1, 1, {11.9632, 11.3502, 7.92661}},
        {0, 1, 1, {10.7161, 10.2076, 5.59559}},
        {0, 1, 1, {10.7482, 10.5714, 4.86691}},
        {0, 1, 1, {10.844, 11.4899, 4.84485}},
        {0, 1, 1, {10.2337, 10.8853, 4.3489}},
        {0, 1, 1, {10.3525, 9.02603, 4.2298}},
        {0, 1, 1, {12.1749, 9.81276, 8.18912}},
        {0, 1, 1, {11.9157, 9.69253, 4.89073}},
        {0, 1, 1, {11.1634, 12.0893, 5.12786}},
        {0, 1, 2, {12.3585, 10.7324, 12.9861}},
        {0, 1, 2, {11.9745, 9.80216, 11.3562}},
        {0, 1, 2, {12.0776, 10.8642, 11.2046}},
        {0, 1, 2, {10.2224, 10.6165, 12.8707}},
        {0, 1, 2, {10.5875, 11.6954, 8.99175}},
        {0, 1, 2, {11.4041, 9.41442, 10.2392}},
        {0, 1, 2, {10.5406, 10.8331, 9.63922}},
        {0, 1, 2, {12.4948, 12.1771, 11.857}},
        {0, 1, 2, {10.5073, 10.6491, 9.14946}},
        {0, 1, 2, {12.3948, 11.4297, 11.0237}},
        {0, 1, 3, {12.3558, 9.1828, 17.2125}},
        {0, 1, 3, {10.994, 10.7765, 13.8498}},
        {0, 1, 3, {10.9092, 9.38125, 15.9993}},
        {0, 1, 3, {11.9821, 8.74853, 16.8558}},
        {0, 1, 3, {11.4048, 10.1273, 17.7131}},
        {0, 1, 3, {10.0082, 10.8719, 13.8087}},
        {0, 1, 3, {11.4852, 12.0672, 16.7088}},
        {0, 1, 3, {11.538, 11.0335, 16.2148}},
        {0, 1, 3, {11.4197, 9.78268, 17.4227}},
        {0, 1, 3, {10.4294, 10.6467, 16.381}},
        {0, 1, 4, {11.9081, 11.8918, 19.8369}},
        {0, 1, 4, {11.6096, 11.3748, 21.9163}},
        {0, 1, 4, {11.5718, 10.6125, 19.7946}},
        {0, 1, 4, {11.0762, 10.619, 18.0673}},
        {0, 1, 4, {10.0719, 11.0595, 17.8649}},
        {0, 1, 4, {10.3803, 10.7543, 18.2874}},
        {0, 1, 4, {11.1633, 9.5825, 21.9651}},
        {0, 1, 4, {12.0509, 10.6954, 19.3658}},
        {0, 1, 4, {11.9724, 10.116, 21.0084}},
        {0, 1, 4, {11.4173, 11.5949, 21.0434}},
        {0, 2, 0, {11.7601, 15.1383, 2.66047}},
        {0, 2, 0, {11.6204, 12.5416, -0.921045}},
        {0, 2, 0, {11.6815, 14.8543, 1.15006}},
        {0, 2, 0, {12.3192, 15.0604, 1.02653}},
        {0, 2, 0, {10.8418, 14.2453, 0.0789302}},
        {0, 2, 0, {10.2842, 12.7732, -0.805988}},
        {0, 2, 0, {12.3331, 13.4982, 1.02246}},
        {0, 2, 0, {11.7407, 13.4521, 0.0999474}},
        {0, 2, 0, {10.5637, 14.4713, 3.20734}},
        {0, 2, 0, {10.4918, 15.5909, 1.26873}},
        {0, 2, 1, {11.243, 13.29, 5.67257}},
        {0, 2, 1, {11.377, 14.2802, 5.16672}},
        {0, 2, 1, {12.1391, 13.5151, 7.61839}},
        {0, 2, 1, {12.2595, 12.9744, 7.2444}},
        {0, 2, 1, {11.9359, 14.7893, 8.19196}},
        {0, 2, 1, {10.2863, 14.8915, 7.68868}},
        {0, 2, 1, {10.4246, 13.4487, 5.21012}},
        {0, 2, 1, {11.8908, 15.7113, 7.83187}},
        {0, 2, 1, {10.0235, 15.5854, 3.94444}},
        {0, 2, 1, {10.6586, 13.4114, 8.15905}},
        {0, 2, 2, {10.4808, 13.9306, 8.88269}},
        {0, 2, 2, {11.2072, 13.0979, 12.1735}},
        {0, 2, 2, {11.6464, 15.2134, 10.9537}},
        {0, 2, 2, {11.2117, 15.3206, 10.2308}},
        {0, 2, 2, {11.6648, 14.1894, 8.52787}},
        {0, 2, 2, {10.3673, 12.5779, 9.28392}},
        {0, 2, 2, {11.3556, 15.537, 9.36099}},
        {0, 2, 2, {10.8973, 15.5371, 9.69015}},
        {0, 2, 2, {11.6161, 15.2028, 11.4387}},
        {0, 2, 2, {10.0452, 12.7984, 8.6384}},
        {0, 2, 3, {12.2326, 13.7613, 15.5806}},
        {0, 2, 3, {10.6551, 14.2717, 15.0631}},
        {0, 2, 3, {11.9733, 12.5868, 13.4197}},
        {0, 2, 3, {10.3546, 13.3282, 15.6225}},
        {0, 2, 3, {12.4061, 14.0079, 17.305}},
        {0, 2, 3, {10.8522, 15.252, 16.6476}},
        {0, 2, 3, {11.1308, 12.4439, 15.2005}},
        {0, 2, 3, {10.2337, 13.4702, 17.4694}},
        {0, 2, 3, {11.7348, 13.022, 16.1255}},
        {0, 2, 3, {10.7081, 13.831, 13.1783}},
        {0, 2, 4, {12.4475, 14.1338, 18.0432}},
        {0, 2, 4, {12.3864, 13.0904, 19.2181}},
        {0, 2, 4, {10.1142, 12.7164, 21.7218}},
        {0, 2, 4, {11.7193, 13.916, 21.2273}},
        {0, 2, 4, {10.7655, 15.4308, 20.1651}},
        {0, 2, 4, {11.6049, 15.0661, 18.6767}},
        {0, 2, 4, {11.9232, 12.6303, 21.1353}},
        {0, 2, 4, {10.3666, 12.6273, 20.4375}},
        {0, 2, 4, {12.2132, 13.4745, 20.0141}},
        {0, 2, 4, {12.3204, 12.353, 18.8366}},
        {1, 0, 0, {13.6892, 7.44778, 0.816648}},
        {1, 0, 0, {13.9378, 7.31209, 3.48165}},
        {1, 0, 0, {13.6385, 6.82611, 2.60935}},
        {1, 0, 0, {14.379, 8.12636, -0.305341}},
        {1, 0, 0, {12.9721, 7.82858, 3.61116}},
        {1, 0, 0, {14.9494, 8.39762, 1.61496}},
        {1, 0, 0, {12.6303, 7.52279, 3.57749}},
        {1, 0, 0, {13.4219, 7.45455, 1.93209}},
        {1, 0, 0, {14.3329, 5.62775, -0.184235}},
        {1, 0, 0, {12.7504, 5.3823, 2.88932}},
        {1, 0, 1, {13.7076, 5.73347, 7.79288}},
        {1, 0, 1, {14.238, 5.95982, 8.36537}},
        {1, 0, 1, {13.4245, 7.28263, 5.74597}},
        {1, 0, 1, {13.1828, 7.23966, 7.68697}},
        {1, 0, 1, {12.8544, 7.25375, 6.16811}},
        {1, 0, 1, {12.8245, 5.23961, 8.31241}},
        {1, 0, 1, {12.7484, 5.47883, 6.46663}},
        {1, 0, 1, {13.3286, 5.02326, 6.74194}},
        {1, 0, 1, {13.2336, 6.48603, 5.80245}},
        {1, 0, 1, {13.357, 6.53422, 4.48303}},
        {1, 0, 2, {12.9847, 7.2414, 11.2227}},
        {1, 0, 2, {14.9886, 5.66718, 12.1901}},
        {1, 0, 2, {13.7244, 5.55151, 11.4888}},
        {1, 0, 2, {14.2984, 6.46193, 12.4663}},
        {1, 0, 2, {13.1865, 5.72103, 11.4342}},
        {1, 0, 2, {14.1304, 6.27227, 10.313}},
        {1, 0, 2, {13.83, 5.17127, 11.4166}},
        {1, 0, 2, {12.9827, 8.59342, 11.3376}},
        {1, 0, 2, {13.1868, 7.34353, 12.7861}},
        {1, 0, 2, {12.8345, 6.00415, 11.7359}},
        {1, 0, 3, {14.7235, 5.40702, 15.3261}},
        {1, 0, 3, {14.1152, 7.55691, 17.1672}},
        {1, 0, 3, {14.446, 5.9597, 17.0841}},
        {1, 0, 3, {12.959, 7.10302, 16.1388}},
        {1, 0, 3, {13.559, 6.40559, 14.4921}},
        {1, 0, 3, {14.5948, 7.65627, 14.8774}},
        {1, 0, 3, {13.8036, 5.88805, 13.8815}},
        {1, 0, 3, {13.5166, 5.61524, 16.3008}},
        {1, 0, 3, {12.8393, 7.97491, 16.3892}},
        {1, 0, 3, {14.2639, 5.85397, 13.8725}},
        {1, 0, 4, {13.6998, 6.53078, 18.168}},
        {1, 0, 4, {13.8298, 7.84516, 22.4135}},
        {1, 0, 4, {14.1824, 6.81632, 22.2294}},
        {1, 0, 4, {14.8261, 5.66913, 21.0376}},
        {1, 0, 4, {13.0332, 5.70268, 19.6455}},
        {1, 0, 4, {14.5582, 5.35377, 21.6254}},
        {1, 0, 4, {13.6296, 8.31398, 21.7918}},
        {1, 0, 4, {14.3667, 7.2357, 19.9004}},
        {1, 0, 4, {12.6221, 8.50223, 22.4776}},
        {1, 0, 4, {14.0198, 6.14536, 19.6787}},
        {1, 1, 0, {13.5724, 11.9282, -0.445788}},
        {1, 1, 0, {13.0556, 10.2894, 3.56002}},
        {1, 1, 0, {14.565, 9.46111, -0.63268}},
        {1, 1, 0, {14.456, 10.9091, 1.06741}},
        {1, 1, 0, {13.7358, 9.56695, 3.34842}},
        {1, 1, 0, {14.2065, 9.44346, 3.05048}},
        {1, 1, 0, {13.7007, 9.19638, 1.79137}},
        {1, 1, 0, {14.1139, 9.99564, 3.40698}},
        {1, 1, 0, {13.6844, 10.4536, 1.47471}},
        {1, 1, 0, {12.9085, 11.8599, -0.945088}},
        {1, 1, 1, {13.1852, 9.13213, 5.06773}},
        {1, 1, 1, {13.5709, 11.7799, 6.36572}},
        {1, 1, 1, {13.9168, 9.25533, 7.6402}},
        {1, 1, 1, {12.7886, 9.74213, 6.35007}},
        {1, 1, 1, {14.1813, 11.4773, 7.45343}},
        {1, 1, 1, {14.038, 11.9192, 5.10848}},
        {1, 1, 1, {13.1571, 10.6594, 4.39068}},
        {1, 1, 1, {13.1668, 11.6073, 4.74725}},
        {1, 1, 1, {14.4956, 9.05263, 3.7083}},
        {1, 1, 1, {14.4592, 11.9984, 4.75905}},
        {1, 1, 2, {12.5457, 10.5092, 9.99336}},
        {1, 1, 2, {13.1044, 11.3914, 11.3652}},
        {1, 1, 2, {13.8643, 10.9672, 12.6059}},
        {1, 1, 2, {12.8372, 9.83079, 9.60202}},
        {1, 1, 2, {14.6069, 10.223, 11.0123}},
        {1, 1, 2, {13.5306, 10.1548, 11.5571}},
        {1, 1, 2, {13.2486, 10.5908, 9.12841}},
        {1, 1, 2, {14.0889, 9.30952, 9.32741}},
        {1, 1, 2, {13.4385, 9.2021, 9.85994}},
        {1, 1, 2, {12.8517, 10.6791, 8.44437}},
        {1, 1, 3, {13.0658, 10.8801, 14.873}},
        {1, 1, 3, {12.5385, 9.78877, 13.193}},
        {1, 1, 3, {13.2942, 9.08515, 14.1088}},
        {1, 1, 3, {14.2847, 9.18333, 14.9801}},
        {1, 1, 3, {13.4414, 9.58538, 14.1068}},
        {1, 1, 3, {14.8783, 10.8967, 17.4656}},
        {1, 1, 3, {12.952, 9.5263, 15.543}},
        {1, 1, 3, {14.4022, 11.9578, 15.9586}},
        {1, 1, 3, {13.0881, 8.83098, 16.5065}},
        {1, 1, 3, {14.7661, 9.90554, 17.5151}},
        {1, 1, 4, {13.9016, 10.1402, 20.8437}},
        {1, 1, 4, {12.8991, 10.8771, 21.0225}},
        {1, 1, 4, {13.7071, 8.9642, 19.0898}},
        {1, 1, 4, {14.7406, 11.2286, 22.1297}},
        {1, 1, 4, {13.0616, 11.3354, 20.3934}},
        {1, 1, 4, {14.8835, 11.1508, 22.4144}},
        {1, 1, 4, {13.2055, 10.6919, 20.0108}},
        {1, 1, 4, {12.6231, 11.3884, 19.2394}},
        {1, 1, 4, {13.633, 11.7459, 18.9119}},
        {1, 1, 4, {14.0088, 10.2267, 20.5229}},
        {1, 2, 0, {13.5769, 13.7531, 1.7403}},
        {1, 2, 0, {13.2017, 13.4424, 2.17763}},
        {1, 2, 0, {13.6153, 13.175, -0.827368}},
        {1, 2, 0, {14.5286, 12.6558, -0.89956}},
        {1, 2, 0, {14.2495, 15.0099, -0.940966}},
        {1, 2, 0, {13.7703, 14.5826, 2.33904}},
        {1, 2, 0, {13.8391, 12.6601, 1.1295}},
        {1, 2, 0, {14.768, 12.7633, 0.0544268}},
        {1, 2, 0, {14.9099, 14.9055, -0.00986495}},
        {1, 2, 0, {14.7741, 12.5527, 0.858866}},
        {1, 2, 1, {14.6007, 14.0721, 7.09596}},
        {1, 2, 1, {13.549, 14.2726, 5.35286}},
        {1, 2, 1, {12.9906, 14.2363, 8.0456}},
        {1, 2, 1, {12.7779, 13.1075, 4.89952}},
        {1, 2, 1, {14.8868, 14.4282, 6.38314}},
        {1, 2, 1, {12.9225, 15.7215, 8.26845}},
        {1, 2, 1, {13.6671, 14.5595, 5.83931}},
        {1, 2, 1, {13.3101, 12.2985, 7.86195}},
        {1, 2, 1, {14.2962, 14.2903, 7.77966}},
        {1, 2, 1, {12.968, 14.026, 5.78151}},
        {1, 2, 2, {12.6797, 12.6176, 10.9151}},
        {1, 2, 2, {14.2985, 13.8143, 10.8393}},
        {1, 2, 2, {14.0757, 13.089, 9.95009}},
        {1, 2, 2, {14.2112, 15.3278, 9.18505}},
        {1, 2, 2, {12.8804, 15.371, 12.7151}},
        {1, 2, 2, {13.4221, 15.6789, 9.33597}},
        {1, 2, 2, {14.0198, 13.7091, 9.5699}},
        {1, 2, 2, {13.6125, 13.3103, 11.589}},
        {1, 2, 2, {14.9974, 12.9018, 8.85346}},
        {1, 2, 2, {14.8155, 14.0444, 10.4763}},
        {1, 2, 3, {13.5427, 14.7643, 15.5513}},
        {1, 2, 3, {14.6329, 12.8778, 14.2175}},
        {1, 2, 3, {13.9745, 12.6909, 13.9961}},
        {1, 2, 3, {14.354, 14.1508, 17.586}},
        {1, 2, 3, {14.4394, 14.1585, 16.8792}},
        {1, 2, 3, {13.1248, 13.8642, 15.285}},
        {1, 2, 3, {14.9487, 15.3082, 17.692}},
        {1, 2, 3, {12.5247, 14.2927, 13.7825}},
        {1, 2, 3, {12.8602, 14.929, 14.7021}},
        {1, 2, 3, {14.8496, 12.7755, 17.5547}},
        {1, 2, 4, {12.7544, 13.6273, 19.1858}},
        {1, 2, 4, {12.7884, 14.8792, 21.6171}},
        {1, 2, 4, {13.7413, 15.3303, 18.5206}},
        {1, 2, 4, {13.5211, 13.5524, 22.0584}},
        {1, 2, 4, {13.3072, 13.1431, 17.9163}},
        {1, 2, 4, {13.2002, 15.6215, 20.7477}},
        {1, 2, 4, {13.9945, 14.8741, 20.388}},
        {1, 2, 4, {13.2556, 15.7734, 21.7067}},
        {1, 2, 4, {14.5738, 12.7427, 22.3427}},
        {1, 2, 4, {13.9393, 12.7971, 21.6605}}};

// End of file
