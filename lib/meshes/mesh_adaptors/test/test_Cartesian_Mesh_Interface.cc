// Cartesian_Mesh_Interface.cc
// Apr 17, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#include "base/soft_equiv.h"
#include "cartesian/Cartesian_Mesh.h"
#include "mesh_adaptors/Cartesian_Mesh_Interface.h"
#include "gtest/gtest.h"
#include <tuple>

using nut::Cartesian_Mesh_Interface;
using nut_mesh::soft_equiv_os;
using Ray = Cartesian_Mesh_Interface::Ray;
using cell_handle_t = Cartesian_Mesh_Interface::cell_handle_t;
using nut_mesh::Cartesian_Face;
using nut_mesh::Cartesian_Mesh;
using nut_mesh::geom_t;
using nut_mesh::index_t;
using mesh_t = Cartesian_Mesh;

TEST(mesh_adaptors_Cartesian_Interface, instantiate) {
  index_t const nx = 5, ny = 3, nz = 2;
  Cartesian_Mesh mesh{nx, ny, nz};
  Cartesian_Mesh_Interface iface(mesh);
  EXPECT_TRUE(true); // failed on 17 Apr 2019
} // TEST(mesh_adaptors_Cartesian_Interface,instantiate)

mesh_t the_usual_cmi() {
  index_t const nx = 5, ny = 3, nz = 2;
  Cartesian_Mesh mesh{nx, ny, nz};
  return mesh;
}

namespace {
mesh_t the_other_cmi() {
  index_t const nx = 3, ny = 5, nz = 7;
  geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
  geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
  Cartesian_Mesh mesh{nx, ny, nz, dx, dy, dz, xmin, ymin, zmin};
  return mesh;
}
} // namespace

TEST(mesh_adaptors_Cartesian_Interface, intersection) {
  // These tests are all the same tests from Cartesian_Mesh::intersection, the
  // point being that the behavior should be the same.
  mesh_t mesh{the_usual_cmi()};
  Cartesian_Mesh_Interface m{mesh};
  cell_handle_t const c{0};

  // Headed toward HIGH_X
  Ray const r1{{0.5, 0.5, 0.5}, {1.0, 0.0, 0.0}};
  auto const [f1, d1] = m.intersection({r1, c});
  EXPECT_TRUE(f1 == Cartesian_Face(3));
  EXPECT_TRUE(soft_equiv_os(d1, 0.5, "distance to face (1x)"));

  const geom_t s2 = std::sqrt(2.0);
  // If headed directly toward X-Y edge, prefers X
  Ray const r2{{0.5, 0.5, 0.5}, {1.0 / s2, 1.0 / s2, 0.0}};
  auto const [f2, d2] = m.intersection({r2, c});
  EXPECT_TRUE(f2 == Cartesian_Face(3));
  EXPECT_TRUE(soft_equiv_os(d2, s2 / 2.0, "distance to face (2)"));

  const geom_t s3 = std::sqrt(3.0);
  // If headed directly toward X-Y-Z corner, prefers X
  Ray const r3{{0.5, 0.5, 0.5}, {1.0 / s3, 1.0 / s3, 1.0 / s3}};
  auto const [f3, d3] = m.intersection({r3, c});
  EXPECT_TRUE(f3 == Cartesian_Face(3));
  EXPECT_TRUE(soft_equiv_os(d3, s3 / 2.0, "distance to face (3)"));

  // If headed directly toward Y-Z edge, prefers Y
  Ray const r4{{0.5, 0.5, 0.5}, {0.0, 1.0 / s2, 1.0 / s2}};
  auto const [f4, d4] = m.intersection({r4, c});
  EXPECT_TRUE(f4 == Cartesian_Face(4));
  EXPECT_TRUE(soft_equiv_os(d4, s2 / 2.0, "distance to face (4)"));

  // Headed toward LOW_Z
  Ray const r5{{0.5, 0.5, 0.5}, {0.0, 0.0, -1.0}};
  auto const [f5, d5] = m.intersection({r5, c});
  EXPECT_TRUE(f5 == Cartesian_Face(2));
  EXPECT_TRUE(soft_equiv_os(d5, 0.5, "distance to face (5)"));

  // Headed toward LOW_Y
  Ray const r6{{0.5, 0.5, 0.5}, {0.0, -1.0, 0.0}};
  auto const [f6, d6] = m.intersection({r6, c});
  EXPECT_TRUE(f6 == Cartesian_Face(1));
  EXPECT_TRUE(soft_equiv_os(d6, 0.5, "distance to face (6)"));

  // Headed toward LOW_X
  Ray const r7{{0.5, 0.5, 0.5}, {-1.0, 0.0, 0.0}};
  auto const [f7, d7] = m.intersection({r7, c});
  EXPECT_TRUE(f7 == Cartesian_Face(0));
  EXPECT_TRUE(soft_equiv_os(d7, 0.5, "distance to face (7)"));

  // Headed toward HIGH_Z
  Ray const r8{{0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}};
  auto const [f8, d8] = m.intersection({r8, c});
  EXPECT_TRUE(f8 == Cartesian_Face(5));
  EXPECT_TRUE(soft_equiv_os(d8, 0.5, "distance to face (8)"));

  // Headed toward HIGH_Y
  Ray const r9{{0.5, 0.5, 0.5}, {0.0, 1.0, 0.0}};
  auto const [f9, d9] = m.intersection({r9, c});
  EXPECT_TRUE(f9 == Cartesian_Face(4));
  EXPECT_TRUE(soft_equiv_os(d9, 0.5, "distance to face (9)"));

  // Headed toward LOW_X, LOW_Y corner
  Ray const r10{{0.5, 0.5, 0.5}, {-1.0 / s2, -1.0 / s2, 0.0}};
  auto const [f10, d10] = m.intersection({r10, c});
  EXPECT_TRUE(f10 == Cartesian_Face(0));
  EXPECT_TRUE(soft_equiv_os(d10, s2 / 2.0, "distance to face (10)"));

  // Headed toward LOW_X, LOW_Y corner
  Ray const r11{{0.5, 0.5, 0.5}, {0.0, -1.0 / s2, -1.0 / s2}};
  auto const [f11, d11] = m.intersection({r11, c});
  EXPECT_TRUE(f11 == Cartesian_Face(1));
  EXPECT_TRUE(soft_equiv_os(d11, s2 / 2.0, "distance to face (11)"));

  return;
}

// Now test intersection through the Mesh_Interface interface.
TEST(mesh_adaptors_Interface, Cartesian_instantiate) {
  index_t const nx = 5, ny = 3, nz = 2;
  Cartesian_Mesh mesh(nx, ny, nz);
  Cartesian_Mesh_Interface m(mesh);
  EXPECT_TRUE(true); // failed 17 Apr 2019
}

// Now test intersection through the Mesh_Interface interface.
TEST(mesh_adaptors_Interface, Cartesian_intersection) {
  // index_t const nx = 5, ny = 3, nz = 2;
  mesh_t mesh{the_usual_cmi()};
  Cartesian_Mesh_Interface m{mesh};

  cell_handle_t const c{0};

  // Headed toward HIGH_X
  Ray const r1{{0.5, 0.5, 0.5}, {1.0, 0.0, 0.0}};
  auto const [f1, d1] = m.intersection({r1, c});
  EXPECT_TRUE(f1 == Cartesian_Face(3));
  EXPECT_TRUE(soft_equiv_os(d1, 0.5, "distance to face (1x)"));

  const geom_t s2 = std::sqrt(2.0);
  // If headed directly toward X-Y edge, prefers X
  Ray const r2{{0.5, 0.5, 0.5}, {1.0 / s2, 1.0 / s2, 0.0}};
  auto const [f2, d2] = m.intersection({r2, c});
  EXPECT_TRUE(f2 == Cartesian_Face(3));
  EXPECT_TRUE(soft_equiv_os(d2, s2 / 2.0, "distance to face (2)"));

  const geom_t s3 = std::sqrt(3.0);
  // If headed directly toward X-Y-Z corner, prefers X
  Ray const r3{{0.5, 0.5, 0.5}, {1.0 / s3, 1.0 / s3, 1.0 / s3}};
  auto const [f3, d3] = m.intersection({r3, c});
  EXPECT_TRUE(f3 == Cartesian_Face(3));
  EXPECT_TRUE(soft_equiv_os(d3, s3 / 2.0, "distance to face (3)"));

  // If headed directly toward Y-Z edge, prefers Y
  Ray const r4{{0.5, 0.5, 0.5}, {0.0, 1.0 / s2, 1.0 / s2}};
  auto const [f4, d4] = m.intersection({r4, c});
  EXPECT_TRUE(f4 == Cartesian_Face(4));
  EXPECT_TRUE(soft_equiv_os(d4, s2 / 2.0, "distance to face (4)"));

  // Headed toward LOW_Z
  Ray const r5{{0.5, 0.5, 0.5}, {0.0, 0.0, -1.0}};
  auto const [f5, d5] = m.intersection({r5, c});
  EXPECT_TRUE(f5 == Cartesian_Face(2));
  EXPECT_TRUE(soft_equiv_os(d5, 0.5, "distance to face (5)"));

  // Headed toward LOW_Y
  Ray const r6{{0.5, 0.5, 0.5}, {0.0, -1.0, 0.0}};
  auto const [f6, d6] = m.intersection({r6, c});
  EXPECT_TRUE(f6 == Cartesian_Face(1));
  EXPECT_TRUE(soft_equiv_os(d6, 0.5, "distance to face (6)"));

  // Headed toward LOW_X
  Ray const r7{{0.5, 0.5, 0.5}, {-1.0, 0.0, 0.0}};
  auto const [f7, d7] = m.intersection({r7, c});
  EXPECT_TRUE(f7 == Cartesian_Face(0));
  EXPECT_TRUE(soft_equiv_os(d7, 0.5, "distance to face (7)"));

  // Headed toward HIGH_Z
  Ray const r8{{0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}};
  auto const [f8, d8] = m.intersection({r8, c});
  EXPECT_TRUE(f8 == Cartesian_Face(5));
  EXPECT_TRUE(soft_equiv_os(d8, 0.5, "distance to face (8)"));

  // Headed toward HIGH_Y
  Ray const r9{{0.5, 0.5, 0.5}, {0.0, 1.0, 0.0}};
  auto const [f9, d9] = m.intersection({r9, c});
  EXPECT_TRUE(f9 == Cartesian_Face(4));
  EXPECT_TRUE(soft_equiv_os(d9, 0.5, "distance to face (9)"));

  // Headed toward LOW_X, LOW_Y corner
  Ray const r10{{0.5, 0.5, 0.5}, {-1.0 / s2, -1.0 / s2, 0.0}};
  auto const [f10, d10] = m.intersection({r10, c});
  EXPECT_TRUE(f10 == Cartesian_Face(0));
  EXPECT_TRUE(soft_equiv_os(d10, s2 / 2.0, "distance to face (10)"));

  // Headed toward LOW_X, LOW_Y corner
  Ray const r11{{0.5, 0.5, 0.5}, {0.0, -1.0 / s2, -1.0 / s2}};
  auto const [f11, d11] = m.intersection({r11, c});
  EXPECT_TRUE(f11 == Cartesian_Face(1));
  EXPECT_TRUE(soft_equiv_os(d11, s2 / 2.0, "distance to face (11)"));

} // TEST(mesh_adaptors_Interface,Cartesian_intersection)

// find_cell test cases from test_Cartesian_Mesh. They should all work the same.
size_t constexpr n_example3_points{300};
extern const std::tuple<index_t, index_t, index_t, Cartesian_Mesh::Point>
    example3_points[n_example3_points];

TEST(mesh_adaptors_Interface, Cartesian_find_cell) {
  using nut_mesh::Vector;
  using geometron_t = Cartesian_Mesh_Interface::geometron_t;
  index_t const nx = 2, ny = 3, nz = 5;
  geom_t const dx = 2.5, dy = 3.6, dz = 4.7;
  geom_t const xmin = 10.0, ymin = 5.0, zmin = -1.0;
  Cartesian_Mesh msh(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
  Cartesian_Mesh_Interface m(msh);
  Vector const dir{1.0, 0.0, 0.0};
  cell_handle_t c{1};
  auto geomo = [&](Vector const &v) { return geometron_t{Ray{v, dir}, c}; };

  for (size_t i = 0; i < n_example3_points; ++i) {
    auto [ix, iy, iz, p] = example3_points[i];
    cell_handle_t expected_cell{
        nut_mesh::cartesian_to_linear(ix, iy, iz, nx, ny)};
    auto const result{m.find_cell(p)};
    EXPECT_TRUE(m.in_cell(geomo(p), result));
    EXPECT_TRUE(m.in_cell(geomo(p), expected_cell));
    EXPECT_EQ(result, expected_cell);
  }
  return;
} // TEST(mesh_adaptors_Interface,Cartesian_find_cell)

TEST(mesh_adaptors_Interface, Cartesian_cell_across) {
  mesh_t mesh{the_other_cmi()};
  Cartesian_Mesh_Interface m{mesh};
  Cartesian_Mesh::Point punused{0.0, 0.0, 0.0};
  {
    cell_handle_t c0{0};
    cell_handle_t xl =
        m.cell_across(c0, Cartesian_Face{Cartesian_Mesh::LOW_X}, punused);
    EXPECT_EQ(xl.as_id(), Cartesian_Mesh::void_cell_idx);
    cell_handle_t xh =
        m.cell_across(c0, Cartesian_Face{Cartesian_Mesh::HIGH_X}, punused);
    EXPECT_EQ(xh, cell_handle_t{1});
    cell_handle_t yl =
        m.cell_across(c0, Cartesian_Face{Cartesian_Mesh::LOW_Y}, punused);
    EXPECT_EQ(yl.as_id(), Cartesian_Mesh::void_cell_idx);
    cell_handle_t yh =
        m.cell_across(c0, Cartesian_Face{Cartesian_Mesh::HIGH_Y}, punused);
    EXPECT_EQ(yh, cell_handle_t{3});
    cell_handle_t zl =
        m.cell_across(c0, Cartesian_Face{Cartesian_Mesh::LOW_Z}, punused);
    EXPECT_EQ(zl.as_id(), Cartesian_Mesh::void_cell_idx);
    cell_handle_t zh =
        m.cell_across(c0, Cartesian_Face{Cartesian_Mesh::HIGH_Z}, punused);
    EXPECT_EQ(zh, cell_handle_t{15});
  }
  {
    cell_handle_t c4{4};
    cell_handle_t xl =
        m.cell_across(c4, Cartesian_Face{Cartesian_Mesh::LOW_X}, punused);
    EXPECT_EQ(xl, cell_handle_t{3});
    cell_handle_t xh =
        m.cell_across(c4, Cartesian_Face{Cartesian_Mesh::HIGH_X}, punused);
    EXPECT_EQ(xh, cell_handle_t{5});
    cell_handle_t yl =
        m.cell_across(c4, Cartesian_Face{Cartesian_Mesh::LOW_Y}, punused);
    EXPECT_EQ(yl, cell_handle_t{1});
    cell_handle_t yh =
        m.cell_across(c4, Cartesian_Face{Cartesian_Mesh::HIGH_Y}, punused);
    EXPECT_EQ(yh, cell_handle_t{7});
    cell_handle_t zl =
        m.cell_across(c4, Cartesian_Face{Cartesian_Mesh::LOW_Z}, punused);
    EXPECT_EQ(zl.as_id(), Cartesian_Mesh::void_cell_idx);
    cell_handle_t zh =
        m.cell_across(c4, Cartesian_Face{Cartesian_Mesh::HIGH_Z}, punused);
    EXPECT_EQ(zh, cell_handle_t{19});
  }
  {
    cell_handle_t c19{19};
    cell_handle_t xl =
        m.cell_across(c19, Cartesian_Face{Cartesian_Mesh::LOW_X}, punused);
    EXPECT_EQ(xl, cell_handle_t{18});
    cell_handle_t xh =
        m.cell_across(c19, Cartesian_Face{Cartesian_Mesh::HIGH_X}, punused);
    EXPECT_EQ(xh, cell_handle_t{20});
    cell_handle_t yl =
        m.cell_across(c19, Cartesian_Face{Cartesian_Mesh::LOW_Y}, punused);
    EXPECT_EQ(yl, cell_handle_t{16});
    cell_handle_t yh =
        m.cell_across(c19, Cartesian_Face{Cartesian_Mesh::HIGH_Y}, punused);
    EXPECT_EQ(yh, cell_handle_t{22});
    cell_handle_t zl =
        m.cell_across(c19, Cartesian_Face{Cartesian_Mesh::LOW_Z}, punused);
    EXPECT_EQ(zl, cell_handle_t{4});
    cell_handle_t zh =
        m.cell_across(c19, Cartesian_Face{Cartesian_Mesh::HIGH_Z}, punused);
    EXPECT_EQ(zh, cell_handle_t{34});
  }
  {
    cell_handle_t c104{104};
    cell_handle_t xl =
        m.cell_across(c104, Cartesian_Face{Cartesian_Mesh::LOW_X}, punused);
    EXPECT_EQ(xl, cell_handle_t{103});
    cell_handle_t xh =
        m.cell_across(c104, Cartesian_Face{Cartesian_Mesh::HIGH_X}, punused);
    EXPECT_EQ(xh, cell_handle_t{Cartesian_Mesh::void_cell_idx});
    cell_handle_t yl =
        m.cell_across(c104, Cartesian_Face{Cartesian_Mesh::LOW_Y}, punused);
    EXPECT_EQ(yl, cell_handle_t{101});
    cell_handle_t yh =
        m.cell_across(c104, Cartesian_Face{Cartesian_Mesh::HIGH_Y}, punused);
    EXPECT_EQ(yh, cell_handle_t{Cartesian_Mesh::void_cell_idx});
    cell_handle_t zl =
        m.cell_across(c104, Cartesian_Face{Cartesian_Mesh::LOW_Z}, punused);
    EXPECT_EQ(zl, cell_handle_t{89});
    cell_handle_t zh =
        m.cell_across(c104, Cartesian_Face{Cartesian_Mesh::HIGH_Z}, punused);
    EXPECT_EQ(zh, cell_handle_t{Cartesian_Mesh::void_cell_idx});
  }
  return;
} // TEST(mesh_adaptors_Interface,Cartesian_cell_across)

/**\brief A class that acts like a random number generator--just returns the
 * next in a sequence of pre-programmed numbers.  */
template <size_t S> struct Buffer_RNG {
  geom_t &random() {
    geom_t &urd = a_[current_];
    current_++;
    if (current_ > a_.size()) {
      current_ = 0;
    }
    return urd;
  }

  Buffer_RNG(geom_t *begin, geom_t *end) : current_(0) {
    std::copy(begin, end, a_.begin());
  }

  std::array<geom_t, S> a_;

  size_t current_;
}; // Buffer_RNG

TEST(mesh_adaptors_Cartesian_Interface, sample_position) {
  constexpr size_t n_urds{3};
  geom_t urds[n_urds]{0.1, 0.2, 0.3};

  mesh_t mesh{the_other_cmi()};
  Cartesian_Mesh_Interface m{mesh};
  /* Run sample_position in a number of cells. Since the buffer RNG is simple,
   * we can easily calculate what we should get. Notionally, could reuse the
   * same RNG, since it should just roll over every three calls. But, this
   * is simpler. */
  {
    Buffer_RNG<3> rng(&urds[0], &urds[n_urds]);
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
    Buffer_RNG<3> rng(&urds[0], &urds[n_urds]);
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
    Buffer_RNG<3> rng(&urds[0], &urds[n_urds]);
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
    Buffer_RNG<3> rng(&urds[0], &urds[n_urds]);
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
    Buffer_RNG<3> rng(&urds[0], &urds[n_urds]);
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
} // TEST(mesh_adaptors_Cartesian_Interface, sample_position)

TEST(mesh_adaptors_Cartesian_Interface, sample_direction_isotropic) {
  constexpr size_t n_urds{2};
  geom_t urds[n_urds]{0.5, 0.5};

  mesh_t mesh{the_other_cmi()};
  Cartesian_Mesh_Interface m{mesh};
  {
    Buffer_RNG<3> rng(&urds[0], &urds[n_urds]);
    nut_mesh::Vector omega = m.sample_direction_isotropic(rng);
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

// ----------------------- test data ----------------------- //

const std::tuple<index_t, index_t, index_t, Cartesian_Mesh::Point>
    example3_points[n_example3_points] = {
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
