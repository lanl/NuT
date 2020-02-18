// T. M. Kelley (c) 2011 LANS LLC

#include "Mesh_3D_Cartesian.hh"
#include "RNG.hh"
#include "expect.hh"
#include "gtest/gtest.h"

#ifdef HAVE_MURMELN

using test_aux::expect;
using Mesh_Iface_T = murmeln::Cartesian_Mesh_Interface;
using Mesh_T = Mesh_Iface_T::mesh_t;
using Cell_T = Mesh_Iface_T::cell_handle_t;
using Index_T = Mesh_Iface_T::index_t;
using Geom_T = Mesh_Iface_T::distance_t;
using Point = Mesh_Iface_T::point_t;

// These are specific to Cartesian_Mesh
using Face = Mesh_Iface_T::face_handle_t;
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
  s << c.as_id();
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
  Mesh_T msh(nx, ny, nz);
  Mesh_Iface_T mesh(msh);
  EXPECT_TRUE(true);
  return;
}

template <typename Index_T>
bool
test_2_core(Index_T const nx, Index_T const ny, Index_T const nz)
{
  // vbd bds;
  // nut::mkReflectBCs<Mesh_T, Cell_T>(bds, nx, ny, nz);
  Mesh_T msh(nx, ny, nz, 1.0, 1.0, 1.0);
  Mesh_Iface_T mesh(msh);
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
  Mesh_T msh(nx, ny, nz, dx, dy, dz);
  Mesh_Iface_T mesh(msh);

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

TEST(nut_mesh_3D_Cartesian, intersection)
{
  Geom_T const dx(2.0);
  Geom_T const dy(3.0);
  Geom_T const dz(5.0);
  Index_T const nx = 1;
  Index_T const ny = 3;
  Index_T const nz = 2;
  Mesh_T msh(nx, ny, nz, dx, dy, dz);
  Mesh_Iface_T mesh(msh);

  Mesh_Iface_T::Vector x, o;
  x = {{0.8, 2.9, 4.9}};
  o = {{1.0, 0.0, 0.0}};

  Mesh_Iface_T::Ray r{x, o};
  Mesh_Iface_T::Intersection d2b = mesh.intersection(r, Cell_T(0));

  Mesh_Iface_T::Geom_T dist = mesh.get_distance(d2b);
  Mesh_Iface_T::Face face = mesh.get_face(d2b);

  bool passed = expect(dist, 1.2, "distance") &&
                expect(face.as_id(), f_h_x.as_id(), "face");
  EXPECT_TRUE(passed);
  return;
}  // test_8

#endif  // HAVE_MURMELN

// End of file
