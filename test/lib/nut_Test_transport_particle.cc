// T. M. Kelley (c) 2011 LANS LLC

#include "Assert.hh"
#include "Census.hh"
#include "Log.hh"
#include "Particle.hh"
#include "RNG.hh"
#include "gtest/gtest.h"
#include "meshes/mesh_adaptors/Spherical_Mesh_Interface.h"
#include "test_aux.hh"
#include "transport.hh"
#include <vector>

using nut::cell_t;
using nut::Null_Log;
using nut::Species;
using std::cerr;
using std::endl;
using test_aux::check_one_changed;
using test_aux::check_same;
using test_aux::check_same_verb;
using test_aux::check_two_changed;
using test_aux::comp_verb;

using fp_t = double;
using rng_t = nut::Buffer_RNG<fp_t>;
using p_t = nut::Particle<fp_t, rng_t, nut::Vector1>;
constexpr size_t dim = 1;
using tally_t = nut::Tally<fp_t, dim>;
using c_t = nut::Census<p_t>;
using mesh_t = murmeln::Spherical_Mesh_Interface;
using mesh_impl_t = murmeln_mesh::Spherical_1D_Mesh;

namespace {
// bits for std particle
fp_t const rns[] = {0.30897681610609407, 0.92436920169905545,
                    0.21932404923057958};
rng_t rng(rns, 3);
fp_t const x = 0.5;
fp_t const omega = 1.0;
fp_t const e = 5.0;
fp_t const t = 1.0;
fp_t const wt = 1.0;
cell_t const cell = 1;
Species const s(nut::nu_e);

p_t
make_std_particle()
{
  p_t p({x}, {omega}, e, t, wt, cell, rng, s);
  return p;
}
}  // namespace

// generate a mesh that's reflective at the inner boundary,
// vacuum at the outer boundary, and transmissive at all others
struct gen_bdy_types {
  nut::bdy_types::descriptor operator()()
  {
    if(ctr++ == 0) return nut::bdy_types::REFLECTIVE;
    if(ctr == nbdy) return nut::bdy_types::VACUUM;
    nut::Insist(ctr <= nbdy, "called gen_bdy_types too often");
    return nut::bdy_types::CELL;
  }
  explicit gen_bdy_types(cell_t const nbdy_) : ctr(0), nbdy(nbdy_) {}
  cell_t ctr;
  cell_t const nbdy;
};  // gen_bdy_types

// generate a uniformly-spaced mesh, with spacing dx
struct gen_bounds {
  fp_t operator()() { return dx * ctr++; }
  explicit gen_bounds(fp_t const dx_) : dx(dx_), ctr(0) {}
  fp_t const dx;
  cell_t ctr;
};  // gen_bdy_types

TEST(nutally_transport, transport_1_particle)
{
  bool passed(true);

  p_t p(make_std_particle());

  // mesh
  using nut::events::Event;

  using nut::decide_boundary_event;
  // generate uniform mesh
  cell_t n_cells(10);
  std::vector<fp_t> bounds(n_cells + 1);
  std::generate(bounds.begin(), bounds.end(), gen_bounds(1.0));

  // reflects from innermost face
  mesh_impl_t m(std::move(bounds));
  mesh_t mesh{m};
  // opacity

  // tally
  tally_t tally(n_cells), ref(n_cells);

  // census
  c_t c, c_ref;

  EXPECT_TRUE(passed);
  return;
}

// End of file
