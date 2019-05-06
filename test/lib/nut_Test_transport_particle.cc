// T. M. Kelley (c) 2011 LANS LLC

#include "Assert.hh"
#include "Census.hh"
#include "Log.hh"
#include "Mesh.hh"
#include "Particle.hh"
#include "RNG.hh"
#include "gtest/gtest.h"
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

typedef double fp_t;
typedef nut::Buffer_RNG<fp_t> rng_t;
typedef nut::Particle<fp_t, rng_t> p_t;
typedef nut::Tally<fp_t> t_t;
typedef nut::Census<p_t> c_t;
typedef nut::Sphere_1D<cell_t, nut::geom_t, nut::bdy_types::descriptor> mesh_t;


namespace{
// bits for std particle
fp_t const rns[] = {0.30897681610609407, 0.92436920169905545,
                    0.21932404923057958};
rng_t rng(rns, 3);
fp_t const x = 0.5, omega = 1.0, e = 5.0, t = 1.0, wt = 1.0;
cell_t const cell = 1;
Species const s(nut::nu_e);

p_t
make_std_particle()
{
  p_t p(x, omega, e, t, wt, cell, rng, s);
  return p;
}
}  // namespace

// generate a mesh that's reflective at the inner boundary,
// vacuum at the outer boundary, and transmissive at all others
struct gen_bdy_types {
  nut::bdy_types::descriptor operator()()
  {
    if(ctr++ == 0) return nut::bdy_types::R;
    if(ctr == nbdy) return nut::bdy_types::V;
    nut::Insist(ctr <= nbdy, "called gen_bdy_types too often");
    return nut::bdy_types::T;
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

TEST(nut_transport,transport_1_particle)
{
  bool passed(true);

  p_t p(make_std_particle());

  // mesh
  using nut::events::Event;

  using nut::decide_boundary_event;
  typedef mesh_t::vb vb;
  typedef mesh_t::vbd vbd;
  // generate uniform mesh
  cell_t n_cells(10);
  vb bounds(n_cells + 1);
  vbd b_types(n_cells + 1);
  std::generate(b_types.begin(), b_types.end(), gen_bdy_types(n_cells + 1));
  std::generate(bounds.begin(), bounds.end(), gen_bounds(1.0));

  // reflects from innermost face
  mesh_t mesh(bounds, b_types);

  // opacity

  // tally
  nut::Tally<fp_t> tally(n_cells), ref(n_cells);

  // census
  c_t c, c_ref;

  EXPECT_TRUE(passed);
  return;
}

// End of file
