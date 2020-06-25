// T. M. Kelley (c) 2011 LANS LLC

#define USING_HFBC_SIGMAS

#include "Boundary_Cond.hh"
// #include "Mesh.hh"
#include "Opacity.hh"
#include "Particle.hh"
#include "RNG.hh"
#include "decision.hh"
#include "detail/Vector.h"
#include "events.hh"
#include "gtest/gtest.h"
#include "meshes/mesh_adaptors/Spherical_Mesh_Interface.h"
#include "soft_equiv.hh"
#include "test_aux.hh"
#include "types.hh"
#include <algorithm>
#include <iterator>

using namespace nut::events;

using nut::Boundary_Cond;
using nut::cell_t;
using nut::geom_t;
using nut::Species;
using nut::bdy_types::descriptor;

using fp_t = double;
using vf = std::vector<fp_t>;
using BRNG = nut::Buffer_RNG<fp_t>;
using mesh_t = murmeln::Spherical_Mesh_Interface;
using mesh_impl_t = murmeln_mesh::Spherical_1D_Mesh;
using p_t = nut::Particle<fp_t, BRNG, nut::Vector1>;
using vector_t = nut::Vector1;
using OpB = nut::Opacity<fp_t, vector_t>;
using cell_data_t = nut::Cell_Data<fp_t, vector_t>;
using vec_cell_data_t = std::vector<cell_data_t>;

using face_t = mesh_t::face_handle_t;
using nut::make_vacuum_boundary_1D;

size_t ncells(10);
vf nullv(ncells, fp_t(0));

fp_t const pmg = 1.67262e-24;  // proton mass in gram

/* This test uses the setup of Python test_2 of decide_event to test
 * decide_scatter_event. We use the same random seeds, starting at the
 * third seed, since the test of decide_event burns the first one for
 * sampling the exponential to come up with d_collision, and the second
 * to choose the interaction channel. */
TEST(nut_decide_event, decide_scatter_event_nucleon_abs)
{
  bool passed(true);

  using nut::decide_scatter_event;

  vec_cell_data_t data{{1e14 / pmg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}}};
  OpB op(std::move(data));
  // from Python test_2 of decide_event
  fp_t const rns[] = {// 0.30897681610609407, 1st seed burned for d_coll
                      // 0.92436920169905545,
                      0.21932404923057958};
  BRNG rng(rns, 1);
  fp_t const e = 5.0;
  cell_t const cell = 1;
  Species const s(nut::nu_e);

  Event event = decide_scatter_event(rng, e, cell, op, s);

  Event event_exp = nucleon_abs;

  passed = event == event_exp;
  if(!passed) {
    std::cout << "event was " << event << ", should have been "
              << nut::events::nucleon_abs << std::endl;
  }
  EXPECT_TRUE(passed);
  return;
}  // test_1

// generate a mesh that's reflective at the inner boundary,
// vacuum at the outer boundary, and transmissive at all others
struct gen_bdy_types {
  nut::bdy_types::descriptor operator()()
  {
    if(ctr++ == 0) return nut::bdy_types::REFLECTIVE;
    if(ctr == nbdy) return nut::bdy_types::VACUUM;
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

mesh_impl_t
mk_mesh_1(double dx = 1.0)
{
  using vb = std::vector<geom_t>;
  // generate uniform mesh
  cell_t n_cells(10);
  vb bounds(n_cells + 1);
  // vbd b_types(n_cells + 1);
  // std::generate(b_types.begin(), b_types.end(), gen_bdy_types(n_cells + 1));
  std::generate(bounds.begin(), bounds.end(), gen_bounds(dx));
  // reflect from innermost face
  mesh_impl_t m(std::move(bounds));
  return m;
}

TEST(nut_decide_event, decide_boundary_event)
{
  using nut::events::Event;

  using nut::decide_boundary_event;

  mesh_impl_t m{mk_mesh_1()};
  mesh_t mesh{m};
  Boundary_Cond<mesh_t::face_handle_t> bcs{make_vacuum_boundary_1D(mesh)};
  cell_t const c1(1);
  face_t const f1(0);
  Event const exp_ev1(reflect);
  Event ev1 = decide_boundary_event(mesh, c1, f1, bcs);
  EXPECT_EQ(exp_ev1, ev1);
  // transmit: high face of cell 1
  cell_t const c2(1);
  face_t const f2(1);
  Event const exp_ev2(cell_boundary);
  Event ev2 = decide_boundary_event(mesh, c2, f2, bcs);
  EXPECT_EQ(exp_ev2, ev2);
  // transmit: low face of cell 2
  cell_t const c3(2);
  face_t const f3(2);
  Event const exp_ev3(cell_boundary);
  Event ev3 = decide_boundary_event(mesh, c3, f3, bcs);
  EXPECT_EQ(exp_ev3, ev3);
  // transmit: low face of cell 10
  cell_t const c4(10);
  face_t const f4(10);
  Event const exp_ev4(cell_boundary);
  Event ev4 = decide_boundary_event(mesh, c4, f4, bcs);
  EXPECT_EQ(exp_ev4, ev4);
  // escape: high face of cell 10
  cell_t const c5(10);
  face_t const f5(11);
  Event const exp_ev5(escape);
  Event ev5 = decide_boundary_event(mesh, c5, f5, bcs);
  EXPECT_EQ(exp_ev5, ev5);
  return;
}  // test_2

/** Reworks the Python test_0 of decide_event. */
TEST(nut_decide_event, decide_event_stream_to_cell_boundary)
{
  bool passed(true);

  using nut::decide_scatter_event;
  using nut::events::Event;

  vec_cell_data_t data{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}}};
  OpB op(std::move(data));
  fp_t const rns[] = {0.30897681610609407, 0.92436920169905545,
                      0.21932404923057958};
  BRNG rng(rns, 3);
  fp_t const x = 0.5, omega = 1.0, e = 5.0, t = 1.0, wt = 1.0;
  cell_t const cell = 1;
  Species const s(nut::nu_e);

  mesh_impl_t m{mk_mesh_1()};
  mesh_t mesh{m};
  Boundary_Cond<mesh_t::face_handle_t> bcs{make_vacuum_boundary_1D(mesh)};

  p_t p({x}, {omega}, e, t, wt, cell, rng, s);

  auto [event, distance, fc_out] = decide_event(p, mesh, op, bcs);

  Event const event_exp = cell_boundary;
  // auto [event, fc_out] = nut::events::decode_face<uint32_t>(e_n_d.first);
  passed = event == event_exp;
  if(!passed) {
    std::cout << "event was " << event << ", should have been "
              << nut::events::nucleon_abs << std::endl;
  }
  EXPECT_TRUE(passed);
  return;
}  // test_3

/** Reworks the Python test_1 of decide_event. */
TEST(nut_decide_event, decide_event_stream_through_10_steps)
{
  bool passed(true);

  using nut::decide_scatter_event;
  using nut::events::Event;

  vec_cell_data_t data{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}},
                       {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}},
                       {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}},
                       {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}},
                       {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}},
                       {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}},
                       {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}},
                       {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}},
                       {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}},
                       {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}}};
  OpB op(std::move(data));
  fp_t const rns[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  BRNG rng(rns, 10);
  fp_t const x = 0.5, omega = 1.0, e = 5.0, t = 1.0, wt = 1.0;
  cell_t const cell = 1;
  Species const s(nut::nu_e);

  mesh_impl_t m{mk_mesh_1()};
  mesh_t mesh{m};
  Boundary_Cond<mesh_t::face_handle_t> bcs{make_vacuum_boundary_1D(mesh)};

  p_t p({x}, {omega}, e, t, wt, cell, rng, s);

  for(size_t i = 0; i < 9; ++i) {
    auto [event, distance, fc_out] = decide_event(p, mesh, op, bcs);
    p.x[0] += distance;
    p.cell += 1;
  }

  auto [event, distance, fc_out] = decide_event(p, mesh, op, bcs);
  Event const event_exp = escape;
  passed = event == event_exp;
  if(!passed) {
    std::cout << __FUNCTION__ << ":" << __LINE__ << ":"
              << "event was " << event << ", AKA "
              << nut::events::event_name(event) << ", should have been "
              << event_exp << ", AKA " << nut::events::event_name(event_exp)
              << std::endl;
  }
  geom_t const d_exp = 1.0;
  bool d_passed{distance == d_exp};
  if(!d_passed) {
    printf("%s:%i distance = %e, expected %e\n", __FUNCTION__, __LINE__,
           distance, d_exp);
  }
  EXPECT_TRUE(d_passed);
  passed = d_passed && passed;
  EXPECT_TRUE(passed);
  return;
}  // test_4

/* This test uses the setup of Python test_3 of decide_event to test
 * decide_scatter_event. We use the same random seeds, starting at the
 * third seed, since the test of decide_event burns the first one for
 * sampling the exponential to come up with d_collision, and the second
 * is used to choose the particle interaction. This test
 * should result in nucleon elastic scatter. */
TEST(nut_decide_event, decide_scatter_event_nucleon_elastic_scatter)
{
  bool passed(true);

  using nut::decide_scatter_event;

  vec_cell_data_t data{{1e14 / pmg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}}};
  OpB op(std::move(data));
  // from Python test_2 of decide_event
  fp_t const rns[] = {// 0.21932404923057958, 1st used for d_coll
                      // 0.20867489035315723,
                      0.91525579001682567};
  BRNG rng(rns, 1);
  fp_t const e = 5.0;
  cell_t const cell = 1;
  Species const s(nut::nu_e);

  Event event = decide_scatter_event(rng, e, cell, op, s);

  Event event_exp = nucleon_elastic_scatter;

  passed = event == event_exp;
  if(!passed) {
    std::cout << "event was " << nut::events::event_name(event)
              << ", should have been " << nut::events::nucleon_elastic_scatter
              << std::endl;
  }
  EXPECT_TRUE(passed);
  return;
}  // test_5

/* Stack the deck to use electron scatter */
TEST(nut_decide_event, decide_scatter_event_electron_scatter)
{
  bool passed(true);

  using namespace nut::events;

  using nut::decide_scatter_event;

  vec_cell_data_t data{
      {nut::tiny, 1e14 / pmg, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, {0}}};
  OpB op(std::move(data));

  // from Python test_2 of decide_event
  fp_t const rns[] = {0.9, 0.1};
  BRNG rng(rns, 2);
  fp_t const e = 5.0;
  cell_t const cell = 1;
  Species const s(nut::nu_e);

  Event event = decide_scatter_event(rng, e, cell, op, s);

  Event event_exp = electron_scatter;

  passed = event == event_exp;
  if(!passed) {
    std::cout << "event was " << event_name(event) << ", should have been "
              << event_name(event_exp) << std::endl;
  }
  EXPECT_TRUE(passed);
  return;
}  // test_6

/* This test uses the setup of Python test_2 of decide_event to test
 * decide_event. We use the same random seeds, skipping the
 * second seed, since the test of decide_event burns the first one for
 * sampling the exponential to come up with d_collision. In the Python code,
 * the second seed is used to select which type of particle to scatter from;
 * NuT avoids this in keeping with McPhD. */
TEST(nut_decide_event, decide_event_nucleon_abs)
{
  bool passed(true);

  using namespace nut::events;

  mesh_impl_t m{mk_mesh_1(1.0e6)};
  mesh_t mesh{m};
  Boundary_Cond<mesh_t::face_handle_t> bcs{make_vacuum_boundary_1D(mesh)};

  // generate opacity
  vec_cell_data_t data{{1e14 / pmg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}}};
  OpB op(std::move(data));

  // create particle
  // from Python test_2 of decide_event
  fp_t const rns[] = {0.30897681610609407, 0.92436920169905545,
                      0.21932404923057958};
  BRNG rng(rns, 3);
  fp_t const e = 5.0;
  cell_t const cell = 1;
  Species const s(nut::nu_e);
  geom_t const x = 0.5;
  geom_t const omega = 1.0;
  fp_t const t = 100.0;
  fp_t const wt = 1.0;
  p_t p({x}, {omega}, e, t, wt, cell, rng, s);

  auto [event, d, fc_out] = decide_event(p, mesh, op, bcs);

  Event event_exp = nucleon_abs;

  passed = event == event_exp;
  if(!passed) {
    std::cout << "event was " << event_name(event) << ", should have been "
              << event_name(nut::events::nucleon_abs) << std::endl;
  }
  geom_t const d_exp = 7343.827;
  geom_t const epsilon = 0.001;
  passed = nut::soft_equiv(d, d_exp, epsilon) && passed;
  if(!passed) {
    std::cout << "distance to event was " << std::setprecision(15) << d
              << "expected distance was " << std::setprecision(15) << d_exp
              << std::endl;
  }

  EXPECT_TRUE(passed);
  return;
}  // test_7

/* This test uses the setup of Python test_3 of decide_event to test
 * decide_scatter_event. We use the same random seeds, skipping the
 * second seed, since the test of decide_event burns the first one for
 * sampling the exponential to come up with d_collision. See comment to
 * test_7 for more on the rationale. This test
 * should result in nucleon elastic scatter. */
TEST(nut_decide_event, decide_event_nucleon_elastic_scatter)
{
  bool passed(true);

  using namespace nut::events;

  mesh_impl_t m{mk_mesh_1(1.0e6)};
  mesh_t mesh{m};
  Boundary_Cond<mesh_t::face_handle_t> bcs{make_vacuum_boundary_1D(mesh)};

  vec_cell_data_t data{{1e14 / pmg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}}};
  OpB op(std::move(data));

  // from Python test_2 of decide_event
  fp_t const rns[] = {0.21932404923057958,  // 1st used for d_coll
                      0.20867489035315723, 0.91525579001682567};
  BRNG rng(rns, 3);
  fp_t const e = 5.0;
  cell_t const cell = 1;
  Species const s(nut::nu_e);
  geom_t const x = 0.5;
  geom_t const omega = 1.0;
  fp_t const t = 100.0;
  fp_t const wt = 1.0;
  p_t p({x}, {omega}, e, t, wt, cell, rng, s);

  auto [event, d, fc_out] = nut::decide_event(p, mesh, op, bcs);

  Event const event_exp = nucleon_elastic_scatter;

  passed = event == event_exp;
  if(!passed) {
    std::cout << "event was " << event_name(event) << ", should have been "
              << event_name(nucleon_elastic_scatter) << std::endl;
  }
  geom_t const d_exp = 9486.756;
  geom_t const epsilon = 0.001;
  passed = nut::soft_equiv(d, d_exp, epsilon) && passed;
  if(!passed) {
    std::cout << "distance to event was " << d << "expected distance was "
              << d_exp << std::endl;
  }
  EXPECT_TRUE(passed);
  return;
}  // test_8

/* Stack the deck to use electron scatter */
TEST(nut_decide_event, decide_event_electron_scatter)
{
  bool passed(true);

  using namespace nut::events;

  mesh_impl_t m{mk_mesh_1(1.0e6)};
  mesh_t mesh{m};
  Boundary_Cond<mesh_t::face_handle_t> bcs{make_vacuum_boundary_1D(mesh)};

  vec_cell_data_t data{
      {nut::tiny, 1e14 / pmg, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, {0}}};
  OpB op(std::move(data));
  // from Python test_2 of decide_event
  fp_t const rns[] = {0.9, 0.9, 0.1};
  BRNG rng(rns, 3);
  fp_t const e = 5.0;
  cell_t const cell = 1;
  Species const s(nut::nu_e);
  geom_t const x = 0.5;
  geom_t const omega = 1.0;
  fp_t const t = 1000.0;
  fp_t const wt = 1.0;
  p_t p({x}, {omega}, e, t, wt, cell, rng, s);

  auto [event, d, fc_out] = nut::decide_event(p, mesh, op, bcs);

  Event const event_exp = electron_scatter;
  passed = event == event_exp;
  if(!passed) {
    std::cout << "event was " << event_name(event) << ", should have been "
              << event_name(event_exp) << std::endl;
  }
  fp_t const d_exp = 38310.5;
  geom_t const epsilon = 0.001;
  passed = nut::soft_equiv(d, d_exp, epsilon) && passed;
  if(!passed) {
    std::cout << "distance to event was " << d << "expected distance was "
              << d_exp << std::endl;
  }
  EXPECT_TRUE(passed);
  return;
}  // test_9

// End of file
