// T. M. Kelley (c) 2011 LANS LLC
// tests of events application

#include "Census.hh"
#include "Mesh.hh"
#include "Particle.hh"
#include "RNG.hh"
#include "Tally.hh"
#include "Vec3D.hh"
#include "Velocity.hh"
#include "apply_event.hh"
#include "expect.hh"
#include "gtest/gtest.h"
#include "test_aux.hh"
#include "types.hh"

using nut::cell_t;
using nut::geom_t;
using nut::Species;
using std::cerr;
using std::endl;
using test_aux::check_one_changed;
using test_aux::check_same;
using test_aux::check_same_v;
using test_aux::check_same_verb;
using test_aux::check_two_changed;
using test_aux::comp_verb;
using test_aux::comp_verb_iter;

using fp_t = double;
using rng_t = nut::Buffer_RNG<fp_t>;

using vector_t = nut::Vec_T<fp_t, 1>;
using p_t = nut::Particle<fp_t, rng_t, vector_t>;
constexpr size_t dim = 1;
using tally_t = nut::Tally<fp_t, dim>;
using c_t = nut::Census<p_t>;
using v_t = nut::Velocity<fp_t, vector_t>;

using mesh_t = nut::Sphere_1D<cell_t, geom_t, nut::bdy_types::descriptor>;

// bits for std particle
fp_t const rns[] = {0.30897681610609407, 0.92436920169905545,
                    0.21932404923057958};
rng_t rng(rns, 3);
vector_t const x = {0.5}, omega = {1.0};
fp_t e = 5.0, t = 1.0, wt = 1.0;
cell_t const cell = 1;
Species const s(nut::nu_e);

p_t
make_std_particle_c1()
{
  p_t p({x}, {omega}, e, t, wt, cell, rng, s);
  return p;
}

p_t
make_std_particle_c2()
{
  vector_t x{1.5};
  cell_t cell{2};
  p_t p({x}, {omega}, e, t, wt, cell, rng, s);
  return p;
}

using test_aux::expect;

TEST(nut_apply_event, apply_nucleon_abs)
{
  bool passed(false);

  p_t p(make_std_particle_c1());

  size_t const n_cells(100);
  tally_t tally(n_cells), ref(n_cells);

  cell_t const idx = cell - 1;
  ref.momentum[idx] = omega * wt * e;
  ref.energy[idx] = wt * e;
  ref.n_nu_e_nucl_abs[idx] = 1;
  ref.ew_nu_e_nucl_abs[idx] = wt;

  nut::apply_nucleon_abs(p, tally);
  // check tally energy deposition, momentum deposition, and counts.
  bool const t_passed =
      check_same_verb(&tally.momentum, &ref.momentum,
                      comp_verb_iter<vector_t>()) &&
      check_same_verb(&tally.energy, &ref.energy, comp_verb<fp_t>("energy")) &&
      check_same_v(&tally.n_nu_e_nucl_abs, &ref.n_nu_e_nucl_abs,
                   "n nu_e nucl abs") &&
      check_same_v(&tally.ew_nu_e_nucl_abs, &ref.ew_nu_e_nucl_abs,
                   "ew nu_e nucl abs");

  // check particle status
  bool const p_passed = expect(p.alive, false, "Particle.alive");
  passed = t_passed && p_passed;
  EXPECT_TRUE(passed);
  return;
}  // test_1

TEST(nut_apply_event, apply_nucleon_elastic_scatter)
{
  bool passed(false);

  p_t p(make_std_particle_c1());
  size_t const n_cells(100);
  tally_t tally(n_cells), ref(n_cells);
  std::vector<vector_t> vs(n_cells);
  v_t vel(vs);

  cell_t const idx = cell - 1;
  vector_t const new_omega = {2 * rns[0] - 1};
  // since we have zero v, there should be zero energy transfer,
  // and it should be as if we did an elastic scatter (cause we did).
  ref.momentum[idx] = wt * e * (omega - new_omega);
  ref.n_nucl_el_scat[idx] = 1;

  nut::apply_nucleon_elastic_scatter<p_t, tally_t, v_t, mesh_t>(p, tally, vel);

  // check tally energy deposition, momentum deposition, and counts.
  bool const t_passed = check_same(&tally.momentum, &ref.momentum) &&
                        check_same(&tally.n_nucl_el_scat, &ref.n_nucl_el_scat);
  if(!t_passed) {
    cerr << "did not tally correctly, line " << __LINE__ << endl;
    std::cout << "actual momentum: ";
    std::copy(tally.momentum.begin(), tally.momentum.end(),
              std::ostream_iterator<vector_t>(std::cout, ","));
    std::cout << std::endl << "expected momentum: ";
    std::copy(ref.momentum.begin(), ref.momentum.end(),
              std::ostream_iterator<vector_t>(std::cout, ","));
    std::cout << std::endl;
  }
  // check particle status
  bool const p_passed = p.omega == new_omega;
  if(!p_passed) { cerr << "incorrect omega: " << p.omega << __LINE__ << endl; }
  passed = t_passed && p_passed;
  EXPECT_TRUE(passed);
  return;
}  // test_2

TEST(nut_apply_event, apply_lepton_scatter_nu_e__electron)
{
  bool passed(false);

  p_t p(make_std_particle_c1());

  size_t const n_cells(100);
  tally_t tally(n_cells), ref(n_cells);
  std::vector<vector_t> vs(n_cells);
  v_t vel(vs);

  cell_t const idx = cell - 1;
  fp_t const new_omega = 2 * rns[0] - 1;

  fp_t const e_e = 3.0;
  fp_t const de = (p.e - e_e) / 4.0;
  fp_t const e_f = p.e - de;
  ref.momentum[idx] = (p.e * p.omega - e_f * new_omega) * wt;
  ref.energy[idx] = wt * de;
  ref.n_nu_e_el_scat[idx] = 1;

  nut::apply_lepton_scatter<p_t, tally_t, v_t, mesh_t>(p, tally, e_e, vel);

  // check tally energy deposition, momentum deposition, and counts.
  bool const t_passed = check_same_verb(&tally.momentum, &ref.momentum,
                                        comp_verb_iter<vector_t>()) &&
                        check_same_verb(&tally.n_nu_e_el_scat,
                                        &ref.n_nu_e_el_scat, comp_verb<fp_t>());
  if(!t_passed) {
    cerr << "did not tally correctly, line " << __LINE__ << endl;
  }
  // check particle status
  bool const p_passed = p.omega == new_omega;
  if(!p_passed) { cerr << "incorrect omega: " << p.omega << __LINE__ << endl; }
  passed = t_passed && p_passed;
  EXPECT_TRUE(passed);
  return;
}  // test_3

TEST(nut_apply_event, apply_lepton_scatter_nu_e_bar__positron)
{
  bool passed(false);

  p_t p(make_std_particle_c1());
  p.species = nut::nu_e_bar;

  size_t const n_cells(100);
  tally_t tally(n_cells), ref(n_cells);
  std::vector<vector_t> vs(n_cells);
  v_t vel(vs);

  cell_t const idx = cell - 1;
  fp_t const new_omega = 2 * rns[0] - 1;

  fp_t const e_e = 3.0;
  fp_t const de = (p.e - e_e) / 4.0;
  fp_t const e_f = p.e - de;
  ref.momentum[idx] = (p.e * p.omega - e_f * new_omega) * wt;
  ref.energy[idx] = wt * de;
  ref.n_nu_e_bar_pos_scat[idx] = 1;

  nut::apply_lepton_scatter<p_t, tally_t, v_t, mesh_t>(p, tally, e_e, vel);

  // check tally energy deposition, momentum deposition, and counts.
  bool const t_passed = check_same(&tally.momentum, &ref.momentum) &&
                        check_same(&tally.n_nu_e_el_scat, &ref.n_nu_e_el_scat);
  if(!t_passed) {
    cerr << "did not tally correctly, line " << __LINE__ << endl;
  }
  // check particle status
  bool const p_passed = p.omega == new_omega;
  if(!p_passed) { cerr << "incorrect omega: " << p.omega << __LINE__ << endl; }
  passed = t_passed && p_passed;
  EXPECT_TRUE(passed);
  return;
}  // test_4

TEST(nut_apply_event, apply_cell_boundary)
{
  size_t constexpr n_cells(3);
  size_t constexpr n_bdys(n_cells + 1);
  geom_t const bdys_in[n_bdys] = {0.0, 1.0, 2.0, 35.0};
  nut::bdy_types::descriptor const bdy_ts[n_bdys] = {
      nut::bdy_types::REFLECTIVE, nut::bdy_types::CELL, nut::bdy_types::CELL,
      nut::bdy_types::VACUUM};
  mesh_t::vb bdys(&bdys_in[0], &bdys_in[n_bdys]);
  mesh_t::vbd bdy_types(&bdy_ts[0], &bdy_ts[n_bdys]);

  mesh_t mesh(bdys, bdy_types);

  {
    bool passed(false);

    p_t p(make_std_particle_c2());

    size_t const n_cells(100);
    tally_t tally(n_cells), ref(n_cells);

    cell_t const idx = 1;  // cell 2, index is 1
    ref.n_cell_bdy[idx] = 1;

    mesh_t::Face f_low{nut::Sphere_1D_Faces::LOW};
    nut::apply_cell_boundary<p_t, tally_t, mesh_t, mesh_t::Face>(mesh, f_low, p,
                                                                 tally);

    // check tally energy deposition, momentum deposition, and counts.
    bool const t_passed = check_same(&tally.n_cell_bdy, &ref.n_cell_bdy) &&
                          check_one_changed(tally, ref, &tally.n_cell_bdy);
    if(!t_passed) {
      cerr << "did not tally correctly, line " << __LINE__ << endl;
    }
    // check particle status
    bool const p_passed = p.cell == 1;
    if(!p_passed) {
      cerr << __LINE__ << ": incorrect cell: " << p.cell << endl;
    }
    passed = t_passed && p_passed;
    EXPECT_TRUE(passed);
  }
  {
    bool passed(false);

    p_t p(make_std_particle_c1());

    size_t const n_cells(100);
    tally_t tally(n_cells), ref(n_cells);

    cell_t const idx = cell - 1;
    ref.n_cell_bdy[idx] = 1;

    mesh_t::Face f_high{nut::Sphere_1D_Faces::HIGH};
    nut::apply_cell_boundary<p_t, tally_t>(mesh, f_high, p, tally);

    // check tally energy deposition, momentum deposition, and counts.
    bool const t_passed = check_same(&tally.n_cell_bdy, &ref.n_cell_bdy) &&
                          check_one_changed(tally, ref, &tally.n_cell_bdy);
    if(!t_passed) {
      cerr << "did not tally correctly, line " << __LINE__ << endl;
    }
    // check particle status
    bool const p_passed = p.cell == cell + 1;
    if(!p_passed) { cerr << "incorrect cell: " << p.cell << __LINE__ << endl; }
    passed = t_passed && p_passed;
    EXPECT_TRUE(passed);
  }
  return;
}  // test_5

TEST(nut_apply_event, apply_escape)
{
  bool passed(false);

  p_t p(make_std_particle_c1());

  size_t const n_cells(100);
  tally_t tally(n_cells), ref(n_cells);

  cell_t const idx = cell - 1;
  ref.n_escape[idx] = 1;
  ref.ew_escaped[idx] = wt;

  nut::apply_escape<p_t, tally_t>(p, tally);

  // check tally energy deposition, momentum deposition, and counts.
  bool const t_passed =
      check_same(&tally.n_escape, &ref.n_escape) &&
      check_same(&tally.ew_escaped, &ref.ew_escaped) &&
      check_two_changed(tally, ref, &tally.n_escape, &tally.ew_escaped);
  if(!t_passed) {
    cerr << "did not tally correctly, line " << __LINE__ << endl;
  }
  // check particle status
  bool const p_passed = p.alive == false;
  if(!p_passed) {
    cerr << "particle not dead: alive = " << p.alive << __LINE__ << endl;
  }
  passed = t_passed && p_passed;
  EXPECT_TRUE(passed);
  return;
}  // test_7

TEST(nut_apply_event, apply_reflect)
{
  bool passed(false);

  p_t p(make_std_particle_c1());

  size_t constexpr n_cells(3);
  size_t constexpr n_bdys(n_cells + 1);
  geom_t const bdys_in[n_bdys] = {0.0, 1.0, 2.0, 35.0};
  nut::bdy_types::descriptor const bdy_ts[n_bdys] = {
      nut::bdy_types::REFLECTIVE, nut::bdy_types::CELL, nut::bdy_types::CELL,
      nut::bdy_types::VACUUM};
  mesh_t::vb bdys(&bdys_in[0], &bdys_in[n_bdys]);
  mesh_t::vbd bdy_types(&bdy_ts[0], &bdy_ts[n_bdys]);

  mesh_t mesh(bdys, bdy_types);

  // size_t const n_cells(100);
  tally_t tally(n_cells), ref(n_cells);

  cell_t const idx = cell - 1;
  ref.n_reflect[idx] = 1;

  // vector_t face_normal{-1.0};
  typename mesh_t::Face face{1};
  nut::apply_reflect<p_t, tally_t>(p, tally, mesh, face);

  // check tally energy deposition, momentum deposition, and counts.
  bool const t_passed = check_same(&tally.n_reflect, &ref.n_reflect) &&
                        check_one_changed(tally, ref, &tally.n_reflect);
  if(!t_passed) {
    cerr << "did not tally correctly, line " << __LINE__ << endl;
  }
  // check particle status
  bool const p_passed = p.omega == -omega;
  if(!p_passed) {
    cerr << "particle not reflected: new omega = " << p.omega << __LINE__
         << endl;
  }
  passed = t_passed && p_passed;
  EXPECT_TRUE(passed);
  return;
}  // test_8

TEST(nut_apply_event, apply_census)
{
  bool passed(false);

  p_t p(make_std_particle_c1());
  p_t p_ref(make_std_particle_c1());

  size_t const n_cells(100);
  tally_t tally(n_cells), ref(n_cells);

  c_t c, c_ref;

  p_ref.alive = false;
  c_ref.append(p_ref);

  cell_t const idx = cell - 1;
  ref.n_census_nu_e[idx] = 1;
  ref.ew_census_nu_e[idx] = wt;

  nut::apply_step_end<p_t, tally_t, c_t>(p, tally, c);

  // check tally energy deposition, momentum deposition, and counts.
  bool const t_passed =
      check_same_verb(&tally.n_census_nu_e, &ref.n_census_nu_e,
                      comp_verb<fp_t>()) &&
      check_two_changed(tally, ref, &tally.n_census_nu_e,
                        &tally.ew_census_nu_e);
  if(!t_passed) {
    cerr << "did not tally correctly, line " << __LINE__ << endl;
  }

  bool const c_passed = check_same(&c.nu_es, &c_ref.nu_es);
  if(!c_passed) { cerr << "census incorrect, line " << __LINE__ << endl; }

  // check particle status
  bool const p_passed = p.alive == false;
  if(!p_passed) {
    cerr << "particle still alive: " << p.alive << __LINE__ << endl;
  }
  passed = t_passed && c_passed && p_passed;
  EXPECT_TRUE(passed);
  return;
}  // test_9

// End of file
