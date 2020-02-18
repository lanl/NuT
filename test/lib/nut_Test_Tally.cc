// T. M. Kelley (c) 2011 LANS LLC

#include "Tally.hh"
#include "gtest/gtest.h"
#include "soft_equiv.hh"
#include "test_aux.hh"
#include "types.hh"
#include <iomanip>

using test_aux::check_one_changed;
using test_aux::check_same;
using test_aux::check_same_v;
using test_aux::check_same_verb;
using test_aux::check_two_changed;
using test_aux::comp_verb;
using test_aux::comp_verb_iter;
using test_aux::tallies_same;

using fp_t = double;
constexpr size_t dim = 1;
using t_t = nut::Tally<fp_t, dim>;

TEST(nut_Tally, inst_init)
{
  bool passed(true);

  size_t n_cells(100);

  nut::Tally<fp_t, dim> tally(n_cells);

  EXPECT_TRUE(passed);
  return;
}  // test_1

TEST(nut_Tally, deposit_inelastic_el_scat)
{
  bool passed(true);

  size_t n_cells(100);

  nut::Tally<fp_t, dim> tally(n_cells), ref(n_cells);

  fp_t const ei(2.0), ef(1.9);
  nut::Vector1 const omega_i{1.0};
  nut::Vector1 const omega_f{-1.0};
  fp_t const wt(0.5);
  nut::Species const s(nut::nu_e);
  nut::cell_t const c(21);

  tally.deposit_inelastic_scat(c, ei, ef, omega_i, omega_f, wt, s);

  ref.momentum[c - 1] = 0.5 * (ei * omega_i - ef * omega_f);
  ref.energy[c - 1] = 0.5 * (ei - ef);

  passed = check_same(&tally.momentum, &ref.momentum) and passed;
  passed = check_same(&tally.energy, &ref.energy) and passed;

  EXPECT_TRUE(passed);
  return;
}  // test_2

TEST(nut_Tally, deposit_energy)
{
  bool passed(true);

  using vf = std::vector<fp_t>;
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  fp_t const wt = 0.2;
  nut::cell_t const c(21);
  fp_t const e = 3.0;
  tally.deposit_energy(c, wt, e);

  passed = check_one_changed<t_t, vf>(tally, ref, &tally.energy) and passed;

  ref.energy[20] = 0.6;

  passed =
      check_same_verb(&tally.energy, &ref.energy, comp_verb<fp_t>()) and passed;

  EXPECT_TRUE(passed);
  return;
}  // test_3

TEST(nut_Tally, deposit_momentum_elastic)
{
  bool passed(true);

  // using vf = std::vector<fp_t>;
  size_t n_cells(100);

  using vv = t_t::vv;

  t_t tally(n_cells), ref(n_cells);

  nut::Vector1 const o = {0.3};
  fp_t const wt = 0.2;
  fp_t const e = 4.0;
  nut::cell_t const c(21);
  tally.deposit_momentum_elastic(c, o, e, wt);

  passed = check_one_changed<t_t, vv>(tally, ref, &tally.momentum) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.momentum[c - 1] = {0.2 * 0.3 * 4.0};
  passed = check_same_verb(&tally.momentum, &ref.momentum,
                           comp_verb_iter<nut::Vector1>()) and
           passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_4

TEST(nut_Tally, count_electron_scatter_nu_e_)
{
  bool passed(true);

  size_t n_cells(100);

  using vc = t_t::vc;

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_e);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  tally.count_lepton_scatter(c, s, n);

  passed =
      check_one_changed<t_t, vc>(tally, ref, &tally.n_nu_e_el_scat) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_nu_e_el_scat[c - 1] = 1;
  passed = check_same(&tally.n_nu_e_el_scat, &ref.n_nu_e_el_scat) and passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_5

TEST(nut_Tally, count_electron_scatter_nu_e_bar_)
{
  bool passed(true);

  size_t n_cells(100);

  using vc = t_t::vc;

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_e_bar);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  tally.count_lepton_scatter(c, s, n);

  passed =
      check_one_changed<t_t, vc>(tally, ref, &tally.n_nu_e_bar_pos_scat) and
      passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_nu_e_bar_pos_scat[c - 1] = 1;
  passed = check_same(&tally.n_nu_e_bar_pos_scat, &ref.n_nu_e_bar_pos_scat) and
           passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_6

TEST(nut_Tally, count_electron_scatter_nu_mu_)
{
  bool passed(true);
  size_t n_cells(100);

  using vc = t_t::vc;

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_mu);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  tally.count_lepton_scatter(c, s, n);

  passed =
      check_one_changed<t_t, vc>(tally, ref, &tally.n_nu_x_el_scat) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_nu_x_el_scat[c - 1] = 1;
  passed = check_same(&tally.n_nu_x_el_scat, &ref.n_nu_x_el_scat) and passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_7

TEST(nut_Tally, count_electron_scatter_nu_mu_bar_)
{
  bool passed(true);
  size_t n_cells(100);

  using vc = t_t::vc;

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_mu_bar);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  tally.count_lepton_scatter(c, s, n);

  passed =
      check_one_changed<t_t, vc>(tally, ref, &tally.n_nu_x_bar_pos_scat) and
      passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_nu_x_bar_pos_scat[c - 1] = 1;
  passed = check_same(&tally.n_nu_x_bar_pos_scat, &ref.n_nu_x_bar_pos_scat) and
           passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_8

TEST(nut_Tally, count_electron_scatter_nu_tau_)
{
  bool passed(true);
  size_t n_cells(100);

  using vc = t_t::vc;

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_tau);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  tally.count_lepton_scatter(c, s, n);

  passed =
      check_one_changed<t_t, vc>(tally, ref, &tally.n_nu_x_el_scat) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_nu_x_el_scat[c - 1] = 1;
  passed = check_same(&tally.n_nu_x_el_scat, &ref.n_nu_x_el_scat) and passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_9

TEST(nut_Tally, count_electron_scatter_nu_tau_bar_)
{
  bool passed(true);
  size_t n_cells(100);

  using vc = t_t::vc;

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_tau_bar);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  tally.count_lepton_scatter(c, s, n);

  passed =
      check_one_changed<t_t, vc>(tally, ref, &tally.n_nu_x_bar_pos_scat) and
      passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_nu_x_bar_pos_scat[c - 1] = 1;
  passed = check_same(&tally.n_nu_x_bar_pos_scat, &ref.n_nu_x_bar_pos_scat) and
           passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_10

TEST(nut_Tally, count_nucleon_abs_nu_e_)
{
  bool passed(true);
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_e);
  nut::cell_t const c(21);
  fp_t const wt(59.2);
  t_t::cntr_t const n(1);
  tally.count_nucleon_abs(c, s, wt, n);

  ref.n_nu_e_nucl_abs[c - 1] = 1;
  ref.ew_nu_e_nucl_abs[c - 1] = 59.2;

  passed = check_same(&tally.n_nu_e_nucl_abs, &ref.n_nu_e_nucl_abs) and passed;
  passed =
      check_same(&tally.ew_nu_e_nucl_abs, &ref.ew_nu_e_nucl_abs) and passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_11

TEST(nut_Tally, count_nucleon_abs_nu_e_bar_)
{
  bool passed(true);
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_e_bar);
  nut::cell_t const c(21);
  fp_t const wt(59.2);
  t_t::cntr_t const n(1);
  tally.count_nucleon_abs(c, s, wt, n);

  ref.n_nu_e_bar_nucl_abs[c - 1] = 1;
  ref.ew_nu_e_bar_nucl_abs[c - 1] = 59.2;

  passed = check_same(&tally.n_nu_e_bar_nucl_abs, &ref.n_nu_e_bar_nucl_abs) and
           passed;
  passed =
      check_same(&tally.ew_nu_e_bar_nucl_abs, &ref.ew_nu_e_bar_nucl_abs) and
      passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_12

TEST(nut_Tally, count_nucleon_elastic_scatter)
{
  bool passed(true);
  size_t n_cells(100);

  using vc = t_t::vc;

  t_t tally(n_cells), ref(n_cells);

  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  tally.count_nucleon_elastic_scatter(c, n);

  passed =
      check_one_changed<t_t, vc>(tally, ref, &tally.n_nucl_el_scat) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_nucl_el_scat[c - 1] = 1;
  passed = check_same(&tally.n_nucl_el_scat, &ref.n_nucl_el_scat) and passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_13

TEST(nut_Tally, count_escape)
{
  bool passed(true);
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  nut::cell_t const c(21);
  fp_t const ew(37.64);
  t_t::cntr_t const n(1);
  tally.count_escape(c, ew, n);

  // two change
  // passed = check_one_changed<t_t,vc>(tally,ref,&tally.n_escape)
  //     and passed;
  // if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_escape[c - 1] = 1;
  ref.ew_escaped[c - 1] = ew;
  passed = check_same(&tally.n_escape, &ref.n_escape) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  passed = check_same(&tally.ew_escaped, &ref.ew_escaped) and passed;

  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_14

TEST(nut_Tally, count_reflect)
{
  bool passed(true);
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  tally.count_reflect(c, n);

  // two change
  // passed = check_one_changed<t_t,vc>(tally,ref,&tally.n_escape)
  //     and passed;
  // if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_reflect[c - 1] = 1;
  passed = check_same(&tally.n_reflect, &ref.n_reflect) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_15

TEST(nut_Tally, count_cell_bdy)
{
  bool passed(true);
  size_t n_cells(100);

  using vc = t_t::vc;

  t_t tally(n_cells), ref(n_cells);

  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  tally.count_cell_bdy(c, n);

  passed = check_one_changed<t_t, vc>(tally, ref, &tally.n_cell_bdy) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_cell_bdy[c - 1] = 1;
  passed = check_same(&tally.n_cell_bdy, &ref.n_cell_bdy) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_16

TEST(nut_Tally, count_cutoff)
{
  bool passed(true);
  size_t n_cells(100);

  using vc = t_t::vc;

  t_t tally(n_cells), ref(n_cells);

  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  tally.count_cutoff(c, n);

  passed = check_one_changed<t_t, vc>(tally, ref, &tally.n_cutoff) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_cutoff[c - 1] = 1;
  passed = check_same(&tally.n_cutoff, &ref.n_cutoff) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_17

TEST(nut_Tally, count_census_nu_e_)
{
  bool passed(true);
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_e);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  fp_t const ew(17.13);
  tally.count_census(c, ew, s, n);

  passed = check_two_changed(tally, ref, &tally.n_census_nu_e,
                             &tally.ew_census_nu_e) and
           passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_census_nu_e[c - 1] = n;
  ref.ew_census_nu_e[c - 1] = ew;
  passed = check_same(&tally.n_census_nu_e, &ref.n_census_nu_e) and passed;
  passed = check_same(&tally.ew_census_nu_e, &ref.ew_census_nu_e) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_18

TEST(nut_Tally, count_census_nu_e_bar_)
{
  bool passed(true);
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_e_bar);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  fp_t const ew(17.13);
  tally.count_census(c, ew, s, n);

  passed = check_two_changed(tally, ref, &tally.n_census_nu_e_bar,
                             &tally.ew_census_nu_e_bar) and
           passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_census_nu_e_bar[c - 1] = n;
  ref.ew_census_nu_e_bar[c - 1] = ew;
  passed =
      check_same(&tally.n_census_nu_e_bar, &ref.n_census_nu_e_bar) and passed;
  passed =
      check_same(&tally.ew_census_nu_e_bar, &ref.ew_census_nu_e_bar) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_19

TEST(nut_Tally, count_census_nu_mu_)
{
  bool passed(true);
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_mu);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  fp_t const ew(17.13);
  tally.count_census(c, ew, s, n);

  passed = check_two_changed(tally, ref, &tally.n_census_nu_x,
                             &tally.ew_census_nu_x) and
           passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_census_nu_x[c - 1] = n;
  ref.ew_census_nu_x[c - 1] = ew;
  passed = check_same(&tally.n_census_nu_x, &ref.n_census_nu_x) and passed;
  passed = check_same(&tally.ew_census_nu_x, &ref.ew_census_nu_x) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_20

TEST(nut_Tally, count_census_nu_mu_bar_)
{
  bool passed(true);
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_mu_bar);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  fp_t const ew(17.13);
  tally.count_census(c, ew, s, n);

  passed = check_two_changed(tally, ref, &tally.n_census_nu_x_bar,
                             &tally.ew_census_nu_x_bar) and
           passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_census_nu_x_bar[c - 1] = n;
  ref.ew_census_nu_x_bar[c - 1] = ew;
  passed =
      check_same(&tally.n_census_nu_x_bar, &ref.n_census_nu_x_bar) and passed;
  passed =
      check_same(&tally.ew_census_nu_x_bar, &ref.ew_census_nu_x_bar) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_21

TEST(nut_Tally, count_census_nu_tau_)
{
  bool passed(true);
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_tau);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  fp_t const ew(17.13);
  tally.count_census(c, ew, s, n);

  passed = check_two_changed(tally, ref, &tally.n_census_nu_x,
                             &tally.ew_census_nu_x) and
           passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_census_nu_x[c - 1] = n;
  ref.ew_census_nu_x[c - 1] = ew;
  passed = check_same(&tally.n_census_nu_x, &ref.n_census_nu_x) and passed;
  passed = check_same(&tally.ew_census_nu_x, &ref.ew_census_nu_x) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_22

TEST(nut_Tally, count_census_nu_tau_bar_)
{
  bool passed(true);
  size_t n_cells(100);

  t_t tally(n_cells), ref(n_cells);

  nut::Species const s(nut::nu_tau_bar);
  nut::cell_t const c(21);
  t_t::cntr_t const n(1);
  fp_t const ew(17.13);
  tally.count_census(c, ew, s, n);

  passed = check_two_changed(tally, ref, &tally.n_census_nu_x_bar,
                             &tally.ew_census_nu_x_bar) and
           passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  ref.n_census_nu_x_bar[c - 1] = n;
  ref.ew_census_nu_x_bar[c - 1] = ew;
  passed =
      check_same(&tally.n_census_nu_x_bar, &ref.n_census_nu_x_bar) and passed;
  passed =
      check_same(&tally.ew_census_nu_x_bar, &ref.ew_census_nu_x_bar) and passed;
  if(!passed) std::cout << "FAILED " << __LINE__ << std::endl;

  EXPECT_TRUE(passed);
  return;
}  // test_23

TEST(nut_Tally, merge)
{
  size_t const n_cells(3);

  using tally_t = nut::Tally<double, dim>;
  tally_t t1(n_cells);
  tally_t t2(n_cells);
  tally_t ref(n_cells);

  for(uint32_t i = 0; i < n_cells; ++i) {
    t1.energy[i] = (double)i;
    t1.momentum[i] = {(double)i};
    t1.n_n[i] = (double)i;
    t1.n_p[i] = (double)i;
    t1.n_e_minus[i] = (double)i;
    t1.n_e_plus[i] = (double)i;
    t1.ew_n[i] = (double)i;
    t1.ew_p[i] = (double)i;
    t1.ew_e_minus[i] = (double)i;
    t1.ew_e_plus[i] = (double)i;
    t1.n_escape[i] = (double)i;
    t1.n_reflect[i] = (double)i;
    t1.n_cell_bdy[i] = (double)i;
    t1.n_cutoff[i] = (double)i;
    t1.n_nucl_el_scat[i] = (double)i;
    t1.n_nu_e_el_scat[i] = (double)i;
    t1.n_nu_e_bar_pos_scat[i] = (double)i;
    t1.n_nu_x_el_scat[i] = (double)i;
    t1.n_nu_x_bar_pos_scat[i] = (double)i;
    t1.ew_escaped[i] = (double)i;
    t1.n_nu_e_nucl_abs[i] = (double)i;
    t1.n_nu_e_bar_nucl_abs[i] = (double)i;
    t1.n_nu_x_nucl_abs[i] = (double)i;
    t1.ew_nu_e_nucl_abs[i] = (double)i;
    t1.ew_nu_e_bar_nucl_abs[i] = (double)i;
    t1.ew_nu_x_nucl_abs[i] = (double)i;
    t1.n_census_nu_e[i] = (double)i;
    t1.n_census_nu_e_bar[i] = (double)i;
    t1.n_census_nu_x[i] = (double)i;
    t1.n_census_nu_x_bar[i] = (double)i;
    t1.ew_census_nu_e[i] = (double)i;
    t1.ew_census_nu_e_bar[i] = (double)i;
    t1.ew_census_nu_x[i] = (double)i;
    t1.ew_census_nu_x_bar[i] = (double)i;

    t2.energy[i] = 2 * (double)i;
    t2.momentum[i] = {2 * (double)i};
    t2.n_n[i] = 2 * (double)i;
    t2.n_p[i] = 2 * (double)i;
    t2.n_e_minus[i] = 2 * (double)i;
    t2.n_e_plus[i] = 2 * (double)i;
    t2.ew_n[i] = 2 * (double)i;
    t2.ew_p[i] = 2 * (double)i;
    t2.ew_e_minus[i] = 2 * (double)i;
    t2.ew_e_plus[i] = 2 * (double)i;
    t2.n_escape[i] = 2 * (double)i;
    t2.n_reflect[i] = 2 * (double)i;
    t2.n_cell_bdy[i] = 2 * (double)i;
    t2.n_cutoff[i] = 2 * (double)i;
    t2.n_nucl_el_scat[i] = 2 * (double)i;
    t2.n_nu_e_el_scat[i] = 2 * (double)i;
    t2.n_nu_e_bar_pos_scat[i] = 2 * (double)i;
    t2.n_nu_x_el_scat[i] = 2 * (double)i;
    t2.n_nu_x_bar_pos_scat[i] = 2 * (double)i;
    t2.ew_escaped[i] = 2 * (double)i;
    t2.n_nu_e_nucl_abs[i] = 2 * (double)i;
    t2.n_nu_e_bar_nucl_abs[i] = 2 * (double)i;
    t2.n_nu_x_nucl_abs[i] = 2 * (double)i;
    t2.ew_nu_e_nucl_abs[i] = 2 * (double)i;
    t2.ew_nu_e_bar_nucl_abs[i] = 2 * (double)i;
    t2.ew_nu_x_nucl_abs[i] = 2 * (double)i;
    t2.n_census_nu_e[i] = 2 * (double)i;
    t2.n_census_nu_e_bar[i] = 2 * (double)i;
    t2.n_census_nu_x[i] = 2 * (double)i;
    t2.n_census_nu_x_bar[i] = 2 * (double)i;
    t2.ew_census_nu_e[i] = 2 * (double)i;
    t2.ew_census_nu_e_bar[i] = 2 * (double)i;
    t2.ew_census_nu_x[i] = 2 * (double)i;
    t2.ew_census_nu_x_bar[i] = 2 * (double)i;

    ref.energy[i] = 3 * (double)i;
    ref.momentum[i] = {3 * (double)i};
    ref.n_n[i] = 3 * (double)i;
    ref.n_p[i] = 3 * (double)i;
    ref.n_e_minus[i] = 3 * (double)i;
    ref.n_e_plus[i] = 3 * (double)i;
    ref.ew_n[i] = 3 * (double)i;
    ref.ew_p[i] = 3 * (double)i;
    ref.ew_e_minus[i] = 3 * (double)i;
    ref.ew_e_plus[i] = 3 * (double)i;
    ref.n_escape[i] = 3 * (double)i;
    ref.n_reflect[i] = 3 * (double)i;
    ref.n_cell_bdy[i] = 3 * (double)i;
    ref.n_cutoff[i] = 3 * (double)i;
    ref.n_nucl_el_scat[i] = 3 * (double)i;
    ref.n_nu_e_el_scat[i] = 3 * (double)i;
    ref.n_nu_e_bar_pos_scat[i] = 3 * (double)i;
    ref.n_nu_x_el_scat[i] = 3 * (double)i;
    ref.n_nu_x_bar_pos_scat[i] = 3 * (double)i;
    ref.ew_escaped[i] = 3 * (double)i;
    ref.n_nu_e_nucl_abs[i] = 3 * (double)i;
    ref.n_nu_e_bar_nucl_abs[i] = 3 * (double)i;
    ref.n_nu_x_nucl_abs[i] = 3 * (double)i;
    ref.ew_nu_e_nucl_abs[i] = 3 * (double)i;
    ref.ew_nu_e_bar_nucl_abs[i] = 3 * (double)i;
    ref.ew_nu_x_nucl_abs[i] = 3 * (double)i;
    ref.n_census_nu_e[i] = 3 * (double)i;
    ref.n_census_nu_e_bar[i] = 3 * (double)i;
    ref.n_census_nu_x[i] = 3 * (double)i;
    ref.n_census_nu_x_bar[i] = 3 * (double)i;
    ref.ew_census_nu_e[i] = 3 * (double)i;
    ref.ew_census_nu_e_bar[i] = 3 * (double)i;
    ref.ew_census_nu_x[i] = 3 * (double)i;
    ref.ew_census_nu_x_bar[i] = 3 * (double)i;
  }
  t1.path_length = 1.0;
  t2.path_length = 2.0;
  ref.path_length = 3.0;

  t2.merge(t1);

  bool passed = tallies_same(t2, ref);

  EXPECT_TRUE(passed);
  return;
}  // test_24

// End of file
