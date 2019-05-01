// T. M. Kelley (c) 2011 LANS LLC

#include "RNG.hh"
// #include "nut_Test_RNG.hh"
// #include "test_aux.hh"
#include "gtest/gtest.h"
#include <algorithm>
#include <functional>
#include <vector>

// make rng callable by STL algos
template <typename fp_t, typename rng_t>
class rng_adaptor {
public:
  explicit rng_adaptor(rng_t & rng) : m_rng(rng) {}
  fp_t operator()() { return m_rng.random(); }
  rng_t & m_rng;
  void dump_state(std::ostream & o) { return m_rng.dump_state(o); }
};

TEST(nut_RNG,init_inst_Buffer_RNG_float)
{
  bool passed(true);

  typedef float fp_t;
  size_t const vals_sz = 4;
  fp_t const in_vals[vals_sz] = {0., 0.11, 0.22, 0.33};

  typedef nut::Buffer_RNG<fp_t> rng_t;
  rng_t rng(&in_vals[0], vals_sz);
  EXPECT_TRUE(passed);
  return;
}  // test_1

TEST(nut_RNG,use_Buffer_RNG_float)
{
  bool passed(true);

  typedef float fp_t;
  size_t const vals_sz = 4;
  fp_t const in_vals[vals_sz] = {0., 0.11, 0.22, 0.33};

  typedef nut::Buffer_RNG<fp_t> rng_t;
  rng_t rng_base(&in_vals[0], vals_sz);
  rng_adaptor<fp_t, rng_t> rng(rng_base);

  std::vector<fp_t> out_vals(vals_sz);

  std::generate(out_vals.begin(), out_vals.end(), rng);
  passed =
      passed and std::equal(&in_vals[0], &in_vals[vals_sz], out_vals.begin());
  EXPECT_TRUE(passed);
  return;
}  // test_2

TEST(nut_RNG,use_Buffer_RNG_double)
{
  bool passed(true);

  typedef double fp_t;
  size_t const vals_sz = 4;
  fp_t const in_vals[vals_sz] = {0., 0.11, 0.22, 0.33};

  typedef nut::Buffer_RNG<fp_t> rng_t;
  rng_t rng_base(&in_vals[0], vals_sz);
  rng_adaptor<fp_t, rng_t> rng(rng_base);

  std::vector<fp_t> out_vals(vals_sz);

  std::generate(out_vals.begin(), out_vals.end(), rng);
  passed =
      passed and std::equal(&in_vals[0], &in_vals[vals_sz], out_vals.begin());
  EXPECT_TRUE(passed);
  return;
}  // test_3

TEST(nut_RNG,Buffer_RNG_double_rolled_over)
{
  bool passed(true);

  typedef double fp_t;
  size_t const vals_sz = 4;
  fp_t const in_vals[vals_sz] = {0., 0.11, 0.22, 0.33};

  typedef nut::Buffer_RNG<fp_t> rng_t;
  rng_t rng_base(&in_vals[0], vals_sz);
  rng_adaptor<fp_t, rng_t> rng(rng_base);

  std::vector<fp_t> out_vals(2 * vals_sz);

  std::generate(out_vals.begin(), out_vals.end(), rng);
  passed =
      passed and std::equal(&in_vals[0], &in_vals[vals_sz], out_vals.begin());
  passed = passed and std::equal(&in_vals[0], &in_vals[vals_sz],
                                 out_vals.begin() + vals_sz);
  EXPECT_TRUE(passed);
  return;
}  // test_4

TEST(nut_RNG,init_inst_LCG_RNG)
{
  bool passed(true);

  typedef nut::LCG_RNG rng_t;
  int32_t const seed(-120211);
  rng_t rng(seed);
  EXPECT_TRUE(passed);
  return;
}  // test_5

TEST(nut_RNG,use_LCG_RNG_print_a_few_values)
{
  bool passed(true);

  typedef nut::LCG_RNG rng_t;
  int32_t const seed(-120211);
  rng_t rng_base(seed);
  rng_adaptor<double, rng_t> rng(rng_base);

  size_t const vals_sz(4);

  std::vector<double> out_vals(vals_sz);

  std::generate(out_vals.begin(), out_vals.end(), rng);

  std::copy(out_vals.begin(), out_vals.end(),
            std::ostream_iterator<double>(std::cout, ","));
  std::cout << std::endl;
  EXPECT_TRUE(passed);
  return;
}  // test_6

namespace {
bool
lessThanEq0(double const d)
{
  return d <= 0.0;
}
bool
greaterThanEq1(double const d)
{
  return d >= 1.0;
}
std::string err_lte0("some values were < 0.0: FAIL");
std::string err_gte1("some values were > 1.0: FAIL");

/* check an iterator range against a unary predicate, print a
 * string if predicate fails anywhere. */
template <class it_t, class pred_t>
bool
checkPred(it_t first, it_t last, pred_t pred, std::string errstr)
{
  it_t it = find_if(first, last, pred);
  bool passed = it == last;
  if(!passed) { std::cerr << errstr << std::endl; }
  return passed;
}

}  // namespace

TEST(nut_RNG,use_LCG_RNG_draw_1e6_values_check_for_0_lt_x_lt_1)
{
  bool passed(true);

  typedef nut::LCG_RNG rng_t;
  int32_t const seed(-120211);
  rng_t rng_base(seed);
  rng_adaptor<double, rng_t> rng(rng_base);

  size_t const vals_sz(1000000);
  std::vector<double> out_vals(vals_sz);
  std::generate(out_vals.begin(), out_vals.end(), rng);

  // check each value greater than 0 and less than 1
  bool lt0_passed =
      checkPred(out_vals.begin(), out_vals.end(), lessThanEq0, err_lte0);
  passed = passed and lt0_passed;

  bool gt1_passed =
      checkPred(out_vals.begin(), out_vals.end(), greaterThanEq1, err_gte1);
  passed = passed and gt1_passed;
  EXPECT_TRUE(passed);
  return;
}  // test_7

TEST(nut_RNG,init_inst_MLCG)
{
  bool passed(true);

  typedef nut::MLCG rng_t;
  int32_t const seed(42);
  rng_t rng(seed);
  EXPECT_TRUE(passed);
  return;
}  // test_8

TEST(nut_RNG,use_MLCG_print_a_few_values_seed_42)
{
  bool passed(true);

  typedef nut::MLCG rng_t;
  int32_t const seed(42);
  rng_t rng_base(seed);
  rng_adaptor<double, rng_t> rng(rng_base);

  size_t const vals_sz(4);

  std::vector<double> out_vals(vals_sz);

  std::generate(out_vals.begin(), out_vals.end(), rng);

  std::copy(out_vals.begin(), out_vals.end(),
            std::ostream_iterator<double>(std::cout, ","));
  std::cout << std::endl;
  EXPECT_TRUE(passed);
  return;
}  // test_9

TEST(nut_RNG,use_MLCG_draw_1e6_values_check_for_0_lt_x_lt_1)
{
  bool passed(true);

  typedef nut::MLCG rng_t;
  int32_t const seed(42);
  rng_t rng_base(seed);
  rng_adaptor<double, rng_t> rng(rng_base);

  size_t const vals_sz(1000000);
  std::vector<double> out_vals(vals_sz);
  std::generate(out_vals.begin(), out_vals.end(), rng);

  // check each value greater than 0 and less than 1
  bool lt0_passed =
      checkPred(out_vals.begin(), out_vals.end(), lessThanEq0, err_lte0);
  passed = passed and lt0_passed;

  bool gt1_passed =
      checkPred(out_vals.begin(), out_vals.end(), greaterThanEq1, err_gte1);
  passed = passed and gt1_passed;
  EXPECT_TRUE(passed);
  return;
}  // test_10

TEST(nut_RNG,use_MLCG_test_splitting)
{
  bool passed(true);

  typedef nut::MLCG rng_t;
  int32_t const seed(42);
  rng_t rng_base1(seed);
  rng_adaptor<double, rng_t> rng1(rng_base1);

  nut::MLCG::new_gens ngens = rng_base1.split();
  rng_adaptor<double, rng_t> rng2(ngens.first);
  rng_adaptor<double, rng_t> rng3(ngens.second);

  // check old generator
  {
    std::cout << "after split, state of RNG 1: ";
    rng1.dump_state(std::cout);
    size_t const vals_sz(4);
    std::vector<double> out_vals(vals_sz);
    std::generate(out_vals.begin(), out_vals.end(), rng1);
    std::cout << "; some values: ";
    std::copy(out_vals.begin(), out_vals.end(),
              std::ostream_iterator<double>(std::cout, ","));
    std::cout << std::endl;
  }
  // check new generators
  {
    std::cout << "after split, state of RNG 2: ";
    rng2.dump_state(std::cout);
    size_t const vals_sz(4);
    std::vector<double> out_vals(vals_sz);
    std::generate(out_vals.begin(), out_vals.end(), rng2);
    std::cout << "; some values: ";
    std::copy(out_vals.begin(), out_vals.end(),
              std::ostream_iterator<double>(std::cout, ","));
    std::cout << std::endl;
  }

  {
    std::cout << "after split, state of RNG 3: ";
    rng3.dump_state(std::cout);
    size_t const vals_sz(4);
    std::vector<double> out_vals(vals_sz);
    std::generate(out_vals.begin(), out_vals.end(), rng3);
    std::cout << "; some values: ";
    std::copy(out_vals.begin(), out_vals.end(),
              std::ostream_iterator<double>(std::cout, ","));
    std::cout << std::endl;
  }
  EXPECT_TRUE(passed);
  return;
}  // test_11

TEST(nut_RNG, init_inst_Philox4x32)
{
  bool passed(true);

  nut::Philox4x32_RNG::ctr_t c;
  nut::Philox4x32_RNG::key_t k;
  c[0] = c[1] = c[2] = 0;
  c[3] = 1;
  k[0] = 0xdeadbeef;
  k[1] = 0xcafecafe;

  nut::Philox4x32_RNG g(c, k);
  EXPECT_TRUE(passed);
  return;
}  // test_12

TEST(nut_RNG,print_a_few_values)
{
  bool passed(true);

  typedef nut::Philox4x32_RNG rng_t;
  nut::Philox4x32_RNG::ctr_t c;
  nut::Philox4x32_RNG::key_t k;
  c[0] = c[1] = c[2] = 0;
  c[3] = 1;
  k[0] = 0xdeadbeef;
  k[1] = 0xcafecafe;
  rng_t rng_base(c, k);
  rng_adaptor<double, rng_t> rng(rng_base);

  size_t const vals_sz(10);

  std::vector<double> out_vals(vals_sz);

  std::generate(out_vals.begin(), out_vals.end(), rng);

  std::copy(out_vals.begin(), out_vals.end(),
            std::ostream_iterator<double>(std::cout, ","));
  std::cout << std::endl;
  EXPECT_TRUE(passed);
  return;
}  // test_13

TEST(nut_RNG,compatible_with_Haskell_McPhD_usage)
{
  bool passed(true);

  typedef nut::Philox4x32_RNG rng_t;
  nut::Philox4x32_RNG::ctr_t c;
  nut::Philox4x32_RNG::key_t k;
  c[0] = c[1] = c[2] = 0;
  c[3] = 1;
  k[0] = 0xdeadbeef;
  k[1] = 0xcafecafe;
  rng_t rng_base(c, k);
  rng_adaptor<double, rng_t> rng(rng_base);

  size_t const vals_sz(12);

  std::vector<double> out_vals(vals_sz);

  std::generate(out_vals.begin(), out_vals.end(), rng);

  // These values generated via the function randoms in the
  // Philo2 module of McPhD/basic, using same counter & key
  // as this test.
  double evals_arr[12] = {
      0.9848038759115404,  0.8800535828289092,  0.8015360604392531,
      0.2025851484105382,  0.27783947633403117, 0.634825162264666,
      0.49563039019191557, 0.769547123123007,   0.7025144994130103,
      0.6585216467165459,  0.677582395892678,   0.8577909275075531};
  std::vector<double> exp_vals(&evals_arr[0], &evals_arr[12]);

  passed =
      passed and std::equal(out_vals.begin(), out_vals.end(), exp_vals.begin());

  // std::copy(out_vals.begin(),out_vals.end(),
  //           std::ostream_iterator<double>(std::cout,","));
  // std::cout << std::endl;
  EXPECT_TRUE(passed);
  return;
}  // test_14

// End of file
