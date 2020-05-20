// T. M. Kelley (c) 2011 LANS LLC

#include "Cell_Data.hh"
#include "detail/Vector.h"
#include "gtest/gtest.h"
#include "test_aux.hh"

using fp_t = float;
using vec3 = nut::Vector;
using Cell_Data_T = nut::Cell_Data<fp_t, vec3>;
using nut::Cell_Data_Iterator;
auto fetch_l_nue = [](Cell_Data_T const & d) { return d.l_nue; };
auto fetch_l_nuebar = [](Cell_Data_T const & d) { return d.l_nueb; };
auto fetch_l_nux = [](Cell_Data_T const & d) { return d.l_nux; };
using NuE_Data_It = Cell_Data_Iterator<Cell_Data_T, decltype(fetch_l_nue)>;
using NuEBar_Data_It =
    Cell_Data_Iterator<Cell_Data_T, decltype(fetch_l_nuebar)>;
using NuX_Data_It = Cell_Data_Iterator<Cell_Data_T, decltype(fetch_l_nux)>;
using vector_data[[maybe_unused]] = std::vector<Cell_Data_T>;

TEST(nut_Cell_Data, instantiate)
{
  Cell_Data_T d1{};
  EXPECT_EQ(d1.rho_p, fp_t(0.0));
  EXPECT_EQ(d1.rho_p, fp_t(0.0));
  EXPECT_EQ(d1.rho_e_minus, fp_t(0.0));
  EXPECT_EQ(d1.rho_e_plus, fp_t(0.0));
  EXPECT_EQ(d1.rho_A, fp_t(0.0));
  EXPECT_EQ(d1.y_e, fp_t(.0));
  EXPECT_EQ(d1.abar, fp_t(0.0));
  EXPECT_EQ(d1.T_p, fp_t(0.0));
  EXPECT_EQ(d1.T_e_minus, fp_t(0.0));
  EXPECT_EQ(d1.T_e_plus, fp_t(0.0));
  EXPECT_EQ(d1.l_nue, fp_t(0.0));
  EXPECT_EQ(d1.l_nueb, fp_t(0.0));
  EXPECT_EQ(d1.l_nux, fp_t(0.0));
  EXPECT_EQ(d1.velocity, (vec3{0.0, 0.0, 0.0}));
  return;
}  // TEST(nut_Cell_Data, instantiate)

TEST(nut_Cell_Data, vector_cell_data)
{
  vector_data data({{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, {1, 1, 1}},
                    {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, {2, 2, 2}},
                    {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, {3, 3, 3}}});

  fp_t f{1.0};
  for(auto const & d : data) {
    EXPECT_EQ(d.rho_p, fp_t(f));
    EXPECT_EQ(d.rho_p, fp_t(f));
    EXPECT_EQ(d.rho_e_minus, fp_t(f));
    EXPECT_EQ(d.rho_e_plus, fp_t(f));
    EXPECT_EQ(d.rho_A, fp_t(f));
    EXPECT_EQ(d.y_e, fp_t(f));
    EXPECT_EQ(d.abar, fp_t(f));
    EXPECT_EQ(d.T_p, fp_t(f));
    EXPECT_EQ(d.T_e_minus, fp_t(f));
    EXPECT_EQ(d.T_e_plus, fp_t(f));
    EXPECT_EQ(d.l_nue, fp_t(f));
    EXPECT_EQ(d.l_nueb, fp_t(f));
    EXPECT_EQ(d.l_nux, fp_t(f));
    EXPECT_EQ(d.velocity, (vec3{f, f, f}));
    f += 1.0;
  }
}  // TEST(nut_Cell_Data, Data_Iterator)

TEST(nut_Cell_Data, Data_Iterator)
{
  vector_data data({{1, 1, 1, 1, 1, 1, 1, 1, 1, 1.1, 1.2, 1.3, {1, 1, 1}},
                    {2, 2, 2, 2, 2, 2, 2, 2, 2, 2.1, 2.2, 2.3, {2, 2, 2}},
                    {3, 3, 3, 3, 3, 3, 3, 3, 3, 3.1, 3.2, 3.3, {3, 3, 3}}});
  {
    fp_t f{1.1};
    NuE_Data_It ne{data.cbegin(), fetch_l_nue};
    NuE_Data_It ne_end{data.cend(), fetch_l_nue};
    for(; ne != ne_end; ne++) {
      EXPECT_EQ(*ne, f);
      f += 1.0;
    }
  }
  {
    fp_t f{1.2};
    NuEBar_Data_It neb{data.cbegin(),fetch_l_nuebar};
    NuEBar_Data_It neb_end{data.cend(),fetch_l_nuebar};
    for(; neb != neb_end; neb++) {
      EXPECT_EQ(*neb, f);
      f += 1.0;
    }
  }
  {
    fp_t f{1.3};
    NuX_Data_It nx{data.cbegin(), fetch_l_nux};
    NuX_Data_It nx_end{data.cend(), fetch_l_nux};
    for(; nx != nx_end; nx++) {
      EXPECT_EQ(*nx, f);
      f += 1.0;
    }
  }
}  // TEST(nut_Cell_Data, Data_Iterator)

// End of file
