// test_Vector.cc
// Apr 08, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#include "base/soft_equiv.h"
#include "base/test_common.h"
#include "mesh_common/Vector.h"
#include "gtest/gtest.h"
#include <cmath>

using nut_mesh::soft_equiv_os;
using nut_mesh::Vector;
using nut_mesh::Vector1;
using Vector2 = nut_mesh::Vector2<double>;

TEST(default_mesh_Vector, instantiate) {
  Vector v1{1.0, 2.0, 3.0};
  v1[0] = 1.0;
  EXPECT_TRUE(true);
}

TEST(default_mesh_Vector, norm) {
  Vector v1{1.0, 2.0, 3.0};
  double const n1 = v1.norm();
  double const exp_n1 = sqrt(14.0);
  bool const n1_ok = soft_equiv_os(n1, exp_n1, "v1 norm");
  EXPECT_TRUE(n1_ok);
} // TEST(default_mesh_Vector,norm){

constexpr size_t n_test_vectors{100u};
extern Vector test_vectors_1[n_test_vectors];
extern Vector test_vectors_2[n_test_vectors];
/* Components of v1s parallel to counterpart in v2s*/
extern Vector test_vectors_par[n_test_vectors];
/* Components of v1s perpendicular to counterpart in v2s*/
extern Vector test_vectors_perp[n_test_vectors];
/* Dot products v1 . v2 */
extern double test_vectors_dots[n_test_vectors];
/* Some unit vectors */
extern Vector test_unit_vectors[n_test_vectors];

TEST(default_mesh_Vector, operator_equal) {
  Vector v1{1.0, 2.0, 3.0};
  Vector v2{1.0, 2.0, 3.0};
  EXPECT_TRUE(v1 == v2);
  EXPECT_TRUE(v2 == v1);
  EXPECT_FALSE(v1 != v2);
  EXPECT_FALSE(v2 != v1);

  Vector v3{2.0, 3.0, 4.0};
  EXPECT_TRUE(v1 != v3);
  EXPECT_TRUE(v3 != v2);

  EXPECT_TRUE(v1 == v1);
  EXPECT_TRUE(v2 == v2);
  EXPECT_TRUE(v3 == v3);
} // TEST(default_mesh_Vector,operator_equal){

TEST(default_mesh_Vector, dot) {
  Vector v1{1.0, 0.0, 0.0};
  Vector v2{1.0, 0.0, 0.0};
  double const v1_dot_v2 = v1.dot(v2);

  EXPECT_EQ(v1_dot_v2, 1.0);

  Vector v3{0.0, 1.0, 0.0};
  double const v1_dot_v3{v3.dot(v1)};
  EXPECT_EQ(v1_dot_v3, 0.0);

  Vector v4{0.0, 0.0, 1.0};
  double const v3_dot_v4{v4.dot(v3)};
  EXPECT_EQ(v3_dot_v4, 0.0);

  // many test cases
  for (size_t i = 0; i < n_test_vectors; ++i) {
    std::string const label(std::to_string(i) + " ");
    Vector const &v1 = test_vectors_1[i];
    Vector const &v2 = test_vectors_2[i];
    double dot_exp = test_vectors_dots[i];
    double dprod = v1.dot(v2);
    bool dot_ok = soft_equiv_os(dprod, dot_exp,
                                ("dot product test case " + label), 5e-14);
    EXPECT_TRUE(dot_ok);
    // dot product commutes
    double dp2 = v2.dot(v1);
    bool commute_ok = soft_equiv_os(dprod, dp2, "dot commutes case " + label);
    EXPECT_TRUE(commute_ok);
  }
} // TEST(default_mesh_Vector,dot)

TEST(default_mesh_Vector, scale) {
  Vector v1{1.0, 2.0, 3.0};
  Vector v2{0.0, 0.0, 0.0};

  Vector v1_s = v1.scale(0.0);
  EXPECT_TRUE(v1 != v1_s);
  EXPECT_TRUE(v1_s == v2);

  v1.scale_this(2.0);
  Vector v1_exp{2.0, 4.0, 6.0};
  EXPECT_TRUE(v1 == v1_exp);
} // TEST(default_mesh_Vector,scale)

TEST(default_mesh_Vector, operator_mult) {
  Vector const v1{1.0, 2.0, 3.0};

  Vector v2 = v1 * 3.0;
  Vector v2_exp{3.0, 6.0, 9.0};
  EXPECT_TRUE(v2 == v2_exp);

  Vector v3 = 3.0 * v1;
  Vector v3_exp{3.0, 6.0, 9.0};
  EXPECT_TRUE(v3 == v3_exp);
}

TEST(default_mesh_Vector, operator_index) {
  Vector const v1{1.0, 2.0, 3.0};

  EXPECT_EQ(1.0, v1[0]);
  EXPECT_EQ(2.0, v1[1]);
  EXPECT_EQ(3.0, v1[2]);

  Vector v2{1.0, 2.0, 3.0};

  EXPECT_EQ(1.0, v2[0]);
  EXPECT_EQ(2.0, v2[1]);
  EXPECT_EQ(3.0, v2[2]);

  v2[0] += 1.0;
  v2[1] += 1.0;
  v2[2] += 1.0;

  Vector v2_exp{2.0, 3.0, 4.0};
  EXPECT_TRUE(v2 == v2_exp);
}

TEST(default_mesh_Vector, parallel) {
  // some easy to consider cases
  {
    Vector const v1{1.0, 0.0, 0.0};
    Vector const v2{0.0, 1.0, 0.0};
    Vector const v2_p_v1{0.0, 0.0, 0.0};
    EXPECT_TRUE(v2_p_v1 == v1.parallel(v2));
  }
  {
    Vector const v1{1.0, 0.0, 0.0};
    Vector const v2{10.0, 0.0, 0.0};
    Vector const v1_p_v2{v1.parallel(v2)};
    EXPECT_TRUE(v1 == v1.parallel(v2));
    EXPECT_TRUE(vec_soft_equiv(v1, v1_p_v2, "parallel test case 1"));
  }
  {
    double const s2 = std::sqrt(2.0);
    Vector const v1{s2, 0.0, s2};
    Vector const v2{10.0, 0.0, 10.0};
    Vector const v1_p_v2{v1.parallel(v2)};
    vec_soft_equiv(v1, v1_p_v2, "parallel test case 2");
  }
  // ok now try a whole list of vectors
  for (size_t i = 0; i < n_test_vectors; ++i) {
    Vector const &v1 = test_vectors_1[i];
    Vector const &v2 = test_vectors_2[i];
    Vector const &exp_par = test_vectors_par[i];
    Vector const par = v1.parallel(v2);
    vec_soft_equiv(par, exp_par,
                   ("parallel auto test case " + std::to_string(i) + " "));
    // projection is idempotent
    Vector const par2 = par.parallel(v2);
    vec_soft_equiv(
        par2, par,
        ("parall. auto test case (idem) " + std::to_string(i) + " "));
    // The projected part has magnitude less than or equal to the original
    double par_norm = par.norm();
    double v1_norm = v1.norm();
    EXPECT_TRUE(par_norm <= v1_norm);
  }
  // EXPECT_TRUE(false); // Failed 20 June 2019
  return;
} // TEST(default_mesh_Vector, parallel)

TEST(default_mesh_Vector, perpendicular) {
  // some easy to consider cases
  {
    Vector const v1{1.0, 0.0, 0.0};
    Vector const v2{0.0, 1.0, 0.0};
    bool v1_p_v2_ok =
        vec_soft_equiv(v1, v1.perpendicular(v2), "perpendicular test case 1");
    EXPECT_TRUE(v1_p_v2_ok);
  }
  {
    Vector const v1{1.0, 0.0, 0.0};
    Vector const v2{10.0, 0.0, 0.0};
    Vector const v1_p_v2{v1.perpendicular(v2)};
    Vector const v1_p_v2_exp{0.0, 0.0, 0.0};
    EXPECT_TRUE(v1_p_v2 == v1_p_v2_exp);
  }
  for (size_t i = 0; i < n_test_vectors; ++i) {
    std::string label(std::to_string(i) + " ");
    Vector const &v1 = test_vectors_1[i];
    Vector const &v2 = test_vectors_2[i];
    Vector const &exp_perp = test_vectors_perp[i];
    Vector const perp = v1.perpendicular(v2);
    vec_soft_equiv(perp, exp_perp, ("perpendicular auto test case " + label));
    // perpendicular projection is idempotent
    Vector const perp2 = perp.perpendicular(v2);
    vec_soft_equiv(perp2, perp, ("perp. auto test case (idem) " + label));
    // The projected part has magnitude less than or equal to the original
    double perp_norm = perp.norm();
    double v1_norm = v1.norm();
    EXPECT_TRUE(perp_norm <= v1_norm);
    // The perpendicular part is perpendicular to the parallel part
    Vector const par = v1.parallel(v2);
    double const par_dot_perp{par.dot(perp)};
    bool p_dot_p =
        soft_equiv_os(par_dot_perp, 0.0, "parallel dot perp comps", 3e-10);
    EXPECT_TRUE(p_dot_p);
    // Adding up the perpendicular and parallel Components gives back original
    Vector const sum{perp + par};
    bool sum_ok = vec_soft_equiv(sum, v1, "sum of Components case " + label);
    EXPECT_TRUE(sum_ok);
  }
  return;
} // TEST(default_mesh_Vector, perpendicular)

TEST(default_mesh_Vector, reflect) {
  {
    double const s2 = std::sqrt(2.0);
    Vector const n1{s2, 0.0, s2};
    Vector const v1{1.0, 1.0, 1.0};
    Vector const r1{v1.reflect(n1)};
    Vector const r1_exp{-1.0, 1.0, -1.0};
    bool r1_ok = vec_soft_equiv(r1, r1_exp, "reflect: first case");
    EXPECT_TRUE(r1_ok);
  }
  for (size_t i = 0; i < n_test_vectors; ++i) {
    std::string label(std::to_string(i) + " ");
    Vector const &v = test_vectors_1[i];
    Vector const &n = test_unit_vectors[i];
    Vector const r = v.reflect(n);
    // parallel part of r should be negative of parallel part v
    Vector v_par = v.parallel(n);
    Vector r_par = r.parallel(n);
    r_par.scale_this(-1.0);
    bool par_ok =
        vec_soft_equiv(r_par, v_par, "reflect parall part auto case " + label);
    EXPECT_TRUE(par_ok);
    // perpendicular part of r should be same as perp. part v
    Vector v_perp = v.perpendicular(n);
    Vector r_perp = r.perpendicular(n);
    bool perp_ok =
        vec_soft_equiv(r_perp, v_perp, "reflect perp part auto case " + label);
    EXPECT_TRUE(perp_ok);
  }
  // EXPECT_TRUE(false); // Failed 20 June 2019
  return;
} // TEST(default_mesh_Vector, reflect)

// -------------------------------- 2D Tests -------------------------------- //
extern Vector2 test_2D_vectors_1[n_test_vectors];
extern Vector2 test_2D_vectors_2[n_test_vectors];
/* Components of v1s parallel to counterpart in v2s*/
extern Vector2 test_2D_vectors_par[n_test_vectors];
/* Components of v1s perpendicular to counterpart in v2s*/
extern Vector2 test_2D_vectors_perp[n_test_vectors];
/* Dot products v1 . v2 */
extern double test_2D_vectors_dots[n_test_vectors];
/* Some unit vectors */
extern Vector2 test_2D_unit_vectors[n_test_vectors];

TEST(default_mesh_Vector, operator_equal_2D) {
  Vector2 v1{1.0, 2.0};
  Vector2 v2{1.0, 2.0};
  EXPECT_TRUE(v1 == v2);
  EXPECT_TRUE(v2 == v1);
  EXPECT_FALSE(v1 != v2);
  EXPECT_FALSE(v2 != v1);

  Vector2 v3{2.0, 3.0};
  EXPECT_TRUE(v1 != v3);
  EXPECT_TRUE(v3 != v2);

  EXPECT_TRUE(v1 == v1);
  EXPECT_TRUE(v2 == v2);
  EXPECT_TRUE(v3 == v3);
} // TEST(default_mesh_Vector2,operator_equal){

TEST(default_mesh_Vector2, dot_2D) {
  Vector2 v1{1.0, 0.0};
  Vector2 v2{1.0, 0.0};
  double const v1_dot_v2 = v1.dot(v2);

  EXPECT_EQ(v1_dot_v2, 1.0);

  Vector2 v3{0.0, 1.0};
  double const v1_dot_v3{v3.dot(v1)};
  EXPECT_EQ(v1_dot_v3, 0.0);

  Vector2 v4{0.0, 0.0};
  double const v3_dot_v4{v4.dot(v3)};
  EXPECT_EQ(v3_dot_v4, 0.0);

  // many test cases
  for (size_t i = 0; i < n_test_vectors; ++i) {
    std::string const label(std::to_string(i) + " ");
    Vector2 const &v1 = test_2D_vectors_1[i];
    Vector2 const &v2 = test_2D_vectors_2[i];
    double dot_exp = test_2D_vectors_dots[i];
    double dprod = v1.dot(v2);
    bool dot_ok = soft_equiv_os(dprod, dot_exp,
                                ("2D dot product test case " + label), 1e-13);
    EXPECT_TRUE(dot_ok);
    // dot product commutes
    double dp2 = v2.dot(v1);
    bool commute_ok =
        soft_equiv_os(dprod, dp2, "2D dot commutes case " + label);
    EXPECT_TRUE(commute_ok);
  }
} // TEST(default_mesh_Vector2,dot)

TEST(default_mesh_Vector, parallel_2D) {
  // some easy to consider cases
  {
    Vector2 const v1{1.0, 0.0};
    Vector2 const v2{0.0, 1.0};
    Vector2 const v2_p_v1{0.0, 0.0};
    EXPECT_TRUE(v2_p_v1 == v1.parallel(v2));
  }
  {
    Vector2 const v1{1.0, 0.0};
    Vector2 const v2{10.0, 0.0};
    Vector2 const v1_p_v2{v1.parallel(v2)};
    EXPECT_TRUE(v1 == v1_p_v2);
  }
  {
    double const s2 = std::sqrt(2.0);
    Vector2 const v1{s2, s2};
    Vector2 const v2{10.0, 10.0};
    Vector2 const v1_p_v2{v1.parallel(v2)};
    vec_soft_equiv(v1, v1_p_v2, "parallel test case 2");
  }
  // ok now try a whole list of vectors
  for (size_t i = 0; i < n_test_vectors; ++i) {
    Vector2 const &v1 = test_2D_vectors_1[i];
    Vector2 const &v2 = test_2D_vectors_2[i];
    Vector2 const &exp_par = test_2D_vectors_par[i];
    Vector2 const par = v1.parallel(v2);
    vec_soft_equiv(par, exp_par,
                   ("2D parallel auto test case " + std::to_string(i) + " "),
                   1e-13);
    // projection is idempotent
    Vector2 const par2 = par.parallel(v2);
    vec_soft_equiv(
        par2, par,
        ("2D parall. auto test case (idem) " + std::to_string(i) + " "));
    // The projected part has magnitude less than or equal to the original
    double par_norm = par.norm();
    double v1_norm = v1.norm();
    EXPECT_TRUE(par_norm <= v1_norm);
  }
  // EXPECT_TRUE(false); // Failed 20 June 2019
  return;
} // TEST(default_mesh_Vector, parallel_2D)

TEST(default_mesh_Vector, perpendicular_2D) {
  // some easy to consider cases
  {
    Vector2 const v1{1.0, 0.0};
    Vector2 const v2{0.0, 1.0};
    bool v1_p_v2_ok =
        vec_soft_equiv(v1, v1.perpendicular(v2), "perpendicular test case 1");
    EXPECT_TRUE(v1_p_v2_ok);
  }
  {
    Vector2 const v1{1.0, 0.0};
    Vector2 const v2{10.0, 0.0};
    Vector2 const v1_p_v2{v1.perpendicular(v2)};
    Vector2 const v1_p_v2_exp{0.0, 0.0};
    EXPECT_TRUE(v1_p_v2 == v1_p_v2_exp);
  }
  for (size_t i = 0; i < n_test_vectors; ++i) {
    std::string label(std::to_string(i) + " ");
    Vector2 const &v1 = test_2D_vectors_1[i];
    Vector2 const &v2 = test_2D_vectors_2[i];
    Vector2 const &exp_perp = test_2D_vectors_perp[i];
    Vector2 const perp = v1.perpendicular(v2);
    vec_soft_equiv(perp, exp_perp, ("2D perpendicular auto test case " + label),
                   1e-12);
    // perpendicular projection is idempotent
    Vector2 const perp2 = perp.perpendicular(v2);
    vec_soft_equiv(perp2, perp, ("2D perp. auto test case (idem) " + label),
                   1e-12);
    // The projected part has magnitude less than or equal to the original
    double perp_norm = perp.norm();
    double v1_norm = v1.norm();
    EXPECT_TRUE(perp_norm <= v1_norm);
    // The perpendicular part is perpendicular to the parallel part
    Vector2 const par = v1.parallel(v2);
    double const par_dot_perp{par.dot(perp)};
    bool p_dot_p =
        soft_equiv_os(par_dot_perp, 0.0, "parallel dot perp comps", 3e-10);
    EXPECT_TRUE(p_dot_p);
    // Adding up the perpendicular and parallel Components gives back original
    Vector2 const sum{perp + par};
    bool sum_ok = vec_soft_equiv(sum, v1, "sum of Components case " + label);
    EXPECT_TRUE(sum_ok);
  }
  return;
} // TEST(default_mesh_Vector2, perpendicular_2D)

TEST(default_mesh_Vector2, reflect_2D) {
  {
    double const s2 = std::sqrt(2.0);
    Vector2 const n1{s2, s2};
    Vector2 const v1{1.0, 1.0};
    Vector2 const r1{v1.reflect(n1)};
    Vector2 const r1_exp{-1.0, -1.0};
    bool r1_ok = vec_soft_equiv(r1, r1_exp, "reflect: first case");
    EXPECT_TRUE(r1_ok);
  }
  for (size_t i = 0; i < n_test_vectors; ++i) {
    std::string label(std::to_string(i) + " ");
    Vector2 const &v = test_2D_vectors_1[i];
    Vector2 const &n = test_2D_unit_vectors[i];
    Vector2 const r = v.reflect(n);
    // parallel part of r should be negative of parallel part v
    Vector2 v_par = v.parallel(n);
    Vector2 r_par = r.parallel(n);
    r_par.scale_this(-1.0);
    bool par_ok = vec_soft_equiv(r_par, v_par,
                                 "2D reflect parall part auto case " + label);
    EXPECT_TRUE(par_ok);
    // perpendicular part of r should be same as perp. part v
    Vector2 v_perp = v.perpendicular(n);
    Vector2 r_perp = r.perpendicular(n);
    bool perp_ok = vec_soft_equiv(
        r_perp, v_perp, "2D reflect perp part auto case " + label, 2e-12);
    EXPECT_TRUE(perp_ok);
  }
  // EXPECT_TRUE(false); // Failed 20 June 2019
  return;
} // TEST(default_mesh_Vector, reflect_2D)

// -------------------------------- 1D Tests -------------------------------- //
extern Vector1 test_1D_vectors_1[n_test_vectors];
extern Vector1 test_1D_vectors_2[n_test_vectors];
/* Components of v1s parallel to counterpart in v2s*/
extern Vector1 test_1D_vectors_par[n_test_vectors];
/* Components of v1s perpendicular to counterpart in v2s*/
extern Vector1 test_1D_vectors_perp[n_test_vectors];
/* Dot products v1 . v2 */
extern double test_1D_vectors_dots[n_test_vectors];
/* Some unit vectors */
extern Vector1 test_1D_unit_vectors[n_test_vectors];

TEST(default_mesh_Vector, operator_equal_1D) {
  Vector1 v1{1.0};
  Vector1 v2{1.0};
  EXPECT_TRUE(v1 == v2);
  EXPECT_TRUE(v2 == v1);
  EXPECT_FALSE(v1 != v2);
  EXPECT_FALSE(v2 != v1);

  Vector1 v3{2.0};
  EXPECT_TRUE(v1 != v3);
  EXPECT_TRUE(v3 != v2);

  EXPECT_TRUE(v1 == v1);
  EXPECT_TRUE(v2 == v2);
  EXPECT_TRUE(v3 == v3);
} // TEST(default_mesh_Vector1,operator_equal){

TEST(default_mesh_Vector1, dot_1D) {
  Vector1 v1{1.0};
  Vector1 v2{1.0};
  double const v1_dot_v2 = v1.dot(v2);

  EXPECT_EQ(v1_dot_v2, 1.0);

  Vector1 v3{0.0};
  double const v1_dot_v3{v3.dot(v1)};
  EXPECT_EQ(v1_dot_v3, 0.0);

  Vector1 v4{0.0};
  double const v3_dot_v4{v4.dot(v3)};
  EXPECT_EQ(v3_dot_v4, 0.0);

  // many test cases
  for (size_t i = 0; i < n_test_vectors; ++i) {
    std::string const label(std::to_string(i) + " ");
    Vector1 const &v1 = test_1D_vectors_1[i];
    Vector1 const &v2 = test_1D_vectors_2[i];
    double dot_exp = test_1D_vectors_dots[i];
    double dprod = v1.dot(v2);
    bool dot_ok = soft_equiv_os(dprod, dot_exp,
                                ("1D dot product test case " + label), 1e-13);
    EXPECT_TRUE(dot_ok);
    // dot product commutes
    double dp2 = v2.dot(v1);
    bool commute_ok =
        soft_equiv_os(dprod, dp2, "1D dot commutes case " + label);
    EXPECT_TRUE(commute_ok);
  }
} // TEST(default_mesh_Vector1,dot)

TEST(default_mesh_Vector, parallel_1D) {
  // some easy to consider cases
  {
    Vector1 const v1{1.0};
    Vector1 const v2{0.0};
    Vector1 const v2_p_v1{0.0};
    EXPECT_TRUE(v2_p_v1 == v1.parallel(v2));
  }
  {
    Vector1 const v1{1.0};
    Vector1 const v2{10.0};
    Vector1 const v1_p_v2{v1.parallel(v2)};
    EXPECT_TRUE(vec_soft_equiv(v1, v1_p_v2, "parallel test case 1"));
  }
  {
    double const s2 = std::sqrt(2.0);
    Vector1 const v1{s2};
    Vector1 const v2{10.0};
    Vector1 const v1_p_v2{v1.parallel(v2)};
    vec_soft_equiv(v1, v1_p_v2, "parallel test case 2");
  }
  // ok now try a whole list of vectors
  for (size_t i = 0; i < n_test_vectors; ++i) {
    Vector1 const &v1 = test_1D_vectors_1[i];
    Vector1 const &v2 = test_1D_vectors_2[i];
    Vector1 const &exp_par = v1;
    Vector1 const par = v1.parallel(v2);
    vec_soft_equiv(par, exp_par,
                   ("1D parallel auto test case " + std::to_string(i) + " "),
                   1e-13);
    // projection is idempotent
    Vector1 const par2 = par.parallel(v2);
    vec_soft_equiv(
        par2, par,
        ("1D parall. auto test case (idem) " + std::to_string(i) + " "));
    // The projected part has magnitude less than or equal to the original
    double par_norm = par.norm();
    double v1_norm = v1.norm();
    EXPECT_TRUE(par_norm == v1_norm);
  }
  // EXPECT_TRUE(false); // Failed 20 June 2019
  return;
} // TEST(default_mesh_Vector, parallel_1D)

TEST(default_mesh_Vector, perpendicular_1D) {
  // In 1D, there is no perpendicular!
  {
    Vector1 const v1{1.0};
    Vector1 const v2{20.0};
    Vector1 perp_exp{0.0};
    bool v1_p_v2_ok = vec_soft_equiv(perp_exp, v1.perpendicular(v2),
                                     "perpendicular test case 1");
    EXPECT_TRUE(v1_p_v2_ok);
  }
  {
    Vector1 const v1{1.0};
    Vector1 const v2{10.0};
    Vector1 const v1_p_v2{v1.perpendicular(v2)};
    Vector1 const v1_p_v2_exp{0.0};
    EXPECT_TRUE(v1_p_v2 == v1_p_v2_exp);
  }
  for (size_t i = 0; i < n_test_vectors; ++i) {
    std::string label(std::to_string(i) + " ");
    Vector1 const &v1 = test_1D_vectors_1[i];
    Vector1 const &v2 = test_1D_vectors_2[i];
    Vector1 exp_perp{0.0};
    Vector1 const perp = v1.perpendicular(v2);
    vec_soft_equiv(perp, exp_perp, ("1D perpendicular auto test case " + label),
                   1e-12);
    // perpendicular projection is idempotent
    Vector1 const perp2 = perp.perpendicular(v2);
    vec_soft_equiv(perp2, perp, ("1D perp. auto test case (idem) " + label),
                   1e-12);
    // The projected part has magnitude less than or equal to the original
    double perp_norm = perp.norm();
    double v1_norm = v1.norm();
    EXPECT_TRUE(perp_norm <= v1_norm);
    // The perpendicular part is perpendicular to the parallel part
    Vector1 const par = v1.parallel(v2);
    double const par_dot_perp{par.dot(perp)};
    bool p_dot_p =
        soft_equiv_os(par_dot_perp, 0.0, "parallel dot perp comps", 3e-10);
    EXPECT_TRUE(p_dot_p);
    // Adding up the perpendicular and parallel Components gives back original
    Vector1 const sum{perp + par};
    bool sum_ok = vec_soft_equiv(sum, v1, "sum of Components case " + label);
    EXPECT_TRUE(sum_ok);
  }
  return;
} // TEST(default_mesh_Vector1, perpendicular_1D)

TEST(default_mesh_Vector1, reflect_1D) {
  {
    double const s2 = std::sqrt(2.0);
    Vector1 const n1{s2};
    Vector1 const v1{1.0};
    Vector1 const r1{v1.reflect(n1)};
    Vector1 const r1_exp{-1.0};
    bool r1_ok = vec_soft_equiv(r1, r1_exp, "reflect: first case");
    EXPECT_TRUE(r1_ok);
  }
  for (size_t i = 0; i < n_test_vectors; ++i) {
    std::string label(std::to_string(i) + " ");
    Vector1 const &v = test_1D_vectors_1[i];
    Vector1 const &n = test_1D_unit_vectors[i];
    Vector1 const r = v.reflect(n);
    // parallel part of r should be negative of parallel part v
    Vector1 v_par = v.parallel(n);
    Vector1 r_par = r.parallel(n);
    r_par.scale_this(-1.0);
    bool par_ok = vec_soft_equiv(r_par, v_par,
                                 "1D reflect parall part auto case " + label);
    EXPECT_TRUE(par_ok);
    // perpendicular part of r should be same as perp. part v
    Vector1 v_perp = v.perpendicular(n);
    Vector1 r_perp = r.perpendicular(n);
    bool perp_ok = vec_soft_equiv(
        r_perp, v_perp, "1D reflect perp part auto case " + label, 2e-12);
    EXPECT_TRUE(perp_ok);
  }
  // EXPECT_TRUE(false); // Failed 20 June 2019
  return;
} // TEST(default_mesh_Vector, reflect_1D)

// -------------------------------- Test Data ------------------------------- //

Vector1 test_1D_vectors_1[n_test_vectors] = {
    {756.408339468241},  {-487.89271010611},  {-699.643954386022},
    {449.201555111982},  {-498.388673145683}, {-16.5169641129119},
    {-234.382913433274}, {448.919888803649},  {-376.407811787764},
    {252.442758591075},  {-248.169235248285}, {64.8487424803307},
    {40.1258835599051},  {-608.512685025808}, {-306.738139450392},
    {-921.203799301609}, {-994.435995871227}, {-683.460637491271},
    {285.70312006108},   {-17.1954731549904}, {-504.018906671527},
    {-140.653972326683}, {-665.90918717602},  {299.594211057293},
    {511.37585924306},   {82.3292290112322},  {-948.670493028952},
    {500.724975692681},  {622.607956467413},  {-623.87525818653},
    {338.824200438727},  {121.700156791775},  {313.656294035737},
    {969.576634415837},  {-574.724714841841}, {-133.03498770352},
    {101.55545926998},   {40.168423361959},   {770.96965668928},
    {918.208565232925},  {-235.549326102525}, {189.991907236797},
    {122.636017554265},  {-509.703197121361}, {222.697435014578},
    {283.010556134157},  {-13.6438095194176}, {-941.753705451955},
    {-82.7400173457218}, {-462.19025311781},  {-131.573366518108},
    {-248.291761541405}, {768.568825457826},  {644.622098468165},
    {631.25716945108},   {-31.7309422197377}, {433.431054210393},
    {-868.733891097583}, {82.5256384049294},  {217.055960967446},
    {-38.9013045199836}, {540.29041163352},   {-371.584159806537},
    {-222.317566371587}, {427.339076909685},  {-321.979287500611},
    {-869.008736126535}, {276.280556725122},  {621.870681141589},
    {-4.29922240813175}, {599.272032338263},  {-656.958631130266},
    {-541.358146005631}, {-585.048528153321}, {-110.35859519752},
    {235.920311001452},  {316.974871983451},  {496.963520115075},
    {-842.45991424019},  {-594.543441808235}, {-136.684469023877},
    {741.593066984948},  {-379.030384769264}, {501.608478547929},
    {7.6985409432541},   {6.1149126934788},   {-108.738825524657},
    {-950.941080465838}, {256.250338073494},  {-26.2907499538915},
    {-620.375650347693}, {42.186162235505},   {-577.814053255571},
    {-784.996199007922}, {-186.34205315812},  {166.515333742088},
    {-595.594503495783}, {-994.201116468102}, {-90.3941044848812},
    {832.99177473809}};
Vector1 test_1D_vectors_2[n_test_vectors] = {
    {-865.916888715532}, {-630.347407573147}, {32.8357334534835},
    {940.46353642121},   {-15.6060578842607}, {953.325883787598},
    {-183.091613608832}, {835.261400259813},  {497.118182344496},
    {171.031222479174},  {334.619246917729},  {0.475403295361502},
    {-36.8349735696697}, {179.903884593837},  {932.386568102157},
    {-127.588850078247}, {702.694861619576},  {-655.524384994853},
    {472.422950666534},  {-263.100478477833}, {651.687062608361},
    {707.176225529272},  {-549.414096884388}, {828.199609761307},
    {386.470505109478},  {-42.1051932850387}, {-538.072349674564},
    {511.435273028844},  {-399.394304523251}, {-108.56612212997},
    {-14.8991743897491}, {-651.157125246177}, {-563.678700014093},
    {762.883996782358},  {-434.500263246387}, {-635.78344850993},
    {-962.233503819653}, {-480.344802038146}, {554.790753416099},
    {-370.912563176076}, {379.379430932625},  {-516.06921955811},
    {-817.422288002159}, {-787.147013767842}, {-334.889106149886},
    {310.296882996297},  {-853.04694708944},  {-614.472856263748},
    {-362.23710975297},  {-395.718143857325}, {440.742756459908},
    {977.844873020445},  {521.378921795897},  {841.503054242494},
    {348.52165734001},   {-330.808530609395}, {-109.763299663103},
    {305.326510143084},  {-484.357861260329}, {-925.618261804323},
    {364.154694719317},  {-54.1626638574776}, {552.573566266948},
    {951.655492985561},  {-923.414565956685}, {724.005117765958},
    {-613.311540055974}, {256.369197860959},  {-664.462995161193},
    {-334.148561044835}, {-269.066762245073}, {-179.925531510341},
    {335.505978491517},  {-174.10231727594},  {554.917433561387},
    {-544.295069529911}, {739.162757874113},  {592.409048758118},
    {-437.035058995651}, {279.008235517435},  {-791.146561546308},
    {253.902201243885},  {866.040276054544},  {-940.310095258853},
    {-344.484773888526}, {386.621540356179},  {-891.892160734755},
    {822.836542448066},  {-298.342979147681}, {352.836337346198},
    {921.385244729557},  {-453.93859680077},  {-309.344282639616},
    {-444.137396430177}, {151.634603064383},  {-901.528732627801},
    {-297.084831450834}, {342.194994773275},  {-606.875658022781},
    {-24.0508885716999}};

double test_1D_vectors_dots[n_test_vectors] = {
    -654986.755910821, 307541.904989224,  -22973.3223985606, 422457.683086522,
    7777.88248197144,  -15746.0494104297, 42913.5458228373,  374965.454926615,
    -187119.167216202, 43175.5936078466,  -83042.2026069298, 30.8293058751986,
    -1478.03586038875, -109473.795860769, -285998.521148192, 117535.333440604,
    -698785.064508256, 448025.114059656,  134972.71099389,   4524.13721473071,
    -328462.600787845, -99467.1452556822, 365859.89467933,   248123.808684397,
    197631.686622459,  -3466.48810052615, 510453.361251015,  256088.414655747,
    -248666.071763945, 67731.7174741455,  -5048.20084980389, -79245.9242385409,
    -176801.372073302, 739674.498049941,  249718.039892984,  84581.4432546199,
    -97720.0654053674, -19294.6933679846, 427726.836695597,  -340575.092460771,
    -89362.5692933392, -98048.9752900504, -100245.41406068,  401211.349522,
    -74578.9449539043, 87817.2934234772,  11638.810057209,   578682.089286032,
    29971.5047442249,  182897.069072727,  -57990.0082359008, -242790.826036478,
    400715.58554314,   542451.464693167,  220006.794904854,  10496.8663705631,
    -47574.8226865901, -265247.487211847, -39971.9417169549, -200910.961304954,
    -14166.0926716578, -29263.5679507246, -205327.584352605, -211569.733224702,
    -394611.128220888, -233114.651965079, 532973.08627586,   70829.8247121985,
    -413210.055394271, 1436.57898128894,  -161244.185445281, 118203.630886419,
    -181628.894489973, 101858.304470371,  -61239.9084184478, -128410.262080054,
    234296.020552082,  294405.686218857,  368184.518321432,  -165882.516637378,
    108137.447685023,  188292.112134682,  -328255.579058633, -471667.516246052,
    -2652.03013610845, 2364.15696469632,  96983.3060529458,  -782469.070722338,
    -76450.4892684466, -9276.33191981578, -571604.970419867, -19149.9272895947,
    178743.473803433,  348646.168034964,  -28255.9032648336, -150118.3577916,
    176942.092684088,  -340210.645853386, 54857.9816406423,  -20034.1923553684};
Vector1 test_1D_unit_vectors[n_test_vectors] = {};

Vector2 test_2D_vectors_1[n_test_vectors] = {
    {-624.062052406581, -631.68950884771},
    {955.253646606514, -269.523170595603},
    {-49.1463935426677, 741.984787722734},
    {-490.798145217042, 848.150576814166},
    {-121.720009635375, 990.329064035879},
    {-525.345113196958, -717.870748398445},
    {-266.368864082288, -154.41162183654},
    {-294.395915574397, 308.065543698056},
    {592.836869146168, 231.024699681672},
    {860.369412268814, -913.957897529725},
    {360.384516499546, 395.722600067216},
    {931.217310389738, 630.50086988515},
    {-775.874207644894, 475.07359575873},
    {767.361924412131, 704.291673597973},
    {300.274174585995, 833.588031645531},
    {-180.711080082616, -501.522210558596},
    {998.880061667901, 384.85572081783},
    {947.229120971249, 907.477900425275},
    {814.014599702244, -61.7876325263114},
    {839.398416653648, -376.183772970336},
    {335.316405273967, 651.545494096984},
    {-892.903802151338, 466.338215217326},
    {6.628438971793, 916.66590423835},
    {-257.964551794032, 907.425257794054},
    {498.518732742258, 852.625762341752},
    {69.2953160777902, -716.182982949217},
    {-288.721630702951, -799.968671386819},
    {395.915560445481, -884.737921511329},
    {-832.661381896589, 19.8753561807907},
    {-663.826722432581, -947.936443616226},
    {-328.144475010379, -578.501315926875},
    {506.652504100696, -807.168223525976},
    {-685.913410249804, -923.730833405429},
    {-886.124889296288, -608.182477738435},
    {-812.852707687717, 185.235233491155},
    {680.214937357026, -783.685943151783},
    {839.280881825171, 641.143610944947},
    {776.312327604101, -401.983108082762},
    {887.990396622385, 373.825909525649},
    {117.257648804507, 607.15269327261},
    {-608.23111709563, -611.860900820653},
    {933.958403926417, -494.314773848741},
    {-825.734172122654, 67.9558127768109},
    {537.794160442425, -256.804721526007},
    {442.509606603569, -472.510238387203},
    {384.293464397806, 932.899001715692},
    {-423.263018918359, 236.171413890373},
    {-720.604770861929, 167.571957536227},
    {294.52332147972, -426.803683746444},
    {90.6966411637704, 560.94343507297},
    {630.88487653779, -870.40540672867},
    {663.371717604423, -857.084446060314},
    {747.282628551656, 992.451862875284},
    {-5.66029397488819, -318.831337089184},
    {-65.4970559568869, 330.611875008452},
    {-647.14142852972, -63.0694946413923},
    {547.184707799658, 864.227535381996},
    {20.6774729265044, 798.08219808242},
    {-603.130121631464, 148.397670912045},
    {527.303202129636, -154.824767824155},
    {-322.438127153225, 655.743882918819},
    {-855.428150132822, 745.300273110262},
    {692.176932948897, -635.568032616548},
    {579.092781730445, -525.22774735116},
    {250.985156593349, 565.815386238578},
    {543.593814362137, 397.895458218894},
    {-788.873424496729, -757.956358769976},
    {-551.694943701401, -701.472358647437},
    {319.755346452898, -906.362605891525},
    {255.097466285864, -7.3249452440964},
    {-741.428336464695, -739.486925802857},
    {-258.946513281643, 240.165770988006},
    {-77.8447018573002, -551.549677415275},
    {641.482654413112, -670.315383720514},
    {-744.343346208358, 111.155304846288},
    {564.704638164186, -276.392263440141},
    {329.998789506632, -355.113945105368},
    {652.162707843485, -666.23645375463},
    {967.796194296767, -373.721613924012},
    {-828.242126840672, -446.053117694004},
    {-551.923762864197, -481.038486923627},
    {-931.578080977751, 667.84055508361},
    {781.03159803924, 233.975601626775},
    {632.813110868925, 577.841830143119},
    {-509.777196244622, -460.137689769736},
    {435.676339009185, 797.415966983637},
    {302.639480495174, -206.613352354078},
    {-116.701896826813, -615.137230958709},
    {865.681340395824, -389.320775056492},
    {-120.201669273614, 449.754530532494},
    {-846.274519929906, -558.981611095936},
    {-66.4023922034057, 567.451181374263},
    {-197.735586649899, -834.380808471257},
    {241.848395972941, 741.088300593496},
    {781.317050275758, -312.348689850855},
    {734.493122517137, -230.582207400764},
    {-263.963755893269, -409.948224079861},
    {-16.8788488901168, -589.312003953741},
    {436.644200553157, -3.78628958594072},
    {586.276204303799, 499.599840558964}};
Vector2 test_2D_vectors_2[n_test_vectors] = {
    {93.1181139999803, -430.521115932395},
    {608.699374276427, -896.335635249896},
    {145.643334119512, 42.7652339961037},
    {441.659413857497, -625.710267324692},
    {283.945375374343, 444.563508530915},
    {-274.590083527309, -174.72530778951},
    {758.62980923211, 828.005321777435},
    {-678.339218467836, 547.524058710369},
    {-606.397809900411, -235.707869741997},
    {-996.94935868341, 150.907981118864},
    {-992.86707710715, 409.37068341394},
    {832.290283801414, -521.64672177599},
    {907.519079694246, -836.669292584074},
    {502.407497508545, 807.771449520236},
    {850.611933566952, -951.799821540059},
    {-953.785192849502, -699.547321397824},
    {954.199444253683, -314.46661348477},
    {-558.548987001188, -229.914221398789},
    {543.086139390033, 231.576087650342},
    {441.674110799905, 687.580313065836},
    {-825.941324606928, 776.977027037087},
    {-487.089392011264, 603.85887782336},
    {-119.206145826258, -491.136277283381},
    {207.179448472721, -305.038554505785},
    {817.985943563376, 200.593595392007},
    {420.772906363426, 632.837666429961},
    {-837.225331256317, 743.366125113818},
    {920.887670795009, 435.478724236188},
    {-338.350213689001, -55.0758181066749},
    {586.475044323889, 975.722259621619},
    {631.800086660831, -209.874403118496},
    {311.64366116112, 642.089132103522},
    {880.524054137638, 274.468826203106},
    {248.231273385447, -99.6055790484975},
    {-122.468789367379, -506.154040648599},
    {348.665302067561, 638.143955736145},
    {943.572500168659, -192.525817331692},
    {46.2681602778589, -615.915224037033},
    {-532.716478270935, -887.00149502009},
    {-247.963165580265, 901.850747061593},
    {-866.284947719297, -339.572018271168},
    {423.87125747145, -448.810609837635},
    {70.382858443576, 774.89488613998},
    {-465.868638963917, -821.448865430092},
    {-918.528567047714, -473.931597435977},
    {662.59954640921, -757.333683426967},
    {-709.775732572232, -556.005738870159},
    {589.369163505657, -300.964816018237},
    {890.638008852276, -174.764317347522},
    {-844.444550129803, -833.599633933351},
    {-113.80918658791, 352.83765834507},
    {729.154754718767, -332.138297858695},
    {-535.021242803558, 275.091799410067},
    {984.460595157485, 989.035089196318},
    {-384.416354103317, -91.1333983976519},
    {648.366653815712, -186.559803667488},
    {-487.484905476699, 750.282554864814},
    {494.009572374028, -273.733510438657},
    {-697.421944539376, 269.090199031425},
    {-360.089565464576, -700.08629178128},
    {-689.97272552039, -365.878060227625},
    {411.755502000562, -463.929524164619},
    {4.69674889595308, 965.628252795312},
    {-768.437229552719, 88.2016316775612},
    {266.66112897063, 30.6849675143781},
    {254.228266237272, -906.185862158572},
    {-237.772962334526, 253.162305235067},
    {-3.42739183411368, 199.364040034838},
    {-490.378802708448, -995.438843900923},
    {706.16779232743, 416.1491861708},
    {-51.8604944258914, -610.829270367646},
    {737.642498432353, -455.967906450608},
    {-412.3515820422, -133.355704701822},
    {826.294779503355, -304.561021375283},
    {-524.336496619138, 229.412606719512},
    {267.727672584279, -457.03226598116},
    {337.496836788531, -306.868626200321},
    {438.103238788802, -998.68234387963},
    {-224.891782045701, 5.69323512226765},
    {-17.5714424259049, 880.218324978195},
    {662.078577059629, -273.991661170785},
    {-563.672151065977, -422.44891871765},
    {323.127535273179, -430.034500304253},
    {-481.837146193602, -136.92475315716},
    {761.806059848061, -569.029092899005},
    {-834.791918413466, -714.261225157433},
    {-300.774903188482, -977.805495101023},
    {278.75073530811, -722.19872814037},
    {-735.614397789123, 685.993145401594},
    {-857.8052931609, 813.046007651823},
    {-500.244464780862, 552.700731824729},
    {508.974956826456, -130.583149148699},
    {832.177595783631, 301.185680877076},
    {162.174173705666, 683.854162318418},
    {-237.281370714415, -291.96624278358},
    {601.117618125552, 917.787071685965},
    {474.375294864091, 589.211331687101},
    {793.3321990661, -476.479901600656},
    {-135.025571357256, 738.38239392667},
    {541.907339534588, 173.084874992223}};
Vector2 test_2D_vectors_par[n_test_vectors] = {
    {102.632861815351, -474.511481193545},
    {426.759880070003, -628.422016461615},
    {155.330171971547, 45.6095789840084},
    {-562.797456910011, 797.329653039938},
    {413.989078119002, 648.168461343754},
    {-699.081494143642, -444.834815828485},
    {-198.470856024377, -216.620706181119},
    {-328.824158719594, 265.411659952123},
    {593.035422919259, 230.513886982539},
    {976.344213526053, -147.788985325071},
    {168.566936802821, -69.5021153496442},
    {384.860103206865, -241.215132613673},
    {-656.155783739975, 604.929866148527},
    {529.905074158863, 851.982090203885},
    {-280.842449112102, 314.251167186019},
    {-356.678099439575, -261.603148103802},
    {786.612571060076, -259.236570337015},
    {1129.40821269259, 464.895677739476},
    {666.485213296413, 284.194397495279},
    {74.1261338518046, 115.39655387063},
    {-147.273791799223, 138.542956386274},
    {-579.8496249985, 718.856434980201},
    {210.479584645647, 867.188171637412},
    {-503.191698265091, 740.86918084592},
    {667.467451024557, 163.682147427734},
    {-308.964246156734, -464.678712890424},
    {235.727815571517, -209.300968689033},
    {-18.3620805439337, -8.6832473310139},
    {-808.01648699226, -131.526942393732},
    {-594.735370693341, -989.465017115597},
    {-122.461846861633, 40.6799675364264},
    {-220.474579977516, -454.250637350403},
    {-887.607148059259, -276.676702825463},
    {-553.039639791951, 221.913350429113},
    {-2.61543463991276, -10.8093892156721},
    {-173.369659231234, -317.309464091754},
    {680.146556029646, -138.77658747106},
    {34.3841896385275, -457.717482981052},
    {400.389451321258, 666.669901154447},
    {-146.961935234513, 534.505722939567},
    {-735.118813638351, -288.156662393225},
    {687.064492715883, -727.489369798496},
    {-0.634628864425061, -6.98707998672776},
    {20.6809324542586, 36.4659199605468},
    {156.930710481081, 80.9712674959167},
    {-295.694060509776, 337.970457762801},
    {-147.651663005325, -115.663537392536},
    {-639.432055613, 326.529725153973},
    {364.247440727404, -71.4739935771876},
    {326.3827421016, 322.191118761126},
    {313.747977091812, -972.69917189038},
    {872.705279934461, -397.527197530962},
    {187.440122799119, -96.3760623639122},
    {-162.23098217942, -162.984821047674},
    {12.1955281865715, 2.89118794513593},
    {-580.899394995107, 167.146901283503},
    {-232.409026450642, 357.698128044636},
    {-322.520068557526, 178.710202170587},
    {-574.815030417658, 221.784089463232},
    {47.3424325819257, 92.0431782782281},
    {19.738665418734, 10.4670001403919},
    {-746.94103484215, 841.586808651494},
    {-3.07491071656831, -632.186375835388},
    {631.064694510342, -72.4339394413675},
    {311.963297151327, 35.8979341148032},
    {-63.8208560384585, 227.48673193929},
    {8.50193656136875, -9.05220609486378},
    {11.8928803122001, -691.783368067729},
    {-296.857799133045, -602.603095280324},
    {186.138334285653, 109.692508168436},
    {-67.640631383396, -796.692703617943},
    {-294.771134379285, 182.210457383671},
    {-231.957201290656, -75.0156356514529},
    {782.275077435465, -288.33595768737},
    {-665.566393647441, 291.204831813326},
    {264.817751762043, -452.064801488695},
    {357.406119984586, -324.971119963351},
    {350.347207646672, -798.637260641052},
    {976.631227288849, -24.7238523085487},
    {8.5709036246859, -429.348157607615},
    {-301.259299787575, 124.671812156386},
    {-276.021345730485, -206.866560333554},
    {169.475746725402, -225.546912908966},
    {737.466457530438, 209.567513540078},
    {-106.602010611425, 79.6260998652696},
    {645.443078999699, 552.251350554371},
    {-31.9008153498561, -103.708096043324},
    {191.511082402865, -496.174691640754},
    {657.207163907706, -612.87491232431},
    {-287.873272087276, 272.852378559903},
    {-102.974509954461, 113.772547260579},
    {-198.896204243557, 51.0290189242107},
    {-441.841987953427, -159.913557702212},
    {179.266622245097, 755.929399767244},
    {157.906317526576, 194.298077852534},
    {114.80815490692, 175.289223141158},
    {-304.060136186249, -377.666543125079},
    {247.709695549812, -148.775874091641},
    {14.7989912954909, -80.9277421352476},
    {676.8033393925, 216.170575385966}};
Vector2 test_2D_vectors_perp[n_test_vectors] = {
    {-726.694914221932, -157.178027654165},
    {528.493766536511, 358.898845866012},
    {-204.476565514215, 696.375208738726},
    {71.9993116929693, 50.8209237742281},
    {-535.709087754377, 342.160602692125},
    {173.736380946684, -273.035932569961},
    {-67.8980080579103, 62.2090843445791},
    {34.4282431451963, 42.6538837459322},
    {-0.198553773090566, 0.51081269913243},
    {-115.974801257239, -766.168912204654},
    {191.817579696726, 465.22471541686},
    {546.357207182873, 871.716002498823},
    {-119.718423904919, -129.856270389797},
    {237.456850253268, -147.690416605912},
    {581.116623698097, 519.336864459512},
    {175.967019356959, -239.919062454794},
    {212.267490607825, 644.092291154845},
    {-182.17909172134, 442.582222685799},
    {147.529386405831, -345.98203002159},
    {765.272282801844, -491.580326840966},
    {482.59019707319, 513.00253771071},
    {-313.054177152838, -252.518219762875},
    {-203.851145673854, 49.4777326009381},
    {245.227146471058, 166.556076948134},
    {-168.948718282299, 688.943614914019},
    {378.259562234524, -251.504270058793},
    {-524.449446274469, -590.667702697786},
    {414.277640989415, -876.054674180315},
    {-24.6448949043295, 151.402298574523},
    {-69.0913517392406, 41.528573499371},
    {-205.682628148745, -619.181283463302},
    {727.127084078212, -352.917586175573},
    {201.693737809455, -647.054130579965},
    {-333.085249504337, -830.095828167549},
    {-810.237273047804, 196.044622706827},
    {853.58459658826, -466.376479060029},
    {159.134325795524, 779.920198416007},
    {741.928137965573, 55.7343748982905},
    {487.600945301127, -292.843991628798},
    {264.21958403902, 72.6469703330429},
    {126.887696542721, -323.704238427428},
    {246.893911210534, 233.174595949754},
    {-825.099543258228, 74.9428927635387},
    {517.113227988166, -293.270641486554},
    {285.578896122489, -553.48150588312},
    {679.987524907581, 594.928543952891},
    {-275.611355913034, 351.834951282909},
    {-81.1727152489299, -158.957767617746},
    {-69.7241192476845, -355.329690169256},
    {-235.686100937829, 238.752316311844},
    {317.136899445978, 102.293765161711},
    {-209.333562330038, -459.557248529352},
    {559.842505752537, 1088.8279252392},
    {156.570688204532, -155.84651604151},
    {-77.6925841434584, 327.720687063317},
    {-66.2420335346128, -230.216395924895},
    {779.5937342503, 506.52940733736},
    {343.19754148403, 619.371995911833},
    {-28.315091213806, -73.3864185511873},
    {479.960769547711, -246.867946102383},
    {-342.176792571959, 645.276882778427},
    {-108.487115290672, -96.2865355412318},
    {695.251843665465, -3.38165678116036},
    {-51.9719127798975, -452.793807909793},
    {-60.9781405579774, 529.917452123775},
    {607.414670400596, 170.408726279604},
    {-797.375361058098, -748.904152675113},
    {-563.587824013601, -9.68899057970816},
    {616.613145585942, -303.759510611202},
    {68.9591320002108, -117.017453412532},
    {-673.787705081299, 57.2057778150863},
    {35.824621097642, 57.955313604335},
    {154.112499433356, -476.534041763822},
    {-140.792423022354, -381.979426033144},
    {-78.7769525609168, -180.049526967037},
    {299.886886402143, 175.672538048554},
    {-27.4073304779542, -30.1428251420177},
    {301.815500196813, 132.400806886422},
    {-8.83503299208235, -348.997761615463},
    {-836.813030465357, -16.7049600863891},
    {-250.664463076622, -605.710299080013},
    {-655.556735247266, 874.707115417164},
    {611.555851313838, 459.52251453574},
    {-104.653346661513, 368.27431660304},
    {-403.175185633197, -539.763789635005},
    {-209.766739990514, 245.164616429267},
    {334.54029584503, -102.905256310754},
    {-308.212979229678, -118.962539317955},
    {208.474176488118, 223.554137267818},
    {167.671602813662, 176.902151972591},
    {-743.300009975445, -672.754158356515},
    {132.493812040152, 516.422162450052},
    {244.106401303528, -674.467250769045},
    {62.5817737278445, -14.8410991737483},
    {623.410732749182, -506.646767703389},
    {619.684967610217, -405.871430541921},
    {40.0963802929803, -32.281680954782},
    {-264.588544439929, -440.5361298621},
    {421.845209257666, 77.1414525493068},
    {-90.5271350887009, 283.429265172999}};
double test_2D_vectors_dots[n_test_vectors] = {
    213844.190932845,    823045.519295036,  24573.3084530083,
    -747462.145288889,   405702.329481442,  269684.745980421,
    -329929.00516906,    368373.592125236,  -413949.318898205,
    -995668.274936069,   -195816.690299721, 446144.407692566,
    -1.10160013616923e6, 954435.090194507,  -537992.143513331,
    523198.071384422,    832106.524513773,  -737715.940889402,
    427771.508153794,    112083.992893117,  229284.205048471,
    716526.441441724,    -450998.0303829,   -330244.642525738,
    578812.583178177,    -424069.976123007, -352944.548555127,
    -20690.7830529907,   280636.504993471,  -1.31424049520494e6,
    -85909.0893654127,   -360378.902790946, -857498.574352936,
    -159385.541786269,   5791.52514247544,  -262937.101213236,
    668485.662283365,    283506.059273161,  -804631.257453583,
    518485.532222297,    734672.302467321,  617731.838201125,
    -5458.91954429538,   -39589.4864334259, -182520.182774722,
    -451883.162036675,   169109.157843386,  -475135.494390788,
    336903.719034611,    -544190.526479968, -378912.300150412,
    768371.211072855,    -126796.711866981, -320907.716291867,
    -4951.64426469261,   -407818.690023628, 381670.557670914,
    -208246.972141502,   460568.541040008,  -81485.6833410296,
    -17448.7864848679,   -697994.048447157, -610471.467622439,
    -491322.397165623,   84290.0119579189,  -222370.325892477,
    -4313.2080042155,    -137957.488647756, 745427.30061049,
    177093.344593333,    490151.099444998,  -300517.836815915,
    105651.681874238,    734205.706392993,  415786.810668145,
    277507.240932954,    220346.876081593,  951073.177737232,
    -219777.09581079,    -378070.719261428, -233616.365443238,
    242976.100419778,    151755.234315935,  -384033.713371207,
    -126519.625053191,   -933262.392386665, 111001.310846113,
    411720.586277843,    -903879.040943744, 468780.753622446,
    114394.598750726,    -107896.716956903, -415854.677009359,
    546017.882785298,    -94196.7072378472, 229891.387423868,
    -366764.023568159,   267404.791348263,  -61753.8622280891,
    404180.554055991};
Vector2 test_2D_unit_vectors[n_test_vectors] = {
    {0.837579537839011, 0.546315401387686},
    {0.566801149046583, -0.823854633682104},
    {0.443186793716114, -0.896429286600806},
    {0.985316726914574, 0.170736486031404},
    {0.750469094984102, -0.660905543533827},
    {0.601864988718647, 0.798597855841539},
    {0.959444331617004, -0.281898163399482},
    {0.945428113944436, 0.32583075570527},
    {0.964651831309407, -0.263527691811331},
    {0.427982458229142, -0.903787040982631},
    {0.825363052287148, -0.564602366200536},
    {0.899928189409094, -0.436038133546678},
    {0.993102395317923, -0.117250298139513},
    {0.891980601857183, 0.452073673100411},
    {0.993777801156305, 0.111380796948756},
    {0.820866233059481, -0.571120501665575},
    {0.912582891120166, -0.40889175442256},
    {0.969267851924649, -0.246007786920206},
    {0.568495388879313, -0.822686448668603},
    {0.228406200039923, -0.973565923696656},
    {0.180134457701833, -0.983641996433899},
    {0.958671083201719, -0.2845167029066},
    {0.674017913557811, 0.738715000662079},
    {0.695971359817567, 0.718069541418996},
    {0.918131857531417, 0.396275020895729},
    {0.71673055870947, 0.697350203421503},
    {0.93115037295469, -0.364635410987938},
    {0.556205549813111, -0.831044755929002},
    {0.277775343343923, 0.96064606314197},
    {0.62373206251763, 0.781638224620255},
    {0.88940671874252, 0.457116712290925},
    {0.999127889402885, -0.0417547676000862},
    {0.8736706497434, 0.486517826782273},
    {0.952878511802119, -0.303352174453687},
    {0.943980512938603, -0.330001198773841},
    {0.964595387571009, -0.263734219009848},
    {0.850429910145746, -0.526088365134127},
    {0.995765526466305, -0.0919294093382628},
    {0.940201668770551, 0.34061829375574},
    {0.983602216520436, 0.180351544645688},
    {0.460028989521536, -0.887903896150814},
    {0.908254638738604, 0.418417866743054},
    {0.935974450538335, -0.352067930859175},
    {0.0106434917098523, -0.999943356437865},
    {0.944189133533958, 0.329403825291075},
    {0.984758923223826, -0.173924877838254},
    {0.980395238856322, 0.197041050615997},
    {0.896468494296951, 0.443107479888297},
    {0.758074305332329, -0.652168189652721},
    {0.793627018676609, -0.608404598294981},
    {0.891573094672948, -0.452876823049383},
    {0.818992022608801, -0.573804903171057},
    {0.811373894561413, 0.584527504249582},
    {0.997860833984577, -0.0653739703521499},
    {0.671284223340055, -0.741200034737411},
    {0.999879970073819, 0.0154934000522737},
    {0.941876796192017, 0.335958480760736},
    {0.968106628152906, -0.250538533017201},
    {0.750476504932986, -0.660897129320116},
    {0.980275422078795, 0.197636274171114},
    {0.887791276590295, 0.460246291902687},
    {0.440342755413766, -0.897829748757866},
    {0.999934117902939, 0.0114786695079041},
    {0.991037535872448, -0.133583690965124},
    {0.913533925598441, -0.406762543482931},
    {0.909624601796285, -0.415431202254898},
    {0.951053307249748, 0.309026870626677},
    {0.539586224786673, -0.841930345112032},
    {0.530360951552902, 0.847771939302016},
    {0.816471723232211, 0.577385421674486},
    {0.895380803605855, 0.44530126491414},
    {0.708048758592216, 0.706163547243853},
    {0.996985169937687, 0.0775923380516548},
    {0.999999995621825, 0.0000935753736750478},
    {0.636132771474851, -0.771579611612259},
    {0.646768997925245, -0.762685953274855},
    {0.903356172237091, -0.428891158781749},
    {0.999681161282385, 0.0252502629115337},
    {0.771181760686304, -0.636615026515062},
    {0.489058824858116, -0.87225080442978},
    {0.982766393055924, 0.184851877674666},
    {0.96656762140353, -0.25641184304615},
    {0.997433972871028, -0.0715923862063406},
    {0.933061370939244, 0.359717219578072},
    {0.996154478095304, -0.0876142441083156},
    {0.715700267018375, 0.698407565673388},
    {0.649947339364229, 0.759979247120183},
    {0.904212049244557, 0.427083797399244},
    {0.418885047095614, -0.908039270802593},
    {0.770039447772223, -0.637996276536665},
    {0.687504469206222, 0.726180146259502},
    {0.951831558827842, -0.306621401111144},
    {0.92002192935565, -0.391866877274295},
    {0.694691561603252, 0.719307746543324},
    {0.682052535685107, 0.73130317828211},
    {0.755460927513899, 0.655193701892685},
    {0.366802396644996, 0.930298877681516},
    {0.995699365225257, -0.0926432625182252},
    {0.384047483915534, -0.923313343398734},
    {0.131125757018221, -0.991365742723844}};

Vector test_unit_vectors[n_test_vectors] = {
    {-0.940191081449973, 0.0526929562443838, 0.336547444982365},
    {0.285154150776666, -0.876601943331571, -0.387628873075977},
    {-0.84479470669357, -0.384714882553555, 0.371909078518298},
    {-0.761223803795772, 0.30124584901417, 0.574272808851706},
    {0.777061924065319, 0.607393701297124, 0.165068645697735},
    {-0.498743044341681, -0.503561559158522, 0.705465188268444},
    {-0.836004753337821, 0.220852277574427, 0.502314964824618},
    {-0.023692565959066, -0.851872732793013, -0.52321268088803},
    {0.563183785165542, 0.808288330403872, 0.171738170071575},
    {0.283937289874856, 0.909158191825664, 0.304648974485079},
    {0.443313732852368, -0.616956401598036, 0.650259742558081},
    {-0.932601793208028, -0.356099270216329, -0.0587129036632308},
    {0.0209946957504972, 0.946770470580609, 0.321224063212761},
    {-0.721410600451064, -0.442178682848418, 0.53295849556159},
    {0.327678074171689, 0.71879885232206, -0.613151930281249},
    {0.0521821199768833, 0.22488112678303, -0.972987926529159},
    {-0.730634227221078, -0.64781600171384, 0.215657260338139},
    {0.513870729669331, 0.781696539782829, 0.353394104196239},
    {-0.000494134438625825, 0.140928892086025, -0.990019597385103},
    {-0.419549111262533, -0.708463719541513, 0.567501278705364},
    {-0.572125012655014, -0.226315133598435, 0.788323810498467},
    {-0.755635560035837, -0.386088222233196, -0.529103756424141},
    {-0.389923427697203, 0.184350960351157, -0.902205322501741},
    {-0.078857873911996, 0.958835590110084, 0.272792497807999},
    {0.43394182728908, -0.654822126518626, -0.618791138552128},
    {0.579808681350739, 0.802036505199273, 0.143385275945805},
    {0.545628549977943, -0.712667105116986, 0.440902577371853},
    {0.0720168706622575, 0.822980979768596, -0.563485472110096},
    {-0.582123673230562, 0.297618280632705, 0.756673898121106},
    {0.245903991949659, -0.969169439723819, 0.0155506864363284},
    {-0.0731738990444161, -0.92907716026367, -0.362575800040534},
    {0.60078797575626, -0.150587803808877, 0.785096886715719},
    {0.926906566433105, 0.368675948686178, -0.0701588338239612},
    {-0.454914617914336, 0.828667180859799, 0.326134011985451},
    {-0.0379495995705982, 0.999010481295923, 0.0231923727402013},
    {0.950330256830006, 0.300549787044449, -0.0808840433024387},
    {0.300315853297785, -0.915413652967132, 0.268007895815394},
    {0.278432026632119, -0.117589747573276, -0.953230432692528},
    {-0.424205837196732, -0.905561119318287, 0.00290978817639287},
    {-0.600248885923854, 0.0191462033122512, -0.79958407803426},
    {-0.435972967202716, -0.838954102367944, -0.325704752787639},
    {0.72991886075342, 0.530128873904871, 0.431487930038353},
    {0.812997050731825, 0.199124819874075, -0.547160946716294},
    {0.827461182937613, -0.466402498082972, 0.312692661432035},
    {-0.51856679173243, 0.297657883207708, -0.80155365826415},
    {-0.48456053520024, 0.362176808785229, 0.796259409303655},
    {-0.830879181162632, -0.511375510191297, -0.219396613207917},
    {-0.170691949475731, 0.364115129018401, 0.915578741127211},
    {0.564929078212686, 0.335248010784533, -0.753965455345786},
    {0.615246619543672, -0.324725943844315, -0.718348563396979},
    {-0.524898165875001, -0.53006096250499, -0.665970938922521},
    {-0.750440996037387, -0.00122236678453737, 0.660936318631272},
    {-0.694033094798214, -0.397995241820529, -0.599931538438372},
    {0.196680774292255, 0.978060932471966, 0.0686548279137518},
    {0.151540233864963, 0.579756925395258, -0.800573210253996},
    {0.960081385229715, 0.0627295793785936, 0.272596283185184},
    {-0.0823627261574278, 0.220433142870447, -0.971918520692025},
    {0.243739695212647, -0.633385459434923, 0.734447970079609},
    {0.160692983744955, -0.718068731900344, 0.677166937499296},
    {-0.319271616382626, 0.804860198134361, -0.500265625874437},
    {-0.410754695383332, -0.907390375750823, -0.0890128429796904},
    {0.293131132842691, -0.861871671083773, -0.413825279002655},
    {0.541856119552942, -0.770399315419703, 0.335971487635306},
    {-0.0988090681976125, -0.207448493517533, 0.973242976126325},
    {0.581102797916053, 0.573722474480777, -0.577201923532648},
    {-0.0309757017243094, -0.984957569497875, -0.169997329954108},
    {0.54473861066991, -0.0576431703410422, 0.836622442299064},
    {0.856657096474634, 0.070195418565432, -0.511088272485365},
    {-0.990212818777598, 0.10772230505777, -0.0887382584997254},
    {-0.43161402485055, 0.274265276636298, -0.859353182099143},
    {0.833258243067225, -0.481298370289039, 0.27208928519447},
    {-0.601615919420646, 0.260277951125553, 0.755191150409971},
    {-0.645985758038829, -0.763082716851565, 0.0201783957101349},
    {-0.955018723825842, -0.295961853140295, -0.0185961992843602},
    {-0.964463423199809, -0.0964651364243642, 0.245977199683904},
    {0.185601370015039, -0.32096600720726, -0.928726522538239},
    {0.00232906002723656, 0.96227932760902, -0.272053434339143},
    {-0.22717638144596, 0.931772313414168, -0.283180591969027},
    {0.0459238510463156, -0.948943209159229, -0.312086183122003},
    {-0.0651781520258082, 0.18367855803644, 0.980823121575014},
    {0.289606810896222, -0.887316196396032, -0.358884191203512},
    {-0.829603863961748, -0.545464485616881, -0.119272477254604},
    {-0.444188210001824, -0.55078786860278, -0.706632548001705},
    {0.441901547808439, 0.448310201366617, 0.777007712571201},
    {0.158957750224097, 0.984273887043946, -0.0770541947404229},
    {-0.809243011111887, -0.233922531343055, -0.538893308825248},
    {0.433453007582328, 0.783746335162492, -0.444814761824736},
    {-0.52691277764331, -0.349208182291332, -0.774865517478353},
    {-0.00416063360936794, -0.0713466251989094, 0.997442904732243},
    {-0.702506869429872, 0.69108579931978, -0.16995445384685},
    {-0.601272045664885, 0.386630229084162, -0.699277479302954},
    {0.0498874487901271, -0.00021308643237741, 0.998754823291174},
    {-0.955485767728463, -0.290714600577509, -0.0503186712802513},
    {0.060250873248179, 0.0164196116082409, 0.998048209570793},
    {-0.85070116166329, 0.421810928510214, -0.313660762822001},
    {-0.135390158730839, -0.948197389774372, 0.287386873297834},
    {0.322975576828363, -0.597826269447782, -0.733682852689454},
    {0.721179520737365, 0.0688953172724899, 0.689313813967882},
    {-0.824776429218411, 0.372264385681751, 0.425632551573205},
    {-0.480333791854488, 0.869036505691313, 0.118553786015136}};

double test_vectors_dots[n_test_vectors] = {
    872073.213594283,    -417268.600740014,   -249286.70880553,
    -900794.283262525,   956105.148565154,    1.82963250702583e6,
    -836596.527470191,   880476.437671783,    525070.061754772,
    -51406.3529127657,   -153188.395224371,   756649.802477262,
    49187.6432264308,    -434125.912651637,   351674.430705344,
    400373.392442367,    428431.396185661,    -781052.50095864,
    281852.217373467,    -352453.962101088,   -642941.359374965,
    207188.356951848,    -681254.176566828,   567605.18871948,
    -434315.048675848,   -1.00000695973681e6, 54886.983719829,
    525174.255876388,    183890.974868801,    420053.242567062,
    -41842.6130917062,   -178458.075757046,   -219092.467444318,
    -284410.717811296,   326893.4981799,      186181.479290782,
    -1.25844302925506e6, 197237.560876716,    -39401.5415337377,
    31655.0333493078,    -237935.890219069,   506155.526254458,
    -531440.762723106,   820968.464881545,    -540108.850063068,
    -58910.5519391676,   -534709.636170739,   618736.731893291,
    -30016.0459335557,   -558498.918705161,   -419214.910680556,
    -149038.546539245,   525405.184814948,    -207460.841624904,
    324112.892110184,    203774.834154969,    -83486.1200708506,
    -284685.600680733,   295163.48773739,     31408.437792288,
    -944308.026308756,   -813072.66592001,    -523934.489405215,
    514651.917266344,    634641.875917128,    565050.284496592,
    -62899.5914294491,   -282515.635332808,   -894877.334131346,
    583074.081738119,    702442.4007937,      337348.510791781,
    258612.948347045,    -1.12876680496077e6, -443818.047932314,
    386536.240962004,    739818.180991229,    145015.252280881,
    -288908.15294173,    1.19116251179855e6,  -77366.8414122704,
    714606.254594206,    -165460.662260136,   343248.647838167,
    221797.734748679,    -51971.5645824749,   -12429.5964039363,
    -1.04649426548902e6, -532707.457457638,   101230.127743256,
    61025.3514233649,    654342.520440778,    323363.029606064,
    1.13028462841734e6,  56425.6890927089,    -11433.5336198783,
    -688931.982781608,   197145.107764405,    -236255.30773955,
    373952.596888678};

Vector test_vectors_perp[n_test_vectors] = {
    {-408.377697688621, 441.17302504049, 1042.87330481857},
    {-405.365911346278, 988.627145631829, 792.547414278057},
    {937.727848909246, 622.402895923916, 371.355070973412},
    {128.741499253745, -673.607903874247, -577.632159588393},
    {20.6422259929184, 140.819356460295, 491.673427427005},
    {188.427830092756, 169.066444086588, 360.641424505715},
    {958.668783065098, -475.209797384555, 545.180066416736},
    {-603.625143617549, 162.490309711147, 341.234806943597},
    {-50.7477795038651, 613.486244106002, 863.713708173206},
    {666.447888589294, 717.859974079491, -162.47352764292},
    {109.572396768264, -124.276166214215, 82.7228656271971},
    {197.007963417961, -318.095572497207, 426.801094866039},
    {-613.097366676253, -129.546864950581, -208.051562262983},
    {432.014574143453, -338.485389837814, -843.10492608124},
    {-172.982210718092, -809.571718986696, 648.57564357874},
    {-555.496490776957, -534.376072790188, -624.600525775925},
    {-628.610326883537, 310.828695986221, -264.424326442599},
    {684.457487411529, -100.435851686672, 388.628822365331},
    {411.024023481066, 305.320301541885, -700.187606157989},
    {64.3833245453543, -606.227742457069, -204.385727919895},
    {104.330098738787, -404.013119703574, 487.661237568957},
    {-632.249248949869, 323.05665676358, -696.40922415236},
    {-226.455751000243, 167.149637952707, -151.605276418416},
    {142.475424220036, 184.743213058667, -74.8673495331142},
    {-351.706262982528, -55.7892493166133, 387.322799856731},
    {-203.508199547425, 249.176824324507, 178.691692365619},
    {591.778313429372, -129.527403727817, 44.5338976617646},
    {609.721912984159, -362.500994645967, -485.520574867091},
    {-205.971592862358, -912.579777567569, 751.667111956061},
    {-683.276883007809, 538.696683067137, -251.908301544944},
    {275.940676914606, 326.688935548135, 366.150565007911},
    {902.129452633543, 487.724878051576, 363.249245706984},
    {-142.759137951756, -375.935038944178, -311.45948430378},
    {371.041644624913, -169.430477175702, 26.7864032192385},
    {814.338205845452, -621.097408665374, 471.115690833946},
    {426.35665590009, 402.617805891505, 466.534942837742},
    {191.296410405189, 163.925567261102, -241.445531165462},
    {-212.125694235698, 558.001214710534, -379.848311203323},
    {638.989659705758, 948.858556725613, -478.530598514333},
    {-86.2714299936597, 537.533587394986, -853.010988612088},
    {-347.499549438443, -162.095193464804, 95.7528626811276},
    {305.956029356022, 128.828399511602, 542.698995728158},
    {-194.574005627682, 840.357416867734, 106.977089334425},
    {345.642545955331, 791.057788070593, 487.506518679426},
    {269.680878979321, 321.491565433663, 994.034649417831},
    {-293.779531714365, -789.147292332311, 943.83832734282},
    {-176.260742213012, -363.005368180402, 554.908303560903},
    {486.065981241204, -146.924792646647, 523.559629504705},
    {250.659658572808, -705.034118376533, -301.041960798922},
    {-768.017181749818, 97.1030267551249, 145.194481207599},
    {-533.652732080419, -660.460546174441, -306.070287334517},
    {366.001600854584, -366.549462739307, -110.681850813767},
    {-31.2300412266863, -149.62734682331, -395.094815989902},
    {169.885798013999, -379.060888566106, -1004.08553869076},
    {-422.050669392612, -660.737196684259, -98.7458218856877},
    {-267.175594748877, 22.4136010487385, 288.343936616289},
    {360.701044046773, -453.268554409309, 220.035893279262},
    {574.98492632679, -790.415308866877, -232.237640689443},
    {-0.614245533659229, 805.932589326686, -6.4376715179074},
    {-877.55266284643, -866.896766346259, -816.913638506201},
    {-72.2003650488606, -475.046381452584, -799.044079739731},
    {957.184264249099, -83.6802495648801, 160.277468698834},
    {-258.709154812169, -361.258618894944, -210.178023951695},
    {455.707673140088, 595.345313548196, 68.6941351518952},
    {-214.751525600962, 102.468439805783, 354.304345464676},
    {989.648782278563, -463.401626522377, 125.595860450194},
    {-280.41403391793, 374.198200421694, -508.201438569916},
    {115.145173772755, 79.5733442311659, 602.3278334966},
    {-134.61653977559, -129.067947508342, 13.8590639958292},
    {327.721756077983, -347.735308572441, 380.71152246526},
    {-72.1181771662457, 206.597259721751, 187.30639858793},
    {1105.72938034451, -330.226169957177, 94.8121647561526},
    {558.355423951927, 347.573201152264, 163.646206303919},
    {-127.87730075875, 194.169524231086, 986.625176597396},
    {383.65900387544, 254.623921560605, -1011.84339605979},
    {364.900854111963, 1087.18025392693, 235.39145705154},
    {122.63191483612, -243.229832370685, 743.312457067795},
    {43.6997145993886, -659.621167746111, 73.5599990981897},
    {-90.6881484816376, -284.04779122735, -686.482441170206},
    {-157.359769391139, 502.619209555219, 657.741059200661},
    {809.861455901173, 732.679131545009, -983.51183520927},
    {-277.606190326059, -330.641314380112, -145.992350784963},
    {187.38450184376, 993.802329967896, -6.36099153642996},
    {-191.425960470862, 701.532551189151, 519.297763147776},
    {406.954807819086, -444.572292793205, -219.950906648917},
    {524.712936056976, 636.6981440502, -268.198408688161},
    {12.6498014474297, 172.846911567105, -791.78816998571},
    {74.291795973226, 716.772212165161, -42.2904659287665},
    {-70.0203037998306, 316.17270724844, 175.804078465356},
    {764.584100714417, -820.479973067382, 1017.76932675583},
    {236.955720016637, -46.1746837560523, -708.982622435074},
    {24.2173968662735, 193.529525449192, 436.349833312144},
    {282.876242877252, 108.085785795468, 39.0285335393357},
    {296.317474441846, -284.40956630172, -337.565436769693},
    {-16.630109999855, 197.302139275305, 326.165073190749},
    {-763.110648488195, 976.876331853261, -63.9210595767069},
    {-551.249069574625, 455.995470283257, 398.150405793844},
    {-31.018397143801, -18.6919905376747, -98.5157679347456},
    {149.038786819234, -653.81146410717, -820.988805429921},
    {-45.9175017195183, -828.356916337882, 400.431044101311}};

Vector test_vectors_par[n_test_vectors] = {
    {-337.556830080855, 347.994939078823, -279.397947702155},
    {-215.472758277993, -197.920586024837, 136.67870343446},
    {-99.1281001554395, 182.23218086168, -55.1134441761401},
    {-891.558146520862, 174.610160306156, -402.330986383383},
    {936.375661300756, -715.695030394048, 165.66857398283},
    {613.837580621119, 657.769660405623, -629.075989175644},
    {-347.512472109573, -314.455679398457, 336.983962475677},
    {-370.188165344673, 318.480912940757, -806.497288824398},
    {288.632964952744, 343.479165861746, -227.010709075023},
    {12.0823928330085, 1.6766683245441, 56.9686545675581},
    {162.782673721193, -95.295591567799, -358.781798286511},
    {-78.4678342868004, 519.835274972168, 423.654226285264},
    {-15.5832609355262, -6.70913804121251, 50.0991385509741},
    {-432.257394605379, -219.113423138983, -133.524070751375},
    {-474.920816580623, 239.013039384116, 171.676728055789},
    {-330.371115076076, 83.6919485835169, 222.217264543907},
    {-32.0020913439104, 383.235724390766, 526.568441873651},
    {131.838565600842, 776.165158514938, -31.6062111280493},
    {-242.527676747365, -72.7784894055928, -174.10398407679},
    {-377.085612146422, -9.75076543917806, -89.8636172677568},
    {-118.569838121511, 351.928627911577, 316.929404856295},
    {90.0494354070956, 188.509960750864, 5.69451070545296},
    {5.24577220245978, -537.14399084347, -600.0540418565},
    {-211.721879105174, -38.3287755277796, -497.495182052693},
    {-195.905644997084, 587.483467614011, -93.2709891293835},
    {-654.336935656091, 64.3150378257196, -834.894709314297},
    {-11.5838015078863, -71.4415963619148, -53.8601402848086},
    {244.68662423915, -92.9529778563236, 376.680933843054},
    {-120.887743126646, 178.53904618858, 183.634058019873},
    {-116.807712019822, 155.992057338406, 650.41291711446},
    {-1.1480735229705, 30.3664742236851, -26.2285296629321},
    {-144.83549482068, -118.941304845721, 519.39818539322},
    {-178.877022828204, -126.208437818694, 234.32422913024},
    {-338.348543310098, -609.583386936684, 830.996071289328},
    {-158.63287543198, 64.7994065430085, 359.630464343882},
    {-185.664607784456, 78.9763213432479, 101.518801119799},
    {762.937890953595, 126.232514349986, 690.176436974295},
    {91.2194307076713, 139.782247737592, 154.400262515587},
    {17.0306588075589, -20.1486370116606, -17.2105854556319},
    {13.9597170974557, -15.3253507077363, -11.0692776838472},
    {-1.69274211039583, 128.698288700305, 211.723663576514},
    {296.427127604079, -486.374093308421, -51.6582325948642},
    {-621.912699168388, -150.34244954326, 49.8550441875622},
    {-786.867349749382, 195.790671048476, 240.187760461586},
    {634.779024161494, 312.317750913365, -273.225171784936},
    {-18.2102183161708, 37.1774976009916, 25.4161453961165},
    {-564.365056661908, 396.754361148602, 80.2809381184957},
    {292.250395291636, 264.357573977641, -197.135698777928},
    {-30.173761940752, -7.48105469377835, -7.60341202700412},
    {-33.1514428524145, 341.611537083467, -403.81970063927},
    {119.365351586425, 106.891915819217, -438.780060171762},
    {54.0698116037803, 93.4050557273382, -130.535722730293},
    {-62.7951823487097, -619.174235782353, 239.452633712788},
    {229.519324814653, -100.800714274542, 76.8875548964105},
    {152.726278923802, -122.836458245018, 169.165525171695},
    {221.730007770667, -79.2470214157614, 211.612071813088},
    {-42.2699229939383, -6.47904166624477, 55.9456883248148},
    {-231.825406359077, -171.960543585465, 11.2993396211691},
    {-38.9550492691053, -2.67809787767646, -331.554411569742},
    {52.0213812813663, -105.037210101617, 55.5810479523133},
    {747.677239726924, -382.899301877474, 160.081729423827},
    {-1.80862245308594, -706.530558614885, -358.075710728166},
    {233.362324324246, -235.136272618468, 116.910583444114},
    {392.260114511648, -318.372259940119, 157.007419126559},
    {650.557132543314, -523.425416630949, 545.696729312073},
    {-153.002069362163, -175.182645761935, 559.241271121706},
    {72.6306907607094, 34.4003285585938, -14.7463650693395},
    {137.172631936416, -199.950141570442, 0.192478076925797},
    {-471.182754655085, 432.851008885649, -545.62131627997},
    {203.109102965708, 741.941543265842, 502.837419659091},
    {253.505103191105, 329.235303756665, -265.53703445402},
    {-108.855754092289, -271.32851335924, 324.486101880809},
    {231.64603464732, -123.486595727059, -528.091609685965},
    {-500.941507194066, 778.204982219084, -218.079513916091},
    {-180.024577637769, 530.982130307993, 65.3586339798718},
    {449.134937485257, -140.885334728695, -45.5495218444687},
    {-840.275383846644, 169.128788082728, 193.971922189968},
    {426.820579437138, 20.4595828300042, -70.0973851527666},
    {181.346114107287, 164.220890454323, -91.9069749579005},
    {1094.47688419908, 346.525700272901, -2.95533239569392},
    {-60.8177337482218, 75.7930541939831, 6.38329962212474},
    {-414.379274695757, 385.470899728847, -85.061669851163},
    {50.5849633413943, -10.7731585345221, -192.98121828637},
    {254.28386911783, -92.746634770236, 219.02909131217},
    {-72.5062990672491, -160.318457483553, 189.889452397731},
    {78.3973762635371, -41.1286706607777, 55.7407080622059},
    {20.657272360243, -24.5488976181432, -5.02898488865074},
    {-1070.5366005096, 126.48808537174, 263.204904771135},
    {-455.775259519712, 28.0988529029689, -232.062946918288},
    {79.1953869136452, 32.1064234156617, -33.6116007482456},
    {-141.311672357206, 23.233273581278, -48.7422357570619},
    {699.716116483646, 155.918466587939, -107.98704648624},
    {98.2638727061121, -170.734852410573, -239.376261134367},
    {685.509647626551, 602.9914246589, 93.7061513821605},
    {109.054820079772, 46.1050804377982, -22.3292987095441},
    {7.10693539629225, 5.06640546964724, -7.41737529536779},
    {-10.6878763453803, 359.268115526314, -426.261816435765},
    {-231.436807010674, 292.67846636608, 17.3378912377333},
    {315.501308196263, 21.3372949862042, 40.2823569818559},
    {-628.213008298047, -164.364907891286, -412.052917910106}};

Vector test_vectors_1[n_test_vectors] = {
    {-745.934527769476, 789.167964119313, 763.475357116419},
    {-620.838669624271, 790.706559606992, 929.226117712517},
    {838.599748753806, 804.635076785596, 316.241626797272},
    {-762.816647267117, -498.997743568092, -979.963145971775},
    {957.017887293674, -574.875673933752, 657.342001409836},
    {802.265410713875, 826.836104492211, -268.434564669929},
    {611.156310955525, -789.665476783012, 882.164028892413},
    {-973.813308962222, 480.971222651904, -465.262481880801},
    {237.885185448878, 956.965409967749, 636.702999098184},
    {678.530281422302, 719.536642404035, -105.504873075362},
    {272.355070489457, -219.571757782014, -276.058932659314},
    {118.54012913116, 201.739702474961, 850.455321151303},
    {-628.680627611779, -136.256002991793, -157.952423712009},
    {-0.242820461926385, -557.598812976797, -976.628996832616},
    {-647.903027298716, -570.55867960258, 820.252371634529},
    {-885.867605853032, -450.684124206671, -402.383261232018},
    {-660.612418227448, 694.064420376988, 262.144115431051},
    {816.296053012371, 675.729306828267, 357.022611237282},
    {168.496346733701, 232.541812136292, -874.291590234779},
    {-312.702287601068, -615.978507896247, -294.249345187652},
    {-14.2397393827241, -52.0844917919967, 804.590642425252},
    {-542.199813542773, 511.566617514444, -690.714713446906},
    {-221.209978797783, -369.994352890763, -751.659318274916},
    {-69.2464548851385, 146.414437530887, -572.362531585807},
    {-547.611907979612, 531.694218297398, 294.051810727347},
    {-857.845135203516, 313.491862150226, -656.203016948678},
    {580.194511921486, -200.969000089732, -9.32624262304398},
    {854.408537223309, -455.453972502291, -108.839641024037},
    {-326.859335989004, -734.040731378988, 935.301169975934},
    {-800.084595027632, 694.688740405543, 398.504615569516},
    {274.792603391636, 357.05540977182, 339.922035344979},
    {757.293957812863, 368.783573205855, 882.647431100204},
    {-321.63616077996, -502.143476762872, -77.1352551735404},
    {32.6931013148151, -779.013864112386, 857.782474508566},
    {655.705330413471, -556.298002122365, 830.746155177827},
    {240.692048115634, 481.594127234753, 568.053743957541},
    {954.234301358783, 290.158081611087, 448.730905808833},
    {-120.906263528027, 697.783462448126, -225.448048687736},
    {656.020318513316, 928.709919713953, -495.741183969965},
    {-72.311712896204, 522.208236687249, -864.080266295935},
    {-349.192291548839, -33.3969047644982, 307.476526257641},
    {602.383156960101, -357.545693796819, 491.040763133293},
    {-816.48670479607, 690.014967324475, 156.832133521987},
    {-441.224803794051, 986.848459119069, 727.694279141012},
    {904.459903140815, 633.809316347028, 720.809477632895},
    {-311.989750030536, -751.96979473132, 969.254472738936},
    {-740.62579887492, 33.7489929682001, 635.189241679399},
    {778.316376532839, 117.432781330994, 326.423930726778},
    {220.485896632056, -712.515173070311, -308.645372825926},
    {-801.168624602233, 438.714563838592, -258.625219431671},
    {-414.287380493994, -553.568630355225, -744.850347506279},
    {420.071412458365, -273.144407011969, -241.21757354406},
    {-94.025223575396, -768.801582605663, -155.642182277114},
    {399.405122828651, -479.861602840649, -927.197983794348},
    {-269.32439046881, -783.573654929277, 70.4197032860075},
    {-45.4455869782096, -56.8334203670229, 499.956008429378},
    {318.431121052834, -459.747596075554, 275.981581604076},
    {343.159519967713, -962.375852452342, -220.938301068274},
    {-39.5692948027645, 803.25449144901, -337.992083087649},
    {-825.531281565063, -971.933976447876, -761.332590553888},
    {675.476874678064, -857.945683330058, -638.962350315903},
    {955.375641796014, -790.210808179765, -197.798242029332},
    {-25.3468304879225, -596.394891513412, -93.2674405075813},
    {847.967787651736, 276.973053608077, 225.701554278454},
    {435.805606942352, -420.956976825166, 900.00107477675},
    {836.6467129164, -638.584272284312, 684.8371315719},
    {-207.783343157221, 408.598528980288, -522.947803639255},
    {252.317805709171, -120.376797339276, 602.520311573526},
    {-605.799294430675, 303.783061377307, -531.762252284141},
    {530.830859043692, 394.206234693401, 883.54894212435},
    {181.386926024859, 535.832563478415, -78.2306358660899},
    {996.873626252225, -601.554683316418, 419.298266636961},
    {790.001458599247, 224.086605425205, -364.445403382046},
    {-628.818807952816, 972.37450645017, 768.545662681304},
    {203.634426237671, 785.606051868598, -946.484762079921},
    {814.03579159722, 946.294919198234, 189.841935207071},
    {-717.643469010524, -74.1010442879569, 937.284379257764},
    {470.520294036527, -639.161584916107, 3.46261394542307},
    {90.657965625649, -119.826900773027, -778.389416128107},
    {937.117114807943, 849.144909828119, 654.785726804968},
    {749.043722152951, 808.472185738992, -977.128535587145},
    {-691.985465021816, 54.8295853487343, -231.054020636126},
    {237.969465185155, 983.029171433374, -199.3422098228},
    {62.8579086469681, 608.785916418915, 738.326854459946},
    {334.448508751836, -604.890750276757, -30.0614542511862},
    {603.110312320513, 595.569473389423, -212.457700625955},
    {33.3070738076726, 148.298013948962, -796.81715487436},
    {-996.244804536375, 843.260297536901, 220.914438842369},
    {-525.795563319543, 344.271560151409, -56.2588684529319},
    {843.779487628062, -788.37354965172, 984.157726007582},
    {95.6440476594312, -22.9414101747743, -757.724858192136},
    {723.93351334992, 349.44799203713, 328.362786825903},
    {381.140115583364, -62.6490666151053, -200.347727595031},
    {981.827122068397, 318.58185835718, -243.859285387532},
    {92.4247100799175, 243.407219713103, 303.835774481205},
    {-756.003713091903, 981.942737322908, -71.3384348720747},
    {-561.936945920006, 815.263585809571, -28.1114106419209},
    {-262.455204154475, 273.986475828405, -81.1778766970124},
    {464.540095015497, -632.474169120966, -780.706448448065},
    {-674.130510017565, -992.721824229168, -11.621873808795}};

Vector test_vectors_2[n_test_vectors] = {
    {-940.167557040267, 969.239910381154, -778.182701473645},
    {862.180425861932, 791.948163236398, -546.898381378426},
    {536.357917942021, -986.013783737798, 298.205373880619},
    {813.494857890366, -159.321596800805, 367.103581376094},
    {632.047877594427, -483.08979361665, 111.825280063864},
    {931.880756770541, 998.57504374245, -955.01453049844},
    {872.5178359243, 789.520408339006, -846.083353197654},
    {-366.677128976043, 315.460292151844, -798.848094232703},
    {599.446469600386, 713.353630243932, -471.466480409756},
    {-182.990847486717, -25.3935591983282, -862.804456333813},
    {-151.770150151618, 88.8486834035157, 334.509601991639},
    {-130.240694307841, 862.82115153867, 703.180112964144},
    {-273.969167611181, -117.953294380707, 880.792495462638},
    {742.638159607509, 376.446976583916, 229.400517848007},
    {-535.055434343185, 269.276942885158, 193.414487404162},
    {-799.080898724065, 202.428827576784, 537.485159448213},
    {-32.2477189042443, 386.177196402807, 530.610044040718},
    {-165.867803907543, -976.503420874605, 39.7641828531409},
    {-723.901213222794, -217.230616661169, -519.66887651918},
    {883.889922718964, 22.855826456564, 210.640616251361},
    {319.829132808818, -949.288871932236, -854.882306686849},
    {427.15965103162, 894.218255631375, 27.0125536572264},
    {-5.50973767516916, 564.172893735883, 630.248557115008},
    {-409.039430944384, -74.0498837292571, -961.14368068365},
    {216.933823243398, -650.542942362347, 103.2825382333},
    {579.40043261344, -56.9494991176184, 739.280223083144},
    {-78.1176873849581, -481.780725186198, -363.216652013093},
    {610.755307241353, -232.017278125308, 940.222540552195},
    {-277.144344962387, 409.314341771651, 420.995043884519},
    {-106.429310559239, 142.132114636629, 592.623527593521},
    {29.8120230155701, -788.526179156102, 681.076180515649},
    {84.7725907814238, 69.6166542267733, -304.005105084687},
    {381.107762588089, 268.894319663489, -499.241217636834},
    {81.7846476626914, 147.346762701901, -200.866007090181},
    {-326.76053933427, 133.477307102166, 740.78619687146},
    {-677.595471379916, 288.229395612953, 370.499691457368},
    {-893.667554507131, -147.862498030393, -808.438400979587},
    {348.009113955909, 533.279979983849, 589.048825845067},
    {-676.29908744833, 800.116130461478, 683.44409747257},
    {800.145079500341, -878.420663886318, -634.470455992802},
    {6.56042121261726, -498.785360174911, -820.559969064393},
    {458.695914032837, -752.622781508641, -79.9367467105467},
    {802.473470291352, 193.991580134009, -64.3295278810338},
    {-903.254187806833, 224.750389775827, 275.714325376408},
    {-596.116589219191, -293.295438791498, 256.583868234164},
    {454.611717186319, -928.123196088678, -634.505160903272},
    {625.605510582961, -439.807021629766, -88.9923936445771},
    {931.335973621152, 842.447854686061, -628.226927712498},
    {884.267712973644, 219.238659656943, 222.824445195702},
    {65.9208505131628, -679.286363791278, 802.985807846148},
    {-229.327668105674, -205.363394556086, 842.995112678081},
    {-280.904621773856, -485.259908858634, 678.161930546956},
    {-74.198482968864, -731.613274619312, 282.93607088718},
    {-692.584200333823, 304.1703880285, -232.011425471716},
    {738.471632098248, -593.946506394125, 817.959701163322},
    {450.820317869711, -161.124638672419, 430.248582223994},
    {711.676735123729, 109.084258810335, -941.928491730718},
    {790.945195351691, 586.697411103935, -38.5512464935264},
    {-103.165537310829, -7.09246713084258, -878.063040841682},
    {97.0930542765682, -196.04215209696, 103.736839981131},
    {-965.509263195671, 494.455097989716, -206.720740467556},
    {2.34384008544748, 915.611017683457, 464.039469927941},
    {-990.693560542068, 998.224506923018, -496.320742916367},
    {721.299339549281, -585.432197488436, 288.709821699094},
    {414.952336902761, -333.862452596894, 348.06802007226},
    {-235.665820151885, -269.83008832476, 861.387387576548},
    {-684.302434063994, -324.108559600555, 138.935392254373},
    {-659.111703956274, 960.756359751708, -0.924853241247092},
    {596.330357252415, -547.81758079155, 690.540031966885},
    {140.221423798001, 512.217808323378, 347.146326255081},
    {732.295370230584, 951.055760308359, -767.051781237123},
    {-192.505413907404, -479.829552462776, 573.835824028247},
    {172.249741391946, -91.8234331606063, -392.683791622662},
    {625.421528430271, -971.582794441829, 272.270556450204},
    {250.761091762223, -739.619336700663, -91.0398048307475},
    {776.258993052821, -243.498109242495, -78.7251736838784},
    {-804.942566743314, 162.017076076031, 185.815578946228},
    {330.095138141986, 15.8230627714902, -54.2120205776123},
    {-767.07153047659, -694.633962228167, 388.755084659394},
    {989.173784431953, 313.185361235576, -2.67099047253259},
    {496.120973139961, -618.282226030301, -52.0717992137461},
    {-904.084725896004, 841.012989798086, -185.58591408344},
    {-209.68246382652, 44.6564013394031, 799.936871544925},
    {719.939769200741, -262.58838620899, 620.124878524882},
    {-239.964349512489, -530.584443831631, 628.45158986313},
    {-372.273192596329, 195.30119838018, -264.687063992326},
    {-243.453890300539, 289.318188941005, 59.2685188079809},
    {909.841604020101, -107.501343185532, -223.696016211806},
    {925.377605924698, -57.0501550610275, 471.166106008848},
    {950.724730356117, 385.431171361511, -403.500523244587},
    {-376.829932317231, 61.9551999144833, -129.978883519932},
    {871.149726795444, 194.119195443202, -134.444360831203},
    {330.619480521278, -574.455968977393, -805.407449613446},
    {919.883981861718, 809.152948706385, 125.744091212729},
    {423.875178271289, 179.201608654147, -86.7897032360188},
    {-619.366016867911, -441.534810800491, 646.420705428553},
    {23.6845578256962, -796.145668432047, 944.605112747485},
    {-327.015723339448, 413.549001219022, 24.4981043314624},
    {-733.512838329998, -49.6073372788646, -93.6529429095194},
    {-397.193159803579, -103.921148183533, -260.52405523152}};

// End of file
