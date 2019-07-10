// lorentz.hh
// May 02, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#pragma once

#include "Assert.hh"
#include "constants.hh"
#include "soft_equiv.hh"
#include <cmath>

namespace nut {

/**\Compute relativistic factor gamma
 *\tparam FP_T floating point type
 */
template <typename FP_T, typename vector_t>
inline FP_T
gamma(vector_t const v)
{
  return 1.0 / std::sqrt(1 - v.dot(v) / (nut::c * nut::c));
}

/* Lorentz tranforms for 3D cartesian */
namespace spec_3D_Cartesian {

template <typename FP_T, typename vector_t>
using EandOmega_T = std::pair<FP_T, vector_t>;

/**\brief Transform massless particle energy and direction to lab frame
 * \tparam FP_T: a floating point type
 * \param v_lab: velocity of comoving material (in the lab frame!)
 * \param e_lab: lab particle energy
 * \param omega_lab: lab particle direction
 * \return pair with {e_comoving, direction_comoving } */
template <typename FP_T, typename vector_t>
inline EandOmega_T<FP_T, vector_t>
LT_to_comoving(vector_t const & v_lab,
               FP_T const & e_lab,
               vector_t const & omega_lab)
{
  nut::Equal(vector_t::dim, 3, "vector dimension", "3");
  nut::Require(soft_equiv(1.0, omega_lab.norm()),
               "direction must be unit vector!");
  nut::LessThan(v_lab.norm(), nut::c, "material speed");
  EandOmega_T<FP_T, vector_t> res;
  FP_T & e = res.first;
  vector_t & omega = res.second;
  FP_T const vdo = v_lab.dot(omega_lab);
  FP_T const gam = gamma<FP_T>(v_lab);
  FP_T const goc = gam / nut::c;
  FP_T const fac = 1 - goc * vdo / (gam + 1);
  e = gam * e_lab * (1 - vdo / nut::c);
  FP_T const eoec = e_lab / e;  // E_lab/E_comving
  omega[0] = eoec * (omega_lab[0] - goc * v_lab[0] * fac);
  omega[1] = eoec * (omega_lab[1] - goc * v_lab[1] * fac);
  omega[2] = eoec * (omega_lab[2] - goc * v_lab[2] * fac);
  return res;
}  // LT_to_comoving

/**\brief Transform massless particle energy and direction to lab frame
 * \tparam FP_T: a floating point type
 * \param v_lab: velocity of comoving material (in the lab frame!)
 * \param e_com: comoving particle energy
 * \param omega_com: comoving particle direction
 * \return pair with {e_lab, direction_lab }
 *
 * Note that these are naive implementations. It might make sense to switch
 * to a rapidity-based imeplementation at some energy scale? */
template <typename FP_T, typename vector_t>
inline EandOmega_T<FP_T, vector_t>
LT_to_lab(vector_t const & v_lab,
          FP_T const & e_com,
          vector_t const & omega_com)
{
  nut::Equal(vector_t::dim, 3, "vector dimension", "3");
  nut::Require(soft_equiv(1.0, omega_com.norm()),
               "direction must be unit vector!");
  nut::LessThan(v_lab.norm(), nut::c, "material speed");
  EandOmega_T<FP_T, vector_t> res;
  FP_T const vdo = v_lab.dot(omega_com);
  FP_T const gam = gamma<FP_T>(v_lab);
  FP_T const goc = gam / nut::c;
  FP_T const fac = 1 + goc * vdo / (gam + 1);
  res.first = gam * e_com * (1 + vdo / nut::c);
  FP_T const eoec = e_com / res.first;  // E/E_com
  res.second[0] = eoec * (omega_com[0] + goc * v_lab[0] * fac);
  res.second[1] = eoec * (omega_com[1] + goc * v_lab[1] * fac);
  res.second[2] = eoec * (omega_com[2] + goc * v_lab[2] * fac);
  return res;
}  // LT_to_lab

}  // namespace spec_3D_Cartesian

// specializations for 1D spherical
namespace spec_1D_Spherical {

template <typename FP_T, typename vector_t>
using EandOmega_T = std::pair<FP_T, vector_t>;

template <typename FP_T, typename vector_t>
inline EandOmega_T<FP_T, vector_t>
LT_to_comoving_sphere1D(vector_t const v_lab,
                        FP_T const e_lab,
                        vector_t const omega_lab)
{
  nut::Equal(vector_t::dim, size_t{1u}, "vector dimension", "1");
  nut::Require(nut::soft_equiv(1.0, omega_lab.norm()),
               "direction must be unit vector!");
  nut::LessThan(v_lab.norm(), nut::c, "material speed");
  FP_T const voc = v_lab[0] / nut::c;
  FP_T const onemovoc = (1.0 - omega_lab[0] * voc);
  FP_T const e_com = e_lab * gamma<FP_T>(v_lab) * onemovoc;
  vector_t omega_com;
  omega_com[0] = (omega_lab[0] - voc) / onemovoc;
  return EandOmega_T<FP_T, vector_t>(e_com, omega_com);
}  // LT_to_comoving_sphere1D

// Yes, you should pass v_lab here--it's the material speed in the lab frame!
template <typename FP_T, typename vector_t>
inline EandOmega_T<FP_T, vector_t>
LT_to_lab_sphere1D(vector_t const v_lab, FP_T const e_com, vector_t const omega_com)
{
  nut::Equal(vector_t::dim, size_t{1u}, "vector dimension", "1");
  nut::Require(soft_equiv(1.0, omega_com.norm()),
               "direction must be unit vector!");
  nut::LessThan(v_lab.norm(), nut::c, "material speed");
  FP_T const voc = v_lab[0] / nut::c;
  FP_T const onepovoc = 1.0 + omega_com[0] * voc;
  FP_T const e_lab = e_com * gamma<FP_T>(v_lab) * onepovoc;
  vector_t omega_lab;
  omega_lab[0] = (voc + omega_com[0]) / onepovoc;
  return EandOmega_T<FP_T, vector_t>(e_lab, omega_lab);
}  // LT_to_lab_sphere1D

}  // namespace spec_1D_Spherical

}  // namespace nut

// End of file
