// lorentz.hh
// T. M. Kelley
// May 02, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#pragma once

#include "base/constants.h"
#include "mesh_common/Vector.h"
#include <cmath>

namespace murmeln_mesh {

/**\Compute relativistic gamma
 *\tparam FP_T floating point type
 */

template <typename FP_T, typename vector_t>
inline FP_T
gamma(vector_t const v)
{
  using murmeln::constants::SPEED_LIGHT;
  return 1.0 / std::sqrt(1 - v.dot(v) / (SPEED_LIGHT * SPEED_LIGHT));
}

namespace spec_3D_Cartesian {

template <typename FloatingPoint_T>
using EandOmega_T = murmeln_mesh::Vector4<FloatingPoint_T>;

/**\brief Transform massless particle energy and direction to lab frame
 * \tparam FP_T: a floating point type
 * \param v_lab: velocity of comoving material (in the lab frame!)
 * \param e_lab: lab particle energy
 * \param omega_lab: lab particle direction
 * \return pair with {e_comoving, direction_comoving } */
template <typename FP_T>
inline EandOmega_T<FP_T>
LT_to_comoving(Vector const & v_lab,
               FP_T const & e_lab,
               Vector const & omega_lab)
{
  using murmeln::constants::SPEED_LIGHT;
  // Require that omega_lab is a unit vector, |v_lab| < c
  EandOmega_T<FP_T> res;
  FP_T & e = res[0];
  Vector omega{res[1], res[2], res[3]};
  FP_T const vdo = v_lab.dot(omega_lab);
  FP_T const gam = gamma<FP_T>(v_lab);
  FP_T const goc = gam / SPEED_LIGHT;
  FP_T const fac = 1 - goc * vdo / (gam + 1);
  e = gam * e_lab * (1 - vdo / SPEED_LIGHT);
  FP_T const eoec = e_lab / e;  // E_lab/E_comving
  omega[0] = eoec * (omega_lab[0] - goc * v_lab[0] * fac);
  omega[1] = eoec * (omega_lab[1] - goc * v_lab[1] * fac);
  omega[2] = eoec * (omega_lab[2] - goc * v_lab[2] * fac);
  return {e, omega[0], omega[1], omega[2]};
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
template <typename FP_T>
EandOmega_T<FP_T> inline LT_to_lab(murmeln_mesh::Vector const & v_lab,
                                   FP_T const & e_com,
                                   murmeln_mesh::Vector const & omega_com)
{
  using murmeln::constants::SPEED_LIGHT;
  EandOmega_T<FP_T> res;
  FP_T const vdo = v_lab.dot(omega_com);
  FP_T const gam = gamma<FP_T>(v_lab);
  FP_T const goc = gam / SPEED_LIGHT;
  FP_T const fac = 1 + goc * vdo / (gam + 1);
  res[0] = gam * e_com * (1 + vdo / SPEED_LIGHT);
  FP_T const eoec = e_com / res[0];  // E/E_com
  res[1] = eoec * (omega_com[0] + goc * v_lab[0] * fac);
  res[2] = eoec * (omega_com[1] + goc * v_lab[1] * fac);
  res[3] = eoec * (omega_com[2] + goc * v_lab[2] * fac);
  return res;
}  // LT_to_lab

}  // namespace spec_3D_Cartesian

// specializations for 1D spherical
namespace spec_1D_Spherical {

template <typename FP_T>
using EandOmega_T = murmeln_mesh::Vector2<FP_T>;

using Vector1 = murmeln_mesh::Vector1;

template <typename FP_T>
EandOmega_T<FP_T> inline LT_to_comoving(Vector1 const v_lab,
                                        FP_T const e_lab,
                                        Vector1 const omega_lab)
{
  using murmeln::constants::SPEED_LIGHT;
  FP_T const voc = v_lab[0] / SPEED_LIGHT;
  FP_T const onemovoc = (1.0 - omega_lab[0] * voc);
  FP_T const e_com = e_lab * gamma<FP_T>(v_lab) * onemovoc;
  Vector1 omega_com;
  omega_com[0] = (omega_lab[0] - voc) / onemovoc;
  return EandOmega_T<FP_T>{e_com, omega_com[0]};
}  // LT_to_comoving_sphere1D

// Yes, you should pass v_lab here--it's the material speed in the lab frame!
template <typename FP_T>
EandOmega_T<FP_T> inline LT_to_lab(Vector1 const v_lab,
                                   FP_T const e_com,
                                   Vector1 const omega_com)
{
  using murmeln::constants::SPEED_LIGHT;
  FP_T const voc = v_lab[0] / SPEED_LIGHT;
  FP_T const onepovoc = 1.0 + omega_com[0] * voc;
  FP_T const e_lab = e_com * gamma<FP_T>(v_lab) * onepovoc;
  Vector1 omega_lab;
  omega_lab[0] = (voc + omega_com[0]) / onepovoc;
  return EandOmega_T<FP_T>(e_lab, omega_lab[0]);
}  // LT_to_lab_sphere1D

}  // namespace spec_1D_Spherical

}  // namespace murmeln_mesh

// End of file
