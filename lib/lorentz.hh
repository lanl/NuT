// lorentz.hh
// T. M. Kelley
// May 02, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#ifndef LORENTZ_HH
#define LORENTZ_HH

#include "constants.hh"
// #include "Vec3D.hh"
#include "types.hh"
#include <cmath>

namespace nut {

template <size_t dim>
using EandOmega = std::pair<geom_t, vec_t<dim> >;

geom_t inline gamma(geom_t const v)
{
  return 1.0 / std::sqrt(1 - v * v / (c * c));
}

template <size_t dim>
geom_t inline gamma(vec_t<dim> const v)
{
  return 1.0 / std::sqrt(1 - dot(v, v) / (c * c));
}

namespace spec_3D_Cartesian {

EandOmega<3> inline LT_to_comoving(vec_t<3> const & v_lab,
                                   geom_t const & e_lab,
                                   vec_t<3> const & omega_lab)
{
  EandOmega<3> res;
  geom_t & e = res.first;
  vec_t<3> & omega = res.second;
  geom_t const vdo = dot(v_lab, omega_lab);
  geom_t const gam = gamma(v_lab);
  geom_t const goc = gam / c;
  geom_t const fac = 1 - goc * vdo / (gam + 1);
  e = gam * e_lab * (1 - vdo / c);
  geom_t const eoec = e_lab / e;  // E_lab/E_comving
  omega.v[0] = eoec * (omega_lab.v[0] - goc * v_lab.v[0] * fac);
  omega.v[1] = eoec * (omega_lab.v[1] - goc * v_lab.v[1] * fac);
  omega.v[2] = eoec * (omega_lab.v[2] - goc * v_lab.v[2] * fac);
  return res;
}  // LT_to_comoving

EandOmega<3> inline LT_to_lab(vec_t<3> const & v_lab,
                              geom_t const & e_com,
                              vec_t<3> const & omega_com)
{
  EandOmega<3> res;
  geom_t const vdo = dot(v_lab, omega_com);
  geom_t const gam = gamma(v_lab);
  geom_t const goc = gam / c;
  geom_t const fac = 1 + goc * vdo / (gam + 1);
  res.first = gam * e_com * (1 + vdo / c);
  geom_t const eoec = e_com / res.first;  // E/E_com
  res.second.v[0] = eoec * (omega_com.v[0] + goc * v_lab.v[0] * fac);
  res.second.v[1] = eoec * (omega_com.v[1] + goc * v_lab.v[1] * fac);
  res.second.v[2] = eoec * (omega_com.v[2] + goc * v_lab.v[2] * fac);
  return res;
}  // LT_to_lab

}  // namespace spec_3D_Cartesian

namespace spec_1D {

EandOmega<1> inline LT_to_comoving_sphere1D(geom_t const v_lab,
                                            geom_t const e_lab,
                                            vec_t<1> const omega_lab)
{
  geom_t const voc = v_lab / c;
  geom_t const onemovoc = (1.0 - omega_lab.v[0] * voc);
  geom_t const e_com = e_lab * gamma(v_lab) * onemovoc;
  vec_t<1> omega_com;
  omega_com.v[0] = (omega_lab.v[0] - voc) / onemovoc;
  return EandOmega<1>(e_com, omega_com);
}  // LT_to_comoving_sphere1D

// Yes, you should pass v_lab here--it's the material speed in the lab frame!
EandOmega<1> inline LT_to_lab_sphere1D(geom_t const v_lab,
                                       geom_t const e_com,
                                       vec_t<1> const omega_com)
{
  geom_t const voc = v_lab / c;
  geom_t const onepovoc = 1.0 + omega_com.v[0] * voc;
  geom_t const e_lab = e_com * gamma(v_lab) * onepovoc;
  vec_t<1> omega_lab;
  omega_lab.v[0] = (voc + omega_com.v[0]) / onepovoc;
  return EandOmega<1>(e_lab, omega_lab);
}  // LT_to_lab_sphere1D

}  // namespace spec_1D

}  // namespace nut

#endif  // LORENTZ_HH

// version
// $Id$

// End of file
