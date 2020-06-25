// Spherical_Mesh.i.h
// Jun 25, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#ifndef IM_ALLOWED_TO_INCLUDE_Spherical_MESH_I_H
#error "Use Spherical_Mesh.hh---do not include Spherical_Mesh.i.hh"
#endif

#include "base/constants.h"

namespace murmeln_mesh {

inline Spherical_1D_Mesh::Geom_T
Spherical_1D_Mesh::volume(Cell const cell) const
{
  // cellOK(cell);
  index_t const index(cell.as_id());
  geom_t const lo = m_cell_bounds.at(index - 1);
  geom_t const hi = m_cell_bounds.at(index);
  geom_t const vol =
      4. / 3. * murmeln::constants::pi * (hi * hi * hi - lo * lo * lo);
  // GreaterThan(vol, geom_t(0.0), "volume");
  return vol;
}  // volume

template <typename RNG_T>
Spherical_1D_Mesh::Vector
Spherical_1D_Mesh::sample_position(RNG_T & rng, Cell const & c) const
{
  auto [lo, hi] = this->get_extents(c);
  // Geom_T const lo = extents.first;
  // Geom_T const hi = extents.second;
  Geom_T const lo3 = lo * lo * lo;
  Geom_T const hi3 = hi * hi * hi;
  Geom_T const r = std::pow(lo3 + (hi3 - lo3) * rng.random(), 1.0 / 3.0);
  return Vector{r};
}  // sample_position

template <typename RNG_T>
Spherical_1D_Mesh::Vector
Spherical_1D_Mesh::sample_direction_isotropic(RNG_T & r)
{
  geom_t const ctheta = geom_t(2) * r.random() - geom_t(1);
  return {ctheta};
}  // sample_direction_isotropic

inline Spherical_1D_Mesh::Intersection_T
Spherical_1D_Mesh::intersection(Ray const & r, Cell const & c) const
{
  auto [rlo, rhi] = this->get_extents(c);
  geom_t const x{r.position()[0]};
  geom_t const omega{r.direction()[0]};
  index_t cidx{c.as_id()};
  return dist_to_bdy_impl(x, omega, rlo, rhi, cidx);
}

inline Spherical_1D_Mesh::Intersection_T
Spherical_1D_Mesh::dist_to_bdy_impl(geom_t x,
                                    geom_t omega,
                                    geom_t rlo,
                                    geom_t rhi,
                                    index_t cell_idx)
{
  using murmeln::soft_equiv;
  using murmeln::constants::Max_Dbl;
  using murmeln::constants::pi;
  geom_t const rhisq = rhi * rhi;
  geom_t const rlosq = rlo * rlo;
  geom_t const xsq = x * x;

  // need to (1-2) check for intersection with the two spheres that
  // bound the current cell, and (3) choose the first intersection
  // along the direction of travel.

  // 1. Compute intersections with outer sphere

  // get Tan(theta) from inverting omega=Cos(theta). We work with a
  // reduced polar angle, "abs_theta", in the range 0 <= theta <= pi/2,
  // then select the correct intercept between the line that the
  // particle travels and the spherical shell using the full value of
  // theta.
  geom_t const theta = std::acos(omega);
  geom_t const abs_theta = (theta > pi / 2) ? pi - theta : theta;
  geom_t const t = std::tan(abs_theta);
  geom_t const tsq = t * t;
  geom_t const tcub = tsq * t;
  geom_t const one_plus_tsq = 1 + tsq;
  geom_t const one_on_one_plus_tsq = 1.0 / one_plus_tsq;
  // the determinant is r^2 + r^2 * t^2 - x^2 * t^2. Include the
  // square root and divide by 1+Tan(theta)^2 for convenience here.
  geom_t const dethi =
      std::sqrt(rhisq * one_plus_tsq - xsq * tsq) * one_on_one_plus_tsq;
  // These parts don't depend on the radius of the sphere.
  geom_t const xterm1 = x * tsq * one_on_one_plus_tsq;
  geom_t const yterm1 = -x * t + x * tcub * one_on_one_plus_tsq;
  // These are the two intersection with the outer sphere
  geom_t const xhip = xterm1 + dethi;
  geom_t const yhip = yterm1 + t * dethi;
  geom_t const xhim = xterm1 - dethi;
  geom_t const yhim = yterm1 - t * dethi;
  // if the polar angle is less than pi/2, we want the solution
  // that adds the determinant.
  geom_t const xhi = (theta < pi / 2) ? xhip : xhim;
  geom_t const yhi = (theta < pi / 2) ? yhip : yhim;
  geom_t const d_hi = std::sqrt((xhi - x) * (xhi - x) + (yhi * yhi));
  // 2. Look for intersection with inner sphere
  // There is if the absolute value of the polar angle is
  // less than atan(1/(b^2 -1)), where b = x/r_lo. In this case, we can
  // just use the "plus" solution, since that's always the closest.
  // Note that you may want to know a real solution exists before
  // computing it--complicates branch removal.
  geom_t d_lo = Max_Dbl;
  if(rlo > 0.0) {
    geom_t const b = x / rlo;
    // if the particle is on the inner sphere, and headed inward...
    if(soft_equiv(b, 1.0, 1e-11)) {
      if(theta > pi / 2) { d_lo = 0.0; }
    }
    else {
      geom_t const tan_theta_lim = std::sqrt(1 / (b * b - 1));
      geom_t const theta_lim = std::atan(tan_theta_lim);
      if(abs_theta <= theta_lim and theta > pi / 2) {
        geom_t const detlo =
            std::sqrt(rlosq * one_plus_tsq - xsq * tsq) * one_on_one_plus_tsq;
        geom_t const xlo = xterm1 + detlo;
        geom_t const ylo = yterm1 + t * detlo;
        d_lo = std::sqrt((xlo - x) * (xlo - x) + ylo * ylo);
      }
    }
  }
  // now select the shorter distance and the corresponding face
  geom_t d_bdy = Max_Dbl;
  Face face{null_face_};
  if(d_hi < d_lo) {
    d_bdy = d_hi;
    face = Face{cell_idx};  // will intersect outer sphere
  }
  else {
    d_bdy = d_lo;
    face = Face{cell_idx - 1};  // will intersect inner sphere
  }
  Intersection_T d2b{face, d_bdy};
  return d2b;

}  // dist_to_bdy_impl

inline Spherical_1D_Mesh::Face
Spherical_1D_Mesh::cell_to_face(Cell const & c, Face_Name fn) const
{
  index_t the_index{c.as_id()};
  if(fn == HIGH) { the_index++; }
  return Face{the_index};
}

inline Spherical_1D_Mesh::Ray
Spherical_1D_Mesh::advance_point(Point const & x,
                                 Vector const & dir,
                                 double const dist) const
{
  geom_t const theta = std::acos(dir[0]);
  geom_t const s = std::sin(theta);
  geom_t const new_x = x[0] + dist * dir[0];
  geom_t const new_y = dist * s;
  geom_t const new_r = std::sqrt(new_x * new_x + new_y * new_y);
  geom_t const new_omega{std::cos(theta - std::asin(dist / new_r * s))};
  Ray coord{{new_r}, {new_omega}};
  return coord;
}

inline Spherical_1D_Mesh::Cell
Spherical_1D_Mesh::find_cell(Point const & p) const
{
  // search through bin boundaries to find the pair that bracket p
  auto const it =
      std::lower_bound(m_cell_bounds.begin(), m_cell_bounds.end(), p[0]);
  index_t const cellidx{static_cast<index_t>(it - m_cell_bounds.begin())};
  Cell c{cellidx};
  // Require(in_cell(p, c));
  return c;
}

inline bool
Spherical_1D_Mesh::in_cell(Vector const & p, Cell const & c) const
{
  auto [lo, hi] = get_extents(c);
  return (lo < p[0] && hi > p[0]);
}

inline bool
Spherical_1D_Mesh::is_low_face(Cell const & c, Face const & f)
{
  index_t face_id = f.as_id();
  index_t cell_id = c.as_id();
  return face_id == cell_id;
}

inline Spherical_1D_Mesh::Cell
Spherical_1D_Mesh::cell_across(Cell const & c,
                               Face const & f,
                               Point const & /*p*/) const
{
  index_t const cidx{c.as_id()};
  bool const face_is_low_face{is_low_face(c, f)};
  if(face_is_low_face) {
    if(cidx == 0) { return null_cell_; }
    Cell c_ret{cidx - 1};
    // Check(c.as_id() < null_cell.as_id());
    return c_ret;
  }
  if(cidx >= m_num_cells - 1) { return null_cell_; }
  Cell c_ret{cidx + 1};
  // Check(c.as_id() < null_cell_.as_id());
  return c_ret;
}  // cell_across

inline bool
Spherical_1D_Mesh::is_boundary(Face f) const
{
  bool const is_bdy_face{f.as_id() == 0ull || f.as_id() == m_num_cells};
  return is_bdy_face;
}

inline Spherical_1D_Mesh::Extents_T
Spherical_1D_Mesh::get_extents(Cell const & c) const
{
  index_t const cidx{c.as_id()};
  return {m_cell_bounds[cidx - 1], m_cell_bounds[cidx]};
}

inline index_t
Spherical_1D_Mesh::num_cells() const
{
  return m_num_cells;
}

inline geom_t
Spherical_1D_Mesh::get_distance(Intersection_T const & i)
{
  return i.second;
}

inline Spherical_1D_Mesh::Face
Spherical_1D_Mesh::get_face(Intersection_T const & i)
{
  return i.first;
}

inline geom_t
Spherical_1D_Mesh::get_low(Extents_T const & e)
{
  return e.first;
}

inline geom_t
Spherical_1D_Mesh::get_high(Extents_T const & e)
{
  return e.second;
}

inline Spherical_1D_Mesh::Vector
Spherical_1D_Mesh::reflect(Face const & /*f*/, Vector const & v)
{
  return Vector{-v[0]};
}

inline Spherical_1D_Mesh::Vector
Spherical_1D_Mesh::get_normal(Cell const & c, Face const & f)
{
  if(is_low_face(c, f)) { return Vector{1.0}; }
  return Vector{-1.0};
}

}  // namespace murmeln_mesh

// End of file
