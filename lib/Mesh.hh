// Mesh.hh
// T. M. Kelley
// Jan 11, 2011
// Header for Mesh
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef MESH_H
#define MESH_H

#include "Assert.hh"
#include "constants.hh"
#include "detail/Vector.h"
#include "meshes/geometry/lorentz.h"
#include "soft_equiv.hh"
#include "types.hh"
#include "utilities_io.hh"
#include <cmath>
#include <string>
#include <vector>

namespace nut {

namespace Sphere_1D_Faces {
enum Faces { LOW = 0, HIGH = 1, NULL_FACE };
}

///*!\brief Mesh functions for 1D spherical geometry.
// * \tparam <cell_t> {cell index type}
// * \tparam <boundary_t> {geometry (numerical) type}
// * \tparam <bdy_descriptor_t> {boundary descriptor type}
// */
//template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
//struct Sphere_1D {
//public:
//  static const size_t dim = 1;
//
//  typedef geometry_t geom_t;
//  typedef bdy_descriptor_t bdy_desc_t;
//  typedef std::vector<geom_t> vb;
//  typedef std::vector<bdy_desc_t> vbd;
//  typedef std::pair<geom_t, geom_t> extents_t;
//
//  using Vector = Vector1;
//  using face_handle_t = id_t;
//  using Cell = cell_t;
//  using vector4_t = Vector2;
//
//  struct coord_t {
//    Vector1 x;
//    Vector1 omega;
//
//    Vector1 const & position() const { return x; }
//
//    Vector1 const & direction() const { return omega; }
//  };
//  struct d_to_b_t {
//    geom_t d;    /** distance to boundary */
//    cell_t face; /** which face will be intersected */
//  };
//
//  using intersect_t = d_to_b_t;
//  using Ray = coord_t;
//
//  using Boundary = bdy_descriptor_t;
//
//  static geom_t get_distance(intersect_t const & i) { return i.d; }
//
//  static face_handle_t get_face(intersect_t const & i)
//  {
//    auto fidx = i.face;
//    return static_cast<face_handle_t>(fidx);
//  }
//
//  static const cell_t null_cell_;
//
//  static cell_t null_cell() { return null_cell_; }
//
//  static const face_handle_t null_face_;
//
//  static face_handle_t null_face() { return null_face_; }
//
//  // ctor
//  Sphere_1D(vb const & bdys_, vbd const & descs_)
//      : m_bdys(bdys_), m_descs(descs_), m_ncells(m_bdys.size() - 1)
//  {
//    nut::Require(m_bdys.size() >= 2, "must have at least two boundaries");
//    nut::Equal(m_bdys.size(), m_descs.size(), "bdys size", "descs size");
//  }
//
//  /*!\brief volume of spherical shell 'cell'. 0 < cell <= n_cells. */
//  geom_t volume(cell_t const cell) const
//  {
//    cellOK(cell);
//    cell_t const index(cell - 1);
//    geom_t const lo = m_bdys.at(index);
//    geom_t const hi = m_bdys.at(index + 1);
//    geom_t const vol = 4. / 3. * pi * (hi * hi * hi - lo * lo * lo);
//    nut::GreaterThan(vol, geom_t(0), "volume");
//    return vol;
//  }  // volume
//
//  /**\brief Compute normal to face */
//  static inline Vector1 get_normal(face_handle_t const & f)
//  {
//    if(Sphere_1D_Faces::LOW == f) { return Vector1{1.0}; }
//    return Vector1{-1.0};
//  }
//
//  /**\brief Reflection in spherical just reverses the direction cosine.*/
//  static inline Vector1 reflect(face_handle_t const & /*f*/,
//                                Vector1 const & direction,
//                                Vector1 const & /*unused*/)
//  {
//    return Vector1{-1.0 * direction[0]};
//  }
//  /**\brief For compatibility*/
//  static inline Vector1 reflect(Vector1 const & direction,
//                                face_handle_t const & /*f*/,
//                                Vector1 const & /*unused*/)
//  {
//    return Vector1{-1.0 * direction[0]};
//  }
//
//  /*!\brief which cell is across a face. Note this just tells what the next
//   * cell is, no account is taken of the boundary type.
//   *
//   * \param cell: 0 < cell <= n_cells
//   * \param face: 0 (left) or 1 (right)
//   */
//  cell_t cell_across(cell_t const cell,
//                     cell_t face,
//                     Vector1 const & /*p*/) const
//  {
//    cellOK(cell);
//    // LessThan(face, cell_t(2), "face");
//    if(cell == 1 && face == Sphere_1D_Faces::LOW) { return null_cell_; }
//    else if(cell == m_ncells && face == Sphere_1D_Faces::HIGH) {
//      return null_cell_;
//    }
//    if(face == Sphere_1D_Faces::LOW) { return cell - 1; }
//    return cell + 1;
//  }  // cell_across_face
//
//  template <typename rng_t>
//  geom_t sample_position(rng_t & rng, cell_t const cell) const
//  {
//    extents_t extents = this->cell_extents(cell);
//    geom_t const lo = extents.first;
//    geom_t const hi = extents.second;
//    geom_t const lo3 = lo * lo * lo;
//    geom_t const hi3 = hi * hi * hi;
//    geom_t const r = std::pow(lo3 + (hi3 - lo3) * rng.random(), 1.0 / 3.0);
//    return r;
//  }
//
//  template <typename RNG_T>
//  Vector1 static sample_direction_isotropic(RNG_T & rng)
//  {
//    geom_t const ctheta = geom_t(2) * rng.random() - geom_t(1);
//    return {ctheta};
//  }
//
//  /*!\brief type of cell boundary */
//  bdy_desc_t get_boundary_type(cell_t const c, cell_t const face) const
//  {
//    cell_t const idx = make_idx(c, m_ncells);
//    cell_t const fidx = face_to_index(face);
//    return m_descs[idx + fidx];
//  }
//
//  cell_t num_cells() const { return m_ncells; }
//
//  /*!\brief get the lower and upper bounds of a cell. *
//   * 0 < cell <= n_cells                              */
//  extents_t cell_extents(cell_t const cell) const
//  {
//    cellOK(cell);
//    return extents_t(m_bdys[cell - 1], m_bdys[cell]);
//  }
//
//  /*!\brief compute the face part of the index into m_descs. */
//  cell_t face_to_index(cell_t const f) const
//  {
//    LessThan(f, 2u, "Mesh1D::face_to_index: face");
//    return f;
//  }
//
//  Ray advance_point(Vector1 const x,
//                    Vector1 const omega,
//                    geom_t const distance) const
//  {
//    return this->new_coordinate(x, omega, distance);
//  }
//  /*!\brief calculate new coordinate and new direction cosine at a given
//   *        distance along direction cosine omega.  */
//  Ray new_coordinate(Vector1 const x,
//                     Vector1 const omega,
//                     geom_t const distance) const
//  {
//    geom_t const theta = std::acos(omega[0]);
//    geom_t const s = std::sin(theta);
//    geom_t const new_x = x[0] + distance * omega[0];
//    geom_t const new_y = distance * s;
//    geom_t const new_r = std::sqrt(new_x * new_x + new_y * new_y);
//    coord_t coord;
//    coord.x[0] = new_r;
//    coord.omega[0] = std::cos(theta - std::asin(distance / new_r * s));
//    return coord;
//  }
//
//  intersect_t intersection(Ray const & r, cell_t const c) const
//  {
//    Vector1 const & x{r.position()};
//    Vector1 const & o{r.direction()};
//    return this->distance_to_bdy(x, o, c);
//  }
//
//  d_to_b_t distance_to_bdy(Vector1 const x,
//                           Vector1 const omega,
//                           cell_t const cell) const
//  {
//    cellOK(cell);
//    extents_t extents = this->cell_extents(cell);
//    return dist_to_bdy_impl(x[0], omega[0], extents.first, extents.second,
//                            cell);
//  }  // distance_to_bdy
//
//  static d_to_b_t dist_to_bdy_impl(geom_t const x,
//                                   geom_t const omega,
//                                   geom_t const rlo,
//                                   geom_t const rhi,
//                                   cell_t const cell)
//  {
//    geom_t const rhisq = rhi * rhi;
//    geom_t const rlosq = rlo * rlo;
//    geom_t const xsq = x * x;
//
//    // need to (1-2) check for intersection with the two spheres that
//    // bound the current cell, and (3) choose the first intersection
//    // along the direction of travel.
//
//    // 1. Compute intersections with outer sphere
//
//    // get Tan(theta) from inverting omega=Cos(theta). We work with a
//    // reduced polar angle, "abs_theta", in the range 0 <= theta <= pi/2,
//    // then select the correct intercept between the line that the
//    // particle travels and the spherical shell using the full value of
//    // theta.
//    geom_t const theta = std::acos(omega);
//    geom_t const abs_theta = (theta > pi / 2) ? pi - theta : theta;
//    geom_t const t = std::tan(abs_theta);
//    geom_t const tsq = t * t;
//    geom_t const tcub = tsq * t;
//    geom_t const one_plus_tsq = 1 + tsq;
//    geom_t const one_on_one_plus_tsq = 1.0 / one_plus_tsq;
//    // the determinant is r^2 + r^2 * t^2 - x^2 * t^2. Include the
//    // square root and divide by 1+Tan(theta)^2 for convenience here.
//    geom_t const dethi =
//        std::sqrt(rhisq * one_plus_tsq - xsq * tsq) * one_on_one_plus_tsq;
//    // These parts don't depend on the radius of the sphere.
//    geom_t const xterm1 = x * tsq * one_on_one_plus_tsq;
//    geom_t const yterm1 = -x * t + x * tcub * one_on_one_plus_tsq;
//    // These are the two intersection with the outer sphere
//    geom_t const xhip = xterm1 + dethi;
//    geom_t const yhip = yterm1 + t * dethi;
//    geom_t const xhim = xterm1 - dethi;
//    geom_t const yhim = yterm1 - t * dethi;
//    // if the polar angle is less than pi/2, we want the solution
//    // that adds the determinant.
//    geom_t const xhi = (theta < pi / 2) ? xhip : xhim;
//    geom_t const yhi = (theta < pi / 2) ? yhip : yhim;
//    geom_t const d_hi = std::sqrt((xhi - x) * (xhi - x) + (yhi * yhi));
//
//    // 2. Look for intersection with inner sphere
//    // There is if the absolute value of the polar angle is
//    // less than atan(1/(b^2 -1)), where b = x/r_lo. In this case, we can
//    // just use the "plus" solution, since that's always the closest.
//    // Note that you may want to know a real solution exists before
//    // computing it--complicates branch removal.
//    geom_t d_lo = huge;
//    if(rlo > 0.0) {
//      geom_t const b = x / rlo;
//      // if the particle is on the inner sphere, and headed inward...
//      if(soft_equiv(b, 1.0, 1e-11)) {
//        if(theta > pi / 2) { d_lo = 0.0; }
//      }
//      else {
//        geom_t const tan_theta_lim = std::sqrt(1 / (b * b - 1));
//        geom_t const theta_lim = std::atan(tan_theta_lim);
//        if(abs_theta <= theta_lim and theta > pi / 2) {
//          geom_t const detlo =
//              std::sqrt(rlosq * one_plus_tsq - xsq * tsq) * one_on_one_plus_tsq;
//          geom_t const xlo = xterm1 + detlo;
//          geom_t const ylo = yterm1 + t * detlo;
//          d_lo = std::sqrt((xlo - x) * (xlo - x) + ylo * ylo);
//        }
//      }
//    }
//    // now select the shorter distance and the corresponding face
//    geom_t d_bdy = huge;
//    cell_t face = -1;
//    if(d_hi < d_lo) {
//      d_bdy = d_hi;
//      face = cell + 1;  // will intersect outer sphere
//    }
//    else {
//      d_bdy = d_lo;
//      face = cell;  // will intersect inner sphere
//    }
//    d_to_b_t d2b;
//    d2b.d = d_bdy;
//    d2b.face = face;
//    return d2b;
//  }  // dist_to_bdy_impl
//
//  static vector4_t LT_to_comoving(Vector1 const v_lab,
//                                  geom_t const & e_lab,
//                                  Vector1 const & omega_lab)
//  {
//    return spec_1D_Spherical::LT_to_comoving_sphere1D(v_lab, e_lab, omega_lab);
//  }  // LT_to_comoving
//
//  static vector4_t LT_to_lab(Vector1 const v_lab,
//                             geom_t const & e_com,
//                             Vector1 const & omega_com)
//  {
//    return spec_1D_Spherical::LT_to_lab_sphere1D(v_lab, e_com, omega_com);
//  }
//
//  vb const m_bdys;
//
//  vbd const m_descs;
//
//  cell_t const m_ncells;
//
//  void cellOK(cell_t const cell_idx) const
//  {
//    nut::InOpenRange(cell_idx, cell_t(0), m_ncells + 1, "cell id");
//  }
//
//};  // Sphere_1D

// using Spherical_1D_Mesh = Sphere_1D<cell_t, geom_t, bdy_types::descriptor>;

// template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
// const cell_t Sphere_1D<cell_t, geometry_t, bdy_descriptor_t>::null_cell_ =
//     0xFFFFFFFF;

// template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
// const typename Sphere_1D<cell_t, geometry_t, bdy_descriptor_t>::face_handle_t
//     Sphere_1D<cell_t, geometry_t, bdy_descriptor_t>::null_face_ =
//         static_cast<typename Sphere_1D<cell_t, geometry_t, bdy_descriptor_t>::
                        // face_handle_t>(0xFFFFFFFF);

}  // namespace nut

#endif

// End of file
