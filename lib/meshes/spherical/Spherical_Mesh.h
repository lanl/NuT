// Spherical_Mesh.h
// T. M. Kelley
// Jun 25, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#include "mesh_common/Cell.h"
#include "mesh_common/Cell_Face_Descriptor.h"
#include "mesh_common/Face.h"
#include "mesh_common/Ray.h"
#include "mesh_common/Vector.h"
#include "mesh_common/boundary.h"
#include "mesh_common/types.h"
#include <map>
#include <vector>

namespace nut_mesh {

class Spherical_1D_Mesh {
  /* Implementation notes:
   * 1. We're relying on the cell boundaries being in a sorted vector.
   *
   */
public:
  // types
  using Ray = nut_mesh::Ray1;
  using Vector = nut_mesh::Vector1;
  using Point = nut_mesh::Vector1;
  using Cell = nut_mesh::Cell;
  using Face = Spherical_1D_Face;
  using Intersection_T = std::pair<Face, geom_t>;

  using Index = index_t;
  using Geom_T = geom_t;

  // convenient alias for the low- and high-extents of a cell in any direction
  using Extents_T = std::pair<geom_t, geom_t>;

  // private types
private:
  enum Face_Name { LOW = 0, HIGH = 1 };

  //
public:
  static const Cell null_cell_;

  static Cell null_cell() { return null_cell_; }

  static const Face null_face_;

  static Face null_face() { return null_face_; }

  /**\brief Interpreting Intersection_T */
  static geom_t get_distance(Intersection_T const & i);

  /**\brief Interpreting Intersection_T */
  static Face get_face(Intersection_T const & i);

  static geom_t get_low(Extents_T const & e);

  static geom_t get_high(Extents_T const & e);

public:
  // INTERFACE
  /**\brief Identify cell across the face at given point.
   *
   * \remark The returned cell might not be valid. For example, if you ask for
   * the cell across from the highest-r cell and face, you will get a cell
   * that is not part of the mesh. */
  Cell cell_across(Cell const & c, Face const & f, Point const & /*p*/) const;

  /**\brief Is this cell a boundary on any face? */
  bool is_boundary(Face f) const;

  /**brief Get {r_inner, r_outer} for Cell c. */
  Extents_T get_extents(Cell const & c) const;

  /**\brief Compute vector v reflected in face f */
  static Vector reflect(Face const & f, Vector const & v);

  /**\brief Get inward-facing normal for f */
  static Vector get_normal(Cell const & c, Face const & f);

  /**\brief Find the cell in which the Point lies.
   *
   * Note that this will always return *something*. It will return junk if you
   * provide a point that is not in the mesh. You might want to check the
   * returned cell value with in_cell() (though you may get consistent junk).
   * So you may really want to check valid_cell.
   */
  Cell find_cell(Point const & p) const;

  bool valid_cell(Cell const & c) const;

  /**\brief Is the point p in the cell c? */
  bool in_cell(Vector const & p, Cell const & c) const;

  /**\brief Compute intersection face and distance of ray r in Cell c.
   * \param r: current position & direction
   * \param c: current cell, must be consistent with r
   * \note The cell is a convenience, the position in the Ray must lie within
   * the Cell.
   */
  Intersection_T intersection(Ray const & r, Cell const & c) const;

  /**\brief Construct a face index from a cell and a face_name. Basically, if
   * it's the low face, the Face will have the same index as the Cell; if it's
   * the HIGH face, it will have index one greater than the cell. */
  Face cell_to_face(Cell const & c, Face_Name const fn) const;

  /**\brief The implementation of distance-to-boundary AKA intersection*/
  static Intersection_T dist_to_bdy_impl(geom_t x,
                                         geom_t omega,
                                         geom_t rlo,
                                         geom_t rhi,
                                         index_t cell_idx);

  Ray advance_point(Point const & init_pnt,
                    Vector const & dir,
                    double const dist) const;

  /**\brief Sample a random point in the Cell.
   * \tparam Any class that supports geom_t random(), returning a URD on [0,1)
   * \param: c: the cell in which to sample a point
   */
  template <typename RNG_T>
  Vector sample_position(RNG_T & r, Cell const & c) const;

  /**\brief Sample a random direction isotropically.
   * \tparam Any class that supports geom_t random(), returning a URD on [0,1)
   */
  template <typename RNG_T>
  static Vector sample_direction_isotropic(RNG_T & r);

  /**\brief Number of cells in this Cartesian mesh */
  index_t num_cells() const;

  /**\brief Is the face the low face in the cell?
   *
   * Note: the face need not be contiguous to the cell. */
  static bool is_low_face(Cell const & c, Face const & f);

  /**\brief Volume of indicated Cell */
  geom_t volume(Cell const cell) const;

  // Constructors
  explicit Spherical_1D_Mesh(index_t n_cells, geom_t dr = 1.0)
      : m_cell_bounds(n_cells + 1), m_num_cells(n_cells)
  {
    for(size_t i = 0; i < m_cell_bounds.size(); ++i) {
      m_cell_bounds[i] = static_cast<geom_t>(i) * dr;
    }
    return;
  }

  explicit Spherical_1D_Mesh(std::vector<geom_t> && bdys)
      : m_cell_bounds(std::move(bdys)), m_num_cells((m_cell_bounds.size() - 1))
  {
  }

  Spherical_1D_Mesh(Spherical_1D_Mesh const & other) = default;

  Spherical_1D_Mesh(Spherical_1D_Mesh && other) = default;

  Spherical_1D_Mesh & operator=(Spherical_1D_Mesh const & other) = default;

  Spherical_1D_Mesh & operator=(Spherical_1D_Mesh && other) = default;

private:
  // state
  std::vector<geom_t> m_cell_bounds;

  index_t m_num_cells;
};  // class Spherical_1D_Mesh

}  // namespace nut_mesh

#define IM_ALLOWED_TO_INCLUDE_Spherical_MESH_I_H
#include "Spherical_Mesh.i.h"
#undef IM_ALLOWED_TO_INCLUDE_Spherical_MESH_I_H

// End of file
