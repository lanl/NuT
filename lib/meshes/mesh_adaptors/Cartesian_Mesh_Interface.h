// Cartesian_Mesh_Interface.h
// Apr 08, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#include "cartesian/Cartesian_Mesh.h"
#include "geometry/lorentz.h"
#include "mesh/mesh_interface.h"
#include <array>

namespace murmeln {

class Cartesian_Mesh_Impl {
  // types
public:
  using mesh_t = murmeln_mesh::Cartesian_Mesh;
  using Ray = mesh_t::Ray;
  using cell_handle_t = mesh_t::cell_handle_t;
  using Vector = mesh_t::Vector;
  using Point = mesh_t::Vector;
  using Intersection = mesh_t::intersect_t;
  using Geom_State = std::pair<Ray, cell_handle_t>;
  using Index_T = mesh_t::Index;
  using Geom_T = mesh_t::Geom_T;
  using vector4_t = murmeln_mesh::Vector4<Geom_T>;

  // types for murmeln interface
  using index_t = Index_T;
  using area_t = Geom_T;
  // TO DO: make the different indices into distinct types
  using cell_index_t = index_t;
  using face_index_t = index_t;
  using vertex_index_t = index_t;
  using boundary_index_t = int;

  using global_vertex_index_t = vertex_index_t;
  using global_cell_index_t = cell_index_t;
  using global_face_index_t = face_index_t;
  using global_boundary_index_t = boundary_index_t;

  using bdy_index_t = index_t;
  using ticket_t = index_t;

  using volume_t = Geom_T;
  using distance_t = Geom_T;
  using vector_t = Vector;
  using point_t = Vector;
  using geometron_t = Geom_State;
  using vertex_handle_t = index_t;
  using intersect_t = Intersection;
  using face_handle_t = mesh_t::face_handle_t;
  using face_list_t = std::array<face_handle_t, 6>;
  using vertex_list_t = std::array<vertex_handle_t, 6>;

public:
  /**\brief Compute the closest face with which the ray will intersect, and
   * the distance to that intersection. */
  Intersection intersection(Ray const &r, cell_handle_t const &c) const {
    return mesh_.intersection(r, c);
  }

  Intersection intersection(Geom_State const &gs) const {
    auto const &r = gs.first;
    auto const &c = gs.second;
    return mesh_.intersection(r, c);
  }

  /**\brief Identify cell across the face at given point. */
  cell_handle_t cell_across(cell_handle_t const &c, face_handle_t const &f,
                            Point const &p) const {
    return mesh_.cell_across(c, f, p);
  }

  /**\brief Find the index of the cell in which the Point is located.*/
  cell_handle_t find_cell(Point const &p) const { return mesh_.find_cell(p); }

  bool in_cell(geometron_t const &g, cell_handle_t const &c) const {
    return mesh_.in_cell(g, c);
  }

  bool is_boundary(face_handle_t f) const { return mesh_.is_boundary(f); }

  Index_T num_cells() const { return mesh_.num_cells(); }

  Geom_T volume(cell_handle_t const &c) const { return mesh_.volume(c); }

  /**\brief Compute a point a given distance along a given direction */
  Ray advance_point(Point const &init_point, Vector const &direction,
                    double const distance) const {
    return mesh_.advance_point(init_point, direction, distance);
  }

  template <typename RNG_T>
  Vector sample_position(RNG_T &r, cell_handle_t const &c) const {
    return mesh_.sample_position(r, c);
  }

  template <typename RNG_T> static Vector sample_direction_isotropic(RNG_T &r) {
    return mesh_t::sample_direction_isotropic(r);
  }

  /**\brief Lorentz transform energy and direction from co-moving to lab
   * frame.
   * \param*/
  static vector4_t LT_to_lab(Vector const &v, Geom_T const &ec,
                             Vector const &oc) {
    return murmeln_mesh::spec_3D_Cartesian::LT_to_lab(v, ec, oc);
  }

  static vector4_t LT_to_comoving(Vector const &v, Geom_T const el,
                                  Vector const &ol) {
    return murmeln_mesh::spec_3D_Cartesian::LT_to_comoving(v, el, ol);
  }

  static cell_handle_t null_cell() { return mesh_t::null_cell(); }

  static face_handle_t null_face() { return mesh_t::null_face(); }

  Geom_T get_distance(Intersection const &i) const {
    return mesh_.get_distance(i);
  }

  face_handle_t get_face(Intersection const &i) const {
    return mesh_.get_face(i);
  }

  Vector reflect(Vector const &v, face_handle_t const &f) const {
    return mesh_.reflect(f, v);
  }

  Vector get_normal(cell_handle_t const &c, face_handle_t const &f) const {
    return mesh_.get_normal(c, f);
  }

  /**\brief Get number cells in X-dimension */
  index_t get_num_x() const { return mesh_.get_num_x(); }

  /**\brief Get number cells in Y-dimension */
  index_t get_num_y() const { return mesh_.get_num_y(); }

  /**\brief Get number cells in Z-dimension */
  index_t get_num_z() const { return mesh_.get_num_z(); }

  /** \brief Ctors, assignments */
  Cartesian_Mesh_Impl(mesh_t const &mesh) : mesh_(mesh) {}

  mesh_t const &get_mesh() const { return mesh_; }

  // Cartesian_Mesh_Impl(Cartesian_Mesh_Impl && rhs) noexcept = default;

  // Cartesian_Mesh_Impl(Cartesian_Mesh_Impl && rhs) noexcept = default;

  // ~Cartesian_Mesh_Impl() = default;

  // Cartesian_Mesh_Impl & operator=(Cartesian_Mesh_Impl const & rhs) = default;

  // Cartesian_Mesh_Impl & operator=(
  //   Cartesian_Mesh_Impl && rhs) noexcept = default;

private:
  // state
  mesh_t const &mesh_;
};

using Cartesian_Mesh_Interface = MeshInterface<Cartesian_Mesh_Impl>;

}; // namespace murmeln

// End of file
