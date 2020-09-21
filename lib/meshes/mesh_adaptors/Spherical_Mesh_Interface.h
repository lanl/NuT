// Spherical_Mesh_Interface.h
// Jul 09, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#include "geometry/lorentz.h"
#include "mesh/mesh_interface.h"
#include "spherical/Spherical_Mesh.h"
#include <array>

namespace nut {

class Spherical_Mesh_Impl {
  // types
public:
  using mesh_t = nut_mesh::Spherical_1D_Mesh;
  using Ray = mesh_t::Ray;
  using Cell = mesh_t::Cell;
  using Face = mesh_t::Face;
  using Vector = mesh_t::Vector;
  using Point = mesh_t::Vector;
  using Intersection = mesh_t::Intersection_T;
  using Geom_State = std::pair<Ray, Cell>;
  using Index_T = mesh_t::Index;
  using Geom_T = mesh_t::Geom_T;
  using EandOmega_T = nut_mesh::spec_1D_Spherical::EandOmega_T<Geom_T>;
  using vector4_t = EandOmega_T;

  // types for mesh interface
  using index_t = Index_T;
  // TO DO: make the different indices into distinct types
  using cell_index_t = index_t;
  using face_index_t = index_t;
  using vertex_index_t = index_t;
  using boundary_index_t = int;

  using global_cell_index_t = cell_index_t;
  using global_face_index_t = face_index_t;
  using global_vertex_index_t = vertex_index_t;
  using global_boundary_index_t = boundary_index_t;

  using bdy_index_t = index_t;
  using ticket_t = index_t;

  using area_t = Geom_T;
  using volume_t = Geom_T;
  using distance_t = Geom_T;
  using vector_t = Vector;
  using point_t = Vector;
  using geometron_t = Geom_State;
  using vertex_handle_t = index_t;
  using face_handle_t = Face;
  using cell_handle_t = Cell;
  using intersect_t = Intersection;

  using face_list_t = std::array<face_handle_t, 2>;
  using vertex_list_t = std::array<vertex_handle_t, 2>;

public:
  /**\brief Compute the closest face with which the ray will intersect, and the
   * distance to that intersection. */
  Intersection intersection(Ray const & r, Cell const & c) const
  {
    return mesh_.intersection(r, c);
  }

  Intersection intersection(Geom_State const & gs) const
  {
    auto const & r = gs.first;
    auto const & c = gs.second;
    return mesh_.intersection(r, c);
  }

  /**\brief Identify cell across the face at given point. */
  Cell cell_across(Cell const & c, Face const & f, Point const & p) const
  {
    return mesh_.cell_across(c, f, p);
  }

  /**\brief Find the index of the cell in which the Point is located.*/
  Cell find_cell(Point const & p) const { return mesh_.find_cell(p); }

  bool in_cell(Point const & p, Cell const & c) const
  {
    return mesh_.in_cell(p, c);
  }

  bool is_boundary(Face f) const { return mesh_.is_boundary(f); }

  Index_T num_cells() const { return mesh_.num_cells(); }

  Geom_T volume(Cell const & c) const { return mesh_.volume(c); }

  /**\brief Compute a point a given distance along a given direction */
  Ray advance_point(Point const & init_point,
                    Vector const & direction,
                    double const distance) const
  {
    return mesh_.advance_point(init_point, direction, distance);
  }

  template <typename RNG_T>
  Vector sample_position(RNG_T & r, Cell const & c) const
  {
    return mesh_.sample_position(r, c);
  }

  template <typename RNG_T>
  static Vector sample_direction_isotropic(RNG_T & r)
  {
    return mesh_t::sample_direction_isotropic(r);
  }

  /**\brief Lorentz transform energy and direction from co-moving to lab
   * frame.
   * \param*/
  static EandOmega_T LT_to_lab(Vector const & v,
                               Geom_T const & ec,
                               Vector const & oc)
  {
    return nut_mesh::spec_1D_Spherical::LT_to_lab(v, ec, oc);
  }

  static EandOmega_T LT_to_comoving(Vector const & v,
                                    Geom_T const el,
                                    Vector const & ol)
  {
    return nut_mesh::spec_1D_Spherical::LT_to_comoving(v, el, ol);
  }

  static Cell null_cell() { return mesh_t::null_cell(); }

  static Face null_face() { return mesh_t::null_face(); }

  Geom_T get_distance(Intersection const & i) const
  {
    return mesh_.get_distance(i);
  }

  Face get_face(Intersection const & i) const { return mesh_.get_face(i); }

  Vector reflect(Vector const & v, Face const & f) const
  {
    return mesh_.reflect(f, v);
  }

  Vector get_normal(Cell const & c, Face const & f) const
  {
    return mesh_.get_normal(c, f);
  }

  mesh_t const & get_mesh() const { return mesh_; }

  /** \brief Explicit ctors */
  explicit Spherical_Mesh_Impl(mesh_t const & mesh) : mesh_(mesh) {}

private:
  // state
  mesh_t const & mesh_;
};

using Spherical_Mesh_Interface = MeshInterface<Spherical_Mesh_Impl>;

};  // namespace nut

// End of file
