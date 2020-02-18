// Boundary.hh
// T. M. Kelley
// Oct 17, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#include "types.hh"
#include <unordered_map>

namespace nut {

/**\brief Store/retrieve boundary conditions */
template <class Face_Index_T, typename descriptor_t = bdy_types::descriptor>
struct Boundary_Cond {
  using face_t = Face_Index_T;
  using desc_t = descriptor_t;
  using bc_map_t = std::unordered_map<Face_Index_T, desc_t>;

  /**\brief Get the boundary type associated with Face f.
   *
   * If no boundary is associated with the Face f, the default boundary type
   * will be associated with f, so don't ask for this if you don't want
   * to waste time on useless faces. */
  desc_t get_boundary_type(face_t f) const
  {
    // Require(is_known_boundary(f));
    return m_bcs.at(f);
  }

  bool is_known_boundary(face_t f) const { return (m_bcs.count(f) != 0); }

  void set_boundary_type(face_t f, desc_t d) { m_bcs[f] = d; }

  size_t size() const { return m_bcs.size(); }

  bc_map_t m_bcs;
};

template <typename Mesh_1D>
Boundary_Cond<typename Mesh_1D::face_handle_t>
make_vacuum_boundary_1D(Mesh_1D const & mesh)
{
  using face_t = typename Mesh_1D::face_handle_t;
  Boundary_Cond<face_t> bcs;
  face_t low{0u};
  bcs.set_boundary_type(low, bdy_types::REFLECTIVE);
  face_t high{mesh.num_cells() + 1};
  bcs.set_boundary_type(high, bdy_types::VACUUM);
  return bcs;
}  // make_vacuum_boundary(Mesh_1D)

/* This specific function is going to make demands on Cartesian meshes;
 * specifically, that they implement and provide boundary face iterators. */
template <typename Mesh_Cartesian_3D>
Boundary_Cond<typename Mesh_Cartesian_3D::face_handle_t>
make_vacuum_boundary_3D(Mesh_Cartesian_3D const & mesh)
{
  using face_t = typename Mesh_Cartesian_3D::face_handle_t;
  // using real_mesh_t = typename Mesh_Cartesian_3D::mesh_t;
  using bdy_face_iterator = typename Mesh_Cartesian_3D::boundary_face_iterator;
  Boundary_Cond<face_t> bcs;
  // lambda to set all the faces from a boundary_face_iterator to VACUUM
  auto set_plane_bcs = [&bcs](bdy_face_iterator it) {
    while(it) {
      bcs.set_boundary_type(*it, bdy_types::VACUUM);
      it++;
    }
    return;
  };

  // auto const & real_mesh(mesh.get_mesh());
  auto i_l_x = mesh.boundary_faces_begin(Mesh_Cartesian_3D::LOW_X);
  auto i_l_y = mesh.boundary_faces_begin(Mesh_Cartesian_3D::LOW_Y);
  auto i_l_z = mesh.boundary_faces_begin(Mesh_Cartesian_3D::LOW_Z);
  auto i_h_x = mesh.boundary_faces_begin(Mesh_Cartesian_3D::HIGH_X);
  auto i_h_y = mesh.boundary_faces_begin(Mesh_Cartesian_3D::HIGH_Y);
  auto i_h_z = mesh.boundary_faces_begin(Mesh_Cartesian_3D::HIGH_Z);

  set_plane_bcs(i_l_x);
  set_plane_bcs(i_l_y);
  set_plane_bcs(i_l_z);
  set_plane_bcs(i_h_x);
  set_plane_bcs(i_h_y);
  set_plane_bcs(i_h_z);

  return bcs;
}  // make_vacuum_boundary(Mesh_3D)

}  // namespace nut

// End of file
