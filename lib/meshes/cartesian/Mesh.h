// mesh.h
// Jan 07, 2019
// (c) Copyright 2019 LANSLLC, all rights reserved

#pragma once

#include "mesh_common/Cell.h"
#include "mesh_common/Edge.h"
#include "mesh_common/Face.h"
#include "mesh_common/Ray.h"
#include "mesh_common/Vertex.h"
#include "mesh_common/types.h"
#include <map>
#include <unordered_set>
#include <vector>

namespace nut_mesh {

/**\class Mesh A fairly simple mesh class: it contains sets of mesh elements,
 * and the relations between them.
 */
template <class Face_T = nut_mesh::Face> class Mesh {
public:
  // types
  using cell_list_t = std::unordered_set<Cell>;
  using vertex_list_t = std::vector<Vertex>;
  using edge_list_t = std::vector<Edge>;
  using face_list_t = std::vector<Face_T>;

  // element to element maps
  using cell_to_vertex_map = std::map<Cell, vertex_list_t>;
  using vertex_to_cell_map = std::map<Vertex, cell_list_t>;
  using cell_to_face_map = std::map<Cell, face_list_t>;

public:
  Mesh() {}

  virtual ~Mesh() {}

  void add_cell(Cell &c, face_list_t &faces, vertex_list_t &vtcs) {
    cells_.insert(c);
    cell_to_vertex_[c] = vtcs;
    cell_to_face_[c] = faces;
  }

  /**\brief get the vertices connected to a cell*/
  vertex_list_t const &get_vertices(Cell const &c) {
    // Require(valid_cell(c));
    auto it = cell_to_vertex_.find(c);
    return (*it).second;
  } // get_vertices

  /**\brief Get the faces that bound a cell */
  face_list_t const &get_faces(Cell const &c) const {
    // Require(valid_cell(c));
    auto it = cell_to_face_.find(c);
    return (*it).second;
  } // get_faces

  cell_list_t const &cells() const { return cells_; }

  // virtual Cell & cell_across(Cell const &c, Face_T const &f) = 0;

private:
  bool valid_cell(Cell const &c) { return cells_.find(c) != cells_.end(); }

  cell_list_t cells_;

  // Not clear at this time whether we want to have a list of all vertices,
  // edges, or faces
  /*
  vertex_list_t vertices_;
  edge_list_t edges_;
  face_list_t faces_;
  */

  cell_to_vertex_map cell_to_vertex_;
  cell_to_face_map cell_to_face_;

  // Again, at some point we may want this?
  /*
  vertex_list_t vertices_;
  */

}; // class Mesh

} // namespace nut_mesh

// End of file
